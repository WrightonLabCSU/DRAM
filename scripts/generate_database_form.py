import argparse
import pandas as pd
import re
from itertools import tee
from KEGG_parser import downloader, parsers
from collections import defaultdict
from functools import partial


def read_modules(modules_loc):
    return [module.strip() for module in open(modules_loc).readlines()]


def get_value_from_complex_key(ko: str, orthology: dict):
    for kos, name in orthology.items():
        if ko in kos:
            return name
    raise KeyError(ko, orthology.keys())


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def first_open_paren_is_all(str_):
    curr_level = 1
    for i, char in enumerate(str_[1:-1]):
        if char == ')':
            curr_level -= 1
        elif char == '(':
            curr_level += 1
        if curr_level == 0:
            return False
    return True


def split_into_steps(definition, split_char=' '):
    """Very fancy split on space"""
    curr_level = 0
    step_starts = [-1]
    for i, char in enumerate(definition):
        if char is '(':
            curr_level += 1
        if char is ')':
            curr_level -= 1
        if (curr_level == 0) and (char in split_char):
            step_starts.append(i)
    step_starts.append(len(definition))
    steps = list()
    for a, b in pairwise(step_starts):
        step = definition[a+1:b]
        if step.startswith('(') and step.endswith(')'):
            if first_open_paren_is_all(step):
                step = step[1:-1]
        steps.append(step)
    return steps


split_into_steps_space = partial(split_into_steps, split_char=' ')
split_into_steps_comma = partial(split_into_steps, split_char=',')


def parse_steps(definition, splitter=split_into_steps_space, flipper=False, level=''):
    steps = splitter(definition)
    ko_dict_list = list()
    if len(steps) > 1:
        for i, step in enumerate(steps):
            if flipper:
                ko_dict_list += parse_steps(step, split_into_steps_space, False, level='%s,%s' % (level, i))
            else:
                ko_dict_list += parse_steps(step, split_into_steps_comma, True, level='%s,%s' % (level, i))
    else:
        level = level[1:]
        parsed_step = steps[0]
        for ko in split_into_steps(parsed_step, '+-'):
            if ',' in ko:
                parsed_substep = split_into_steps(ko, ',')
                ko_loc = parsed_step.find(ko)
                if ko_loc == 0:
                    required = True
                elif ko_loc == -1:
                    required = True
                elif parsed_step[ko_loc - 2] == '+':
                    required = True
                elif parsed_step[ko_loc - 2] == '-':
                    required = False
                else:
                    raise ValueError('%s from %s in %s at %s not expected' % (parsed_step[ko_loc - 2], ko,
                                                                              parsed_step, ko_loc))
                for i, substep in enumerate(parsed_substep):
                    ko_dict_list += [{'ko': substep, 'step': '%s(%s)' % (level, i), 'required': required}]
            else:
                ko_loc = parsed_step.find(ko)
                if ko_loc == 0:
                    required = True
                elif ko_loc == -1:
                    required = True
                elif parsed_step[ko_loc-1] == '+':
                    required = True
                elif parsed_step[ko_loc-1] == '-':
                    required = False
                else:
                    raise ValueError('%s from %s in %s at %s not expected' % (parsed_step[ko_loc-1], ko, parsed_step,
                                                                              ko_loc))
                ko_dict_list += [{'ko': ko, 'step': level, 'required': required}]

    return ko_dict_list


def get_rows(modules_dict):
    """substrate, product, substrate_id, product_id, module, gene"""
    all_ko_dicts = list()
    for module, module_dict in modules_dict.items():
        print(module)
        if ' ' in module_dict['DEFINITION']:
            kos_dict = parse_steps(module_dict['DEFINITION'])  # gets step and required
        else:
            kos_dict = parse_steps(module_dict['DEFINITION'], split_into_steps_comma, True, '0,')
        for ko_dict in kos_dict:
            ko_dict['module'] = module
            ko = ko_dict['ko']
            protein_name = get_value_from_complex_key(ko, module_dict['ORTHOLOGY'])
            ko_dict['gene'] = protein_name
            if 'REACTION' in module_dict:
                rn_name_loc = protein_name.find('RN:')
                if rn_name_loc != -1:
                    reaction = protein_name[rn_name_loc + 3:rn_name_loc + 9]
                    reaction_entry = get_value_from_complex_key(reaction, module_dict['REACTION'])
                    substrate_ids = reaction_entry[0]
                    ko_dict['substrate_ids'] = ','.join(substrate_ids)
                    ko_dict['substrate_names'] = ','.join([module_dict['COMPOUND'][substrate_id]
                                                           for substrate_id in substrate_ids])
                    product_ids = reaction_entry[1]
                    ko_dict['product_ids'] = ','.join(product_ids)
                    ko_dict['product_names'] = ','.join([module_dict['COMPOUND'][product_id]
                                                         for product_id in product_ids])
        all_ko_dicts += kos_dict
    return all_ko_dicts


def main(modules_list: list=None, modules_loc=None, output='empty_form.tsv'):
    if modules_list is None:
        modules_list = list()
    if modules_loc is not None:
        modules_list += read_modules(modules_loc)
    modules_dict = downloader.get_kegg_record_dict(modules_list, parsers.parse_module)
    row_dict = get_rows(modules_dict)
    form = pd.DataFrame.from_dict(row_dict)
    form.to_csv(output, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--module_list', nargs='*', help='Space separated list of KEGG Module IDs')
    parser.add_argument('--modules_loc', help='File with one KEGG module per line')
    parser.add_argument('--output', help='Name of output tsv', default='form.tsv')

    args = parser.parse_args()
    main(args.module_list, args.modules_loc, args.output)
