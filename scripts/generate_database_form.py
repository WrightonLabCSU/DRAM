import argparse
import pandas as pd
import re
from itertools import tee
from KEGG_parser import downloader, parsers


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
    """Very fancy split on string of chars"""
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


def parse_definition(definition, path_tuple_base=()):
    steps = split_into_steps(definition, split_char=' ')
    ko_paths = list()
    for i, step in enumerate(steps):
        for j, path in enumerate(split_into_steps(step, ',')):
            path_tuple = path_tuple_base+(i,j)
            if ' ' in path:
                ko_paths += parse_definition(path, path_tuple)
            else:
                ko_paths.append((path, ','.join(map(str, path_tuple))))
    return ko_paths


def parse_complex(complex):
    # TODO: get if required or not and if member of optional subcomplex, optional subcomplexes, 'or' subcomplexs, more
    return re.findall('K\d\d\d\d\d', complex)


def get_rows(modules_dict):
    """substrate, product, substrate_id, product_id, module, gene"""
    all_ko_dicts = list()
    for module, module_dict in modules_dict.items():
        print(module)
        complex_list = parse_definition(module_dict['DEFINITION'])
        kos_dict = []
        for complex, path in complex_list:
            for ko in parse_complex(complex):
                kos_dict.append({'ko': ko, 'path': path})
        for ko_dict in kos_dict:
            ko_dict['module'] = module
            ko_dict['module_name'] = module_dict['NAME']
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
    form.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--module_list', nargs='*', help='Space separated list of KEGG Module IDs')
    parser.add_argument('--modules_loc', help='File with one KEGG module per line')
    parser.add_argument('--output', help='Name of output tsv', default='form.tsv')

    args = parser.parse_args()
    main(args.module_list, args.modules_loc, args.output)
