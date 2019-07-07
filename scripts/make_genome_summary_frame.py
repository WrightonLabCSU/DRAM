import pandas as pd
import re

from mag_annotator.summarize_genomes import FRAME_COLUMNS

pengfei_master_frame_loc = '/Users/shafferm/lab/AMG/genome_frames/Metabolism_AMG_v4_16June2019.xlsx'
cazy_frame_loc = '/Users/shafferm/lab/AMG/genome_frames/cazymes_for_new_annotation_lms.xlsx'
merops_frame_loc = '/Users/shafferm/lab/AMG/genome_frames/MEROPS_database_edited_MB.xlsx'
trna_frame = None

OUTPUT_FILE = 'genome_summary_frame.tsv'

"""
Make a frame for all data from all tables to act as frame for full genome summary output
Headers: gene_id, gene_description, module_id, module_description, module_family, key_gene
Source dbs: Pengfei's metabolism sheet, Lindsey's dbCAN/CAZY sheet, LB's MEROPs peptidases sheet,
            tRNAs from tRNA scan
For each db classifications are made to keep a uniform column set across dbs. Data from multiple
columns of the source sheets are concatenated to make them fit the common format.
"""

CAZY_TYPE_DICT = {'GH': 'Glycoside Hydrolases', 'GT': 'GlycosylTransferases', 'PL': 'Polysaccharide Lyases',
                  'CE': 'Carbohydrate Esterases', 'AA': 'Auxiliary Activities', 'CBM': 'Carbohydrate-Binding Modules'}


def parse_pengfei_frame(pengfei_master_frame):
    master_frame = pd.read_excel(pengfei_master_frame)
    new_frame = master_frame[['KO', 'gene', 'module', 'module/pathway_name', 'metabolism', 'key gene']]
    new_frame.columns = FRAME_COLUMNS
    new_frame['key_gene'] = [not pd.isna(i) for i in new_frame['key_gene']]
    return new_frame


def make_cazy_gene_description(acts_on, activity_summary, long_description):
    description_list = []
    if acts_on is not None:
        description_list.append(acts_on)
    if activity_summary is not None:
        description_list.append(activity_summary)
    if long_description is not None:
        description_list.append(long_description)
    return '; '.join(description_list)


def split_to_numbers_and_letters(str_):
    """
    Takes string that starts with letters and ends with numbers and splits into a
    letter string and a number string"""
    match = re.match(r"([a-z]+)([0-9]+)", str_, re.I)
    if match:
        items = match.groups()
        return items
    else:
        raise ValueError('could not split %s' % str_)

def parse_cazy_frame(cazy_frame_loc):
    cazy_frame = pd.read_excel(cazy_frame_loc)
    new_frame_dict = dict()
    new_frame_dict[FRAME_COLUMNS[0]] = cazy_frame['CAZyme/KO']

    # Get gene_description by merging Polysaccardies(s), activity summary and long description
    acts_on_list = ['Acts on %s' % i if not pd.isna(i) else None for i in cazy_frame['Polysaccharide(s)']]
    activity_summary_list = [i if not pd.isna(i) else None for i in cazy_frame['activity summary']]
    long_description_list = [i if not pd.isna(i) else None for i in cazy_frame['long description']]

    if not (len(acts_on_list) == len(activity_summary_list)) or not (len(acts_on_list) == len(long_description_list)):
        raise ValueError("Why aren't these all the same")

    new_frame_dict[FRAME_COLUMNS[1]] = [make_cazy_gene_description(acts_on_list[i],
                                                                   activity_summary_list[i],
                                                                   long_description_list[i])
                                        for i in range(len(acts_on_list))]

    cazy_family = [split_to_numbers_and_letters(i)[0] for i in cazy_frame['CAZyme/KO']]
    new_frame_dict[FRAME_COLUMNS[2]] = cazy_family
    new_frame_dict[FRAME_COLUMNS[3]] = [CAZY_TYPE_DICT[i] for i in cazy_family]
    new_frame_dict[FRAME_COLUMNS[4]] = cazy_frame['category']
    new_frame_dict[FRAME_COLUMNS[5]] = ['' for _ in range(cazy_frame.shape[0])]

    new_frame = pd.DataFrame.from_dict(new_frame_dict)
    new_frame = new_frame[FRAME_COLUMNS]
    return new_frame


def make_merops_gene_description(catlytic_type, general_specificity, additional_info):
    description_list = []
    if catlytic_type is not None:
        description_list.append(catlytic_type)
    if general_specificity is not None:
        description_list.append(general_specificity)
    if additional_info is not None:
        description_list.append(additional_info)
    return '; '.join(description_list)


def parse_merops_frame(merops_frame_loc):
    merops_frame = pd.read_excel(merops_frame_loc)
    new_frame_dict = dict()

    new_frame_dict[FRAME_COLUMNS[0]] = [i.Subfamily if i.Subfamily != 'No subfamily' else i.Family for _, i in merops_frame.iterrows()]

    # Make gene description by merging catlytic type, general specificity and additional info
    catalytic_type_list = ['Catlytic type: %s' % i if not pd.isna(i) else None for i in merops_frame['Catalytic Type']]
    general_specificity_list = [i if not pd.isna(i) else None for i in merops_frame['General Specificity']]
    additional_info_list = [i if not pd.isna(i) else None for i in merops_frame['Additional Info']]

    if not (len(catalytic_type_list) == len(general_specificity_list)) or not (len(catalytic_type_list) == len(additional_info_list)):
        raise ValueError("Why aren't these all the same")

    new_frame_dict[FRAME_COLUMNS[1]] = [make_merops_gene_description(catalytic_type_list[i],
                                                                     general_specificity_list[i],
                                                                     additional_info_list[i])
                                        for i in range(len(catalytic_type_list))]

    new_frame_dict[FRAME_COLUMNS[2]] = merops_frame['Family']
    new_frame_dict[FRAME_COLUMNS[3]] = merops_frame['Content of Family']
    new_frame_dict[FRAME_COLUMNS[4]] = ['Peptidases' for _ in range(merops_frame.shape[0])]
    new_frame_dict[FRAME_COLUMNS[5]] = ['' for _ in range(merops_frame.shape[0])]

    new_frame = pd.DataFrame.from_dict(new_frame_dict)
    new_frame = new_frame[FRAME_COLUMNS]
    return new_frame


def main(master_frame_loc, cazy_frame_loc, merops_frame_loc, output):
    master_frame = parse_pengfei_frame(master_frame_loc)
    cazy_frame = parse_cazy_frame(cazy_frame_loc)
    merops_frame = parse_merops_frame(merops_frame_loc)
    #TODO: Add tRNAs

    genome_summary_frame = pd.concat([master_frame, cazy_frame, merops_frame])
    genome_summary_frame.to_csv(output, sep='\t', index=None)


if __name__ == '__main__':
    main(pengfei_master_frame_loc, cazy_frame_loc, merops_frame_loc, output=OUTPUT_FILE)
