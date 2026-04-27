#!/usr/bin/env python
"""This is the script that distills the genomes"""
import click
import os
from pathlib import Path
import polars as pl
from xlsxwriter import Workbook
from utils.logger import get_logger
from utils.click_utils import validate_comma_separated
from utils.click_utils import validate_comma_separated
from utils.excel import write_summarized_genomes_to_xlsx
from rule_parser.src.rules import evaluate_rules_on_anno, ID_EXPR_DICT

logger = get_logger(filename=Path(__file__).stem)

COL_GENE_ID, COL_GENE_DESCRIPTION, COL_MODULE, COL_SHEET, COL_HEADER, COL_SUBHEADER, RULES_PARENT, RULES = 'gene_id', 'gene_description', 'pathway', 'topic_ecosystem','category', 'subcategory', 'parent', 'rules'
OPTIONAL_COLUMNS = [RULES_PARENT, RULES]
RRNA_COLUMNS = [COL_GENE_ID, COL_GENE_DESCRIPTION, COL_SHEET, COL_HEADER, COL_SUBHEADER]
TRNA_COLUMNS = RRNA_COLUMNS + ['AA_type']
CORE_COLUMNS = RRNA_COLUMNS + [COL_MODULE]
FRAME_COLUMNS = CORE_COLUMNS + OPTIONAL_COLUMNS
RRNA_TYPES = ['5S rRNA', '16S rRNA', '23S rRNA']
TAXONOMY_LEVELS = ['d', 'p', 'c', 'o', 'f', 'g', 's']
CONSTANT_DISTILLATE_COLUMNS = [COL_GENE_ID, COL_GENE_DESCRIPTION, COL_MODULE, COL_HEADER, COL_SUBHEADER]
DISTILATE_SORT_ORDER_COLUMNS = [COL_HEADER, COL_SUBHEADER, COL_MODULE, COL_GENE_ID]
EXCEL_MAX_CELL_SIZE = 32767

DISTILL_DIR = Path(__file__).parent / "assets/forms/distill_sheets"

    
def check_columns(data, logger):
    functions = [i for i in ID_EXPR_DICT if i in data.columns]
    missing = [i for i in ID_EXPR_DICT if i not in data.columns]
    logger.info("Note: the following id fields "
          f"were not in the annotations file and are not being used: {missing},"
          f" but these are {functions}")


def get_ids_from_annotations_by_row(data):
    functions = {i:j for i,j in FUNCTION_DICT.items() if i in data.columns}
    out = data.apply(lambda x: {i for k, v in functions.items() if not pd.isna(x[k])
                          for i in v(str(x[k])) if not pd.isna(i)}, axis=1)
    return out


def get_ids_from_annotations_all(data):
    data =  get_ids_from_annotations_by_row(data)
    data.apply(list)
    out = Counter(chain(*data.values))
    return out


def fill_genome_summary_frame(annotations, genome_summary_frame, groupby_column, logger):
    genome_summary_id_sets = [set([str(k).strip() for k in j.split(',')]) for j in genome_summary_frame[COL_GENE_ID]]
    logger.info(f"Genome summary ID sets: {genome_summary_id_sets}")
    
    def fill_a_frame(frame: pd.DataFrame):
        id_dict = get_ids_from_annotations_all(frame)
        logger.info(f"ID dictionary for {frame.name}: {id_dict}")
        
        counts = list()
        for set_ in genome_summary_id_sets:
            identifier_count = 0
            for gene_id in set_:
                # Try matching with and without '.hmm'
                matching_keys = [key for key in id_dict.keys() if gene_id == key or key.startswith(gene_id + ".")]
                for key in matching_keys:
                    identifier_count += id_dict[key]
            counts.append(identifier_count)
        # logger.info(f"Counts for {frame.name}: {counts}")
        
        return pd.Series(counts, index=genome_summary_frame.index)
    
    counts = annotations.groupby(groupby_column, sort=False)[annotations.columns].apply(fill_a_frame)
    genome_summary_frame = pd.concat([genome_summary_frame, counts.T], axis=1)
    
    return genome_summary_frame


def fill_genome_summary_frame_gene_names(annotations, genome_summary_frame, groupby_column, logger):
    genome_summary_id_sets = [set([k.strip() for k in j.split(',')]) for j in genome_summary_frame[COL_GENE_ID]]
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        # make dict of identifiers to gene names
        id_gene_dict = defaultdict(list)
        for gene, ids in get_ids_from_annotations_by_row(frame).items():
            for id_ in ids:
                id_gene_dict[id_].append(gene)
        # fill in genome summary_frame
        values = list()
        for id_set in genome_summary_id_sets:
            this_value = list()
            for id_ in id_set:
                this_value += id_gene_dict[id_]
            values.append(','.join(this_value))
        genome_summary_frame[genome] = values
    return genome_summary_frame


def summarize_rrnas(rrnas_df, groupby_column="input_fasta"):
    genome_rrna_dict = dict()
    for genome, frame in rrnas_df.groupby(groupby_column):
        genome_rrna_dict[genome] = Counter(frame['type'])
    row_list = list()
    for rna_type in RRNA_TYPES:
        row = [rna_type, '%s ribosomal RNA gene' % rna_type.split()[0], 'rRNA', 'rRNA', '', '']
        for genome, rrna_dict in genome_rrna_dict.items():
            row.append(genome_rrna_dict[genome].get(rna_type, 0))
        row_list.append(row)
    rrna_frame = pd.DataFrame(row_list, columns=FRAME_COLUMNS + list(genome_rrna_dict.keys()))
    return rrna_frame


def make_genome_summary(annotations, genome_summary_frame, logger, groupby_column="input_fasta"):
    summary = genome_summary_frame.collect()

    if RULES not in summary.columns:
        summary = summary.with_columns(pl.lit(None, dtype=pl.Utf8).alias(RULES))
    summary = summary.with_columns(
        pl.when(pl.col(RULES).is_null() | (pl.col(RULES).cast(pl.Utf8) == ""))
          .then(pl.col(COL_GENE_ID))
          .otherwise(pl.col(RULES))
          .alias(RULES)
    )

    rule_hits = evaluate_rules_on_anno(
        annotations=annotations,
        sample_col="query_id",
        rules=summary.lazy(),
        label_col=COL_GENE_ID,
        parent_col=None,
        rules_col=RULES,
    )

    counts = (
        rule_hits.join(
            annotations.select([pl.col("query_id"), pl.col(groupby_column)]),
            on="query_id",
        )
        .drop("query_id")
        .group_by(groupby_column)
        .agg(pl.exclude(groupby_column).sum())
    )

    counts = counts.select(pl.exclude(groupby_column)).transpose(
        include_header=True,
        header_name=COL_GENE_ID,
        column_names=counts[groupby_column],
    )

    return summary.drop(RULES).join(counts, on=COL_GENE_ID, how="left")


def split_column_str(names):
    if len(names) < EXCEL_MAX_CELL_SIZE:
        return [names]
    out = ['']
    name_list = names.split(',')
    j = 0
    for i in name_list:
        if len(out[j]) + len(i) + 1 < EXCEL_MAX_CELL_SIZE:
            out[j] = ','.join([out[j], i])
        else:
            j += 1
            out += ['']
    return out


# TODO: add assembly stats like N50, longest contig, total assembled length etc
def make_genome_stats(annotations, rrna_frame=None, trna_frame=None, quast_frame=None, groupby_column="input_fasta"):
    rows = list()
    columns = ['genome']
    if 'scaffold' in annotations.columns:
        columns.append('number of scaffolds')
    if 'bin_taxonomy' in annotations.columns:
        columns.append('taxonomy')
    if 'bin_completeness' in annotations.columns:
        columns.append('completeness score')
    if 'bin_contamination' in annotations.columns:
        columns.append('contamination score')
    for genome, frame in annotations.group_by(groupby_column):
        row = [genome[0]]
        if 'scaffold' in frame.columns:
            row.append(len(set(frame['scaffold'])))
        if 'bin_taxonomy' in frame.columns:
            row.append(frame['bin_taxonomy'][0])
        if 'bin_completeness' in frame.columns:
            row.append(frame['bin_completeness'][0])
        if 'bin_contamination' in frame.columns:
            row.append(frame['bin_contamination'][0])
        rows.append(row)
    genome_stats = pl.DataFrame(rows, schema=columns, orient='row')
    if rrna_frame is not None:
        meta_cols = RRNA_COLUMNS
        sample_cols = [c for c in rrna_frame.columns if c not in meta_cols]

        df_rrna = (
            rrna_frame
            .group_by("gene_id")
            .agg([pl.col(c).sum().alias(c) for c in sample_cols])
            .transpose(include_header=True, header_name="genome", column_names="gene_id")
        )
        genome_stats = genome_stats.join(df_rrna, on="genome", how="left")
    if trna_frame is not None:
        meta_cols = TRNA_COLUMNS

        sample_cols = [c for c in trna_frame.columns if c not in meta_cols]

        df_trna = (
            trna_frame
            .filter(~pl.col("AA_type").is_in(["Undet", "Sup"]))
            .group_by("AA_type")
            .agg([pl.col(c).sum().alias(c) for c in sample_cols])
            .select([(pl.col(c) != 0).cast(pl.Int64).sum().alias(c) for c in sample_cols])
            .transpose(include_header=True, header_name="genome", column_names=["tRNA count"])
        )
        genome_stats = genome_stats.join(df_trna, on="genome", how="inner")
        
    if quast_frame is not None:
        quast_frame = (
            quast_frame
            .rename({groupby_column: "genome"})
            .drop("no. contigs")
        )

        genome_stats = genome_stats.join(quast_frame, on="genome", how="inner")
        assert genome_stats.shape[0] == quast_frame.shape[0], "genomes from annotation file don't map to quast file"

    return genome_stats


@click.command()
@click.option("-i", "--input_file", required=True, help="Annotations path")
# @click.option("-o", "--output_dir", required=True, help="Directory to write summarized genomes")
@click.option("--rrna_path", help="rRNA output from annotation", default=None, type=click.Path(exists=True))
@click.option("--trna_path", help="tRNA output from annotation", default=None, type=click.Path(exists=True))
@click.option("--quast_path", help="Quast summary TSV from the quast step", default=None, type=click.Path(exists=True))
@click.option("--groupby_column", help="Column from annotations to group as organism units",
                            default="input_fasta", type = click.STRING)
@click.option("--distil_topics", default="default", help="Default distillates topics to run.")
@click.option("--distil_ecosystem", default="eng_sys,ag", help="Default distillates ecosystems to run.")
@click.option("--custom_distillate", default="", callback=validate_comma_separated, help="Custom distillate forms to add your own modules, comma separated. ")
@click.option("--distillate_gene_names", is_flag=True,
    show_default=True, default=False,
                            help="Give names of genes instead of counts in genome metabolism summary")
def distill(input_file, rrna_path, trna_path, quast_path, groupby_column, distil_topics, distil_ecosystem,
                      custom_distillate, distillate_gene_names):
    """Summarize metabolic content of annotated genomes"""

    # read in data
    try:
        annotations = pl.read_csv(input_file, separator="\t", infer_schema_length=10_000)
    except Exception as e:
        annotations = pl.read_csv(input_file, separator="\t", infer_schema_length=None)
    if 'bin_taxonomy' in annotations.columns:
        annotations = annotations.sort('bin_taxonomy')

    # Check the columns are present
    check_columns(annotations, logger)

    trna_frame = None
    rrna_frame = None
    if trna_path is not None and rrna_path is not None:
        try:
            trna_frame = pl.read_csv(trna_path, separator='\t')
            rrna_frame = pl.read_csv(rrna_path, separator='\t')
        except Exception:
            trna_frame = rrna_frame = None
        if trna_frame is None or rrna_frame is None or trna_frame.is_empty() or rrna_frame.is_empty():
            trna_frame = rrna_frame = None

    quast_frame = None
    if quast_path is not None:
        try:
            quast_frame = pl.read_csv(quast_path, separator='\t')
        except Exception:
            quast_frame = None
        if quast_frame is not None and quast_frame.is_empty():
            quast_frame = None

    distil_sheets_names = []
    if "default" in distil_topics:
        distil_sheets_names = [
            DISTILL_DIR / "distill_carbon.tsv",
            DISTILL_DIR / "distill_energy.tsv",
            DISTILL_DIR / "distill_misc.tsv",
            DISTILL_DIR / "distill_nitrogen.tsv",
            DISTILL_DIR / "distill_transport.tsv",
            DISTILL_DIR / "distill_metals.tsv"
        ]
    else:
        if 'carbon' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_carbon.tsv")
        if 'energy' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_energy.tsv")
        if 'misc' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_misc.tsv")
        if 'nitrogen' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_nitrogen.tsv")
        if 'transport' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_transport.tsv")
        if "metals" in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_metals.tsv")
    
        
    if "ag" in distil_ecosystem:
        distil_sheets_names.append(DISTILL_DIR / "distill_ag.tsv")
    if "eng_sys" in distil_ecosystem:
        distil_sheets_names.append(DISTILL_DIR / "distill_eng_sys.tsv")
    
    if "camper_id" in annotations and ("default" in distil_topics or "camper" in distil_topics):
        distil_sheets_names.append(DISTILL_DIR / "distill_camper.tsv")
        
    logger.info(f"Distillate dir: {DISTILL_DIR}")
    logger.info(f"Distillate sheets to be used: {distil_sheets_names}")
    if custom_distillate:
        for custom_sheet in custom_distillate:
            distil_sheets_names.append(custom_sheet)
    
    genome_summary_form = pl.concat(
        [
            pl.scan_csv(s, separator="\t")
            .select([c for c in FRAME_COLUMNS if c in pl.scan_csv(s, separator="\t", n_rows=0).columns])
            for s in distil_sheets_names
        ],
        how="diagonal",
    )
    
    logger.info('Retrieved distillate genome summary form')

    # make genome stats
    genome_stats = make_genome_stats(annotations, rrna_frame, trna_frame, quast_frame, groupby_column=groupby_column)
    genome_stats.write_csv('genome_stats.tsv', separator='\t')
    logger.info('Calculated genome statistics')

    # make genome metabolism summary
    genome_summary = 'metabolism_summary.xlsx'
    logger.info(f'Giving counts for genome metabolism summary')
    summarized_genomes = make_genome_summary(annotations, genome_summary_form, logger, groupby_column)
    summarized_genomes.write_csv('summarized_genomes.tsv', separator='\t')
    kw = {"extra_frames": []}
    if rrna_frame is not None:
        kw["extra_frames"].append(rrna_frame)
    if trna_frame is not None:
        kw["extra_frames"].append(trna_frame)
    write_summarized_genomes_to_xlsx(
        df=summarized_genomes,
        output_file=genome_summary,
        group_by=COL_SHEET,
        sort_order_columns=DISTILATE_SORT_ORDER_COLUMNS,
        **kw
    )
    logger.info('Generated genome metabolism summary')

    
if __name__ == "__main__":
    distill()
