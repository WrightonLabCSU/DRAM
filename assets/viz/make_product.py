"""
===================
DRAM Visualizations
===================

Script that generates a product visualization from the DRAM output.
"""
import logging
import re
from collections import Counter, defaultdict
from itertools import tee, chain
from math import pi
from pathlib import Path
from typing import Optional

import networkx as nx
import numpy as np
import pandas as pd
import panel as pn
from bokeh.models import ColorBar, LinearColorMapper
from bokeh.palettes import BuGn, Cividis256
from bokeh.plotting import figure
from bokeh.resources import INLINE as INLINE_RESOURCES
from bokeh.transform import factor_cmap, linear_cmap

__authors__ = ["Madeline Scyphers"]
__copyright__ = "Copyright 2024, Wrighton Lab"
__license__ = "NA"

MIN_DB_SET: list[set[str]] = [{"kegg"}, {"kofam"}]

logger = logging.getLogger("dram2_log.viz")

# from dram2.annotate import  get_ids_from_annotations_all

__version__ = "3.0.0"
DRAM_DATAFOLDER_TAG = "dram_data_folder"
DBSETS_COL = "db_id_sets"
DEFAULT_GROUPBY_COLUMN = "fasta"
DEFAULT_GENOMES_PER_PRODUCT = 1000
GENOMES_PRODUCT_LIMIT = 2000
DEFAULT_MAKE_BIG_HTML = False
DRAM_SHEET_TAG = "dram_sheets"
GENOME_SUMMARY_FORM_TAG = "genome_summary_form"
MODULE_STEPS_FORM_TAG = "module_step_form"
FUNCTION_HEATMAP_FORM_TAG = "function_heatmap_form"
ETC_MODULE_DF_TAG = "etc_module_database"

FILES_NAMES: dict[str, Path] = {
    GENOME_SUMMARY_FORM_TAG: Path(__file__).parent.resolve() / "data/genome_summary_form.tsv",
    MODULE_STEPS_FORM_TAG: Path(__file__).parent.resolve() / "data/module_step_form.tsv",
    FUNCTION_HEATMAP_FORM_TAG: Path(__file__).parent.resolve() / "data/function_heatmap_form.tsv",
    ETC_MODULE_DF_TAG: Path(__file__).parent.resolve() / "data/etc_module_database.tsv",
}

ID_FUNCTION_DICT = {
    'kegg_genes_id': lambda x: [x],
    'ko_id': lambda x: [j for j in x.split(',')],
    'kegg_id': lambda x: [j for j in x.split(',')],
    'kegg_hit': lambda x: [i[1:-1] for i in
                           re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', x)],
    'peptidase_family': lambda x: [j for j in x.split(';')],
    'cazy_best_hit': lambda x: [x.split('_')[0]],
    'pfam_hits': lambda x: [j[1:-1].split('.')[0]
                            for j in re.findall(r'\[PF\d\d\d\d\d.\d*\]', x)],
    'camper_id': lambda x: [x],
    'fegenie_id': lambda x: [x],
    'sulfur_id': lambda x: [x],
    'methyl_id': lambda x: [i.split(' ')[0].strip() for i in x.split(',')]
}

# TODO: add RBH information to output
# TODO: add flag to output table and not xlsx


PALETTE_CATEGORICAL = BuGn
PALETTE_CONTINUOUS = Cividis256

HEATMAP_MODULES = [
    "M00001",
    "M00004",
    "M00008",
    "M00009",
    "M00012",
    "M00165",
    "M00173",
    "M00374",
    "M00375",
    "M00376",
    "M00377",
    "M00422",
    "M00567",
]
HEATMAP_CELL_HEIGHT = 15
HEATMAP_CELL_WIDTH = 15
KO_REGEX = r"^K\d\d\d\d\d$"
ETC_COVERAGE_COLUMNS = [
    "module_id",
    "module_name",
    "complex",
    "genome",
    "path_length",
    "path_length_coverage",
    "percent_coverage",
    "genes",
    "missing_genes",
    "complex_module_name",
]
TAXONOMY_LEVELS = ["d", "p", "c", "o", "f", "g", "s"]
LOCATION_TAG = "location"


def get_distillate_sheet(form_tag: str, dram_config: dict, logger: logging.Logger):
    """
    Paths in the config can be complicated. Here is a function that will get
    you the absolute path, the relative path, or whatever. Specifically for
    distillate sheets This should be more
    formalized and the config

    should actually be managed in its own structure. With a data file class
    that can use this function.
    """
    if (
            (dram_sheets := dram_config.get(DRAM_SHEET_TAG)) is None
            or dram_sheets.get(form_tag) is None
            or (sheet_path_str := dram_sheets[form_tag].get(LOCATION_TAG)) is None
    ):
        sheet_path: Path = FILES_NAMES[form_tag]
        logger.debug(
            f"""
            Using the default distillation sheet for {form_tag} with its location at
            {sheet_path}. This information is only important if you intended to use a
            custom distillation sheet.
            """
        )
    else:
        sheet_path = Path(sheet_path_str)
        if not sheet_path.is_absolute():
            dram_data_folder: Optional[str] = dram_config.get(DRAM_DATAFOLDER_TAG)
            if dram_data_folder is None:
                raise DramUsageError(
                    f"""
                    In the DRAM2 config File, the path {form_tag} is a relative path
                    and the {DRAM_DATAFOLDER_TAG} is not set!
                    """
                )
            sheet_path = Path(dram_data_folder) / sheet_path
        if not sheet_path.exists():
            raise DramUsageError(
                f"""
                The file {form_tag} is not at the path {sheet_path}. Most likely you
                moved the DRAM data but forgot to update the config file to point to
                it. The easy fix is to set the {DRAM_DATAFOLDER_TAG} variable in the
                config like:\n
                {DRAM_DATAFOLDER_TAG}: the/path/to/my/file Iyou are using full paths
                and not the {DRAM_DATAFOLDER_TAG} you may want to revue the Configure
                Dram section othe documentation to make sure your config will work with
                dram. rememberer that the config must be a valid yaml file to work.
                Also you can always use db_builder to remake your databases and the
                config file iyou don't feel up to editing it yourself.
                """
            )

    return pd.read_csv(sheet_path, sep="\t")


def build_module_net(module_df):
    """Starts with a data from including a single module"""
    # build net from a set of module paths
    num_steps = max([int(i.split(",")[0]) for i in set(module_df["path"])])
    module_net = nx.DiGraph(
        num_steps=num_steps,
        module_id=list(module_df["module"])[0],
        module_name=list(module_df["module_name"])[0],
    )
    # go through all path/step combinations
    for module_path, frame in module_df.groupby("path"):
        split_path = [int(i) for i in module_path.split(",")]
        step = split_path[0]
        module_net.add_node(module_path, kos=set(frame["ko"]))
        # add incoming edge
        if step != 0:
            module_net.add_edge("end_step_%s" % (step - 1), module_path)
        # add outgoing edge
        module_net.add_edge(module_path, "end_step_%s" % step)
    return module_net


def get_module_step_coverage(kos, module_net):
    # prune network based on what kos were observed
    pruned_module_net = module_net.copy()
    module_kos_present = set()
    for node, data in module_net.nodes.items():
        if "kos" in data:
            ko_overlap = data["kos"] & kos
            if len(ko_overlap) == 0:
                pruned_module_net.remove_node(node)
            else:
                module_kos_present = module_kos_present | ko_overlap
    # count number of missing steps
    missing_steps = 0
    for node, data in pruned_module_net.nodes.items():
        if ("end_step" in node) and (pruned_module_net.in_degree(node) == 0):
            missing_steps += 1
    # get statistics
    num_steps = pruned_module_net.graph["num_steps"] + 1
    num_steps_present = num_steps - missing_steps
    coverage = num_steps_present / num_steps
    return num_steps, num_steps_present, coverage, sorted(module_kos_present)


def make_module_coverage_df(annotation_df, module_nets):
    kos_to_genes = defaultdict(list)
    ko_id: Optional[str] = None
    ko_id_names: list[str] = ["kegg_id", "kofam_id", "ko_id"]
    for id in ko_id_names:
        if id in annotation_df:
            ko_id = id
            break
    if ko_id is None:
        raise ValueError(
            f"""
            No KEGG or KOfam id column could be found.
            These names were tried: {', '.join(ko_id_names)}.
            """
        )
    for gene_id, ko_list in annotation_df[ko_id].items():
        if type(ko_list) is str:
            for ko in ko_list.split(","):
                kos_to_genes[ko].append(gene_id)
    coverage_dict = {}
    for i, (module, net) in enumerate(module_nets.items()):
        (
            module_steps,
            module_steps_present,
            module_coverage,
            module_kos,
        ) = get_module_step_coverage(set(kos_to_genes.keys()), net)
        module_genes = sorted([gene for ko in module_kos for gene in kos_to_genes[ko]])
        coverage_dict[module] = [
            net.graph["module_name"],
            module_steps,
            module_steps_present,
            module_coverage,
            len(module_kos),
            ",".join(module_kos),
            ",".join(module_genes),
        ]
    coverage_df = pd.DataFrame.from_dict(
        coverage_dict,
        orient="index",
        columns=[
            "module_name",
            "steps",
            "steps_present",
            "step_coverage",
            "ko_count",
            "kos_present",
            "genes_present",
        ],
    )
    return coverage_df


def make_module_coverage_frame(
        annotations, module_nets, groupby_column=DEFAULT_GROUPBY_COLUMN
):
    # go through each scaffold to check for modules
    module_coverage_dict = dict()
    for group, frame in annotations.groupby(groupby_column, sort=False):
        module_coverage_dict[group] = make_module_coverage_df(frame, module_nets)
    module_coverage = pd.concat(module_coverage_dict)
    module_coverage.index = module_coverage.index.set_names(["genome", "module"])
    return module_coverage.reset_index()


def pairwise(iterable):
    """s -> (s0, s1), (s1, s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def first_open_paren_is_all(str_):
    """Go through string and return true"""
    curr_level = 1
    for char in str_[1:-1]:
        if char == ")":
            curr_level -= 1
        elif char == "(":
            curr_level += 1
        if curr_level == 0:
            return False
    return True


def split_into_steps(definition, split_char=" "):
    """Very fancy split on string of chars"""
    curr_level = 0
    step_starts = [-1]
    for i, char in enumerate(definition):
        if char == "(":
            curr_level += 1
        if char == ")":
            curr_level -= 1
        if (curr_level == 0) and (char in split_char):
            step_starts.append(i)
    step_starts.append(len(definition))
    steps = list()
    for a, b in pairwise(step_starts):
        step = definition[a + 1: b]
        if step.startswith("(") and step.endswith(")"):
            if first_open_paren_is_all(step):
                step = step[1:-1]
        steps.append(step)
    return steps


def is_ko(ko):
    return re.match(KO_REGEX, ko) is not None


def make_module_network(
        definition, network: nx.DiGraph = None, parent_nodes=("start",)
):
    # TODO: Figure out how to add 'end' step to last step at end
    if network is None:
        network = nx.DiGraph()
    last_steps = []
    for step in split_into_steps(definition, ","):
        prev_steps = parent_nodes
        for substep in split_into_steps(step, "+"):
            if is_ko(substep):
                for prev_step in prev_steps:
                    network.add_edge(prev_step, substep)
                prev_steps = [substep]
            else:
                network, prev_steps = make_module_network(substep, network, prev_steps)
        last_steps += prev_steps
    return network, last_steps


def get_module_coverage(module_net: nx.DiGraph, genes_present: set):
    max_coverage = -1
    max_coverage_genes = list()
    max_coverage_missing_genes = list()
    max_path_len = 0
    for net_path in nx.all_simple_paths(module_net, source="start", target="end"):
        net_path = set(net_path[1:-1])
        overlap = net_path & genes_present
        coverage = len(overlap) / len(net_path)
        if coverage > max_coverage:
            max_coverage = coverage
            max_coverage_genes = overlap
            max_coverage_missing_genes = net_path - genes_present
            max_path_len = len(net_path)
    return (
        max_path_len,
        len(max_coverage_genes),
        max_coverage,
        max_coverage_genes,
        max_coverage_missing_genes,
    )


def make_etc_coverage_df(
        etc_module_df,
        annotation_ids_by_row: pd.DataFrame,
        groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    etc_coverage_df_rows = list()
    for _, module_row in etc_module_df.iterrows():
        definition = module_row["definition"]
        # remove optional subunits
        definition = re.sub(r"-K\d\d\d\d\d", "", definition)
        module_net, _ = make_module_network(definition)
        # add end node
        no_out = [
            node for node in module_net.nodes() if module_net.out_degree(node) == 0
        ]
        for node in no_out:
            module_net.add_edge(node, "end")
        # go through each genome and check pathway coverage
        for group, frame in annotation_ids_by_row.groupby(groupby_column):
            # get annotation genes
            grouped_ids = set(get_all_annotation_ids(frame).keys())
            (
                path_len,
                path_coverage_count,
                path_coverage_percent,
                genes,
                missing_genes,
            ) = get_module_coverage(module_net, grouped_ids)
            complex_module_name = "Complex %s: %s" % (
                module_row["complex"].replace("Complex ", ""),
                module_row["module_name"],
            )
            etc_coverage_df_rows.append(
                [
                    module_row["module_id"],
                    module_row["module_name"],
                    module_row["complex"].replace("Complex ", ""),
                    group,
                    path_len,
                    path_coverage_count,
                    path_coverage_percent,
                    ",".join(sorted(genes)),
                    ",".join(sorted(missing_genes)),
                    complex_module_name,
                ]
            )
    return pd.DataFrame(etc_coverage_df_rows, columns=ETC_COVERAGE_COLUMNS)


def make_functional_df(
        annotation_ids_by_row,
        function_heatmap_form,
        groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    # clean up function heatmap form
    function_heatmap_form = function_heatmap_form.apply(
        lambda x: x.str.strip() if x.dtype == "object" else x
    )
    function_heatmap_form = function_heatmap_form.fillna("")
    # build dict of ids per genome
    genome_to_id_dict = dict()
    for genome, frame in annotation_ids_by_row.groupby(groupby_column, sort=False):
        id_list = get_all_annotation_ids(frame).keys()
        genome_to_id_dict[genome] = set(id_list)
    # build long from data frame
    rows = list()
    for function, frame in function_heatmap_form.groupby("function_name", sort=False):
        for bin_name, id_set in genome_to_id_dict.items():
            presents_in_bin = list()
            functions_present = set()
            for _, row in frame.iterrows():
                function_id_set = set(
                    [i.strip() for i in row.function_ids.strip().split(",")]
                )
                present_in_bin = id_set & function_id_set
                functions_present = functions_present | present_in_bin
                presents_in_bin.append(len(present_in_bin) > 0)
            function_in_bin = np.all(presents_in_bin)
            row = frame.iloc[0]
            rows.append(
                [
                    row.category,
                    row.subcategory,
                    row.function_name,
                    ", ".join(functions_present),
                    "; ".join(get_ordered_uniques(frame.long_function_name)),
                    "; ".join(get_ordered_uniques(frame.gene_symbol)),
                    bin_name,
                    function_in_bin,
                    "%s: %s" % (row.category, row.function_name),
                ]
            )
    return pd.DataFrame(
        rows,
        columns=list(function_heatmap_form.columns)
                + ["genome", "present", "category_function_name"],
    )


# TODO: refactor this to handle splitting large numbers of genomes into multiple heatmaps here
def fill_product_dfs(
        annotations,
        module_nets,
        etc_module_df,
        function_heatmap_form,
        annotation_ids_by_row: pd.DataFrame,
        groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    module_coverage_frame = make_module_coverage_frame(
        annotations, module_nets, groupby_column
    )

    # make ETC frame
    etc_coverage_df = make_etc_coverage_df(
        etc_module_df, annotation_ids_by_row, groupby_column
    )

    # make functional frame
    function_df = make_functional_df(
        annotation_ids_by_row,
        function_heatmap_form,
        groupby_column,
    )

    return module_coverage_frame, etc_coverage_df, function_df


def make_product_df(module_coverage_frame, etc_coverage_df, function_df):
    liquor_df = pd.concat(
        [
            module_coverage_frame.pivot(
                index="genome", columns="module_name", values="step_coverage"
            ),
            etc_coverage_df.pivot(
                index="genome",
                columns="complex_module_name",
                values="percent_coverage",
            ),
            function_df.pivot(
                index="genome",
                columns="category_function_name",
                values="present",
            ),
        ],
        axis=1,
        sort=False,
    )
    return liquor_df


def get_phylum_and_most_specific(taxa_str):
    taxa_ranks = [i[3:] for i in taxa_str.split(";")]
    phylum = taxa_ranks[1]
    most_specific_rank = TAXONOMY_LEVELS[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    most_specific_taxa = taxa_ranks[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    if most_specific_rank == "d":
        return "d__%s;p__" % most_specific_taxa
    if most_specific_rank == "p":
        return "p__%s;c__" % most_specific_taxa
    else:
        return "p__%s;%s__%s" % (phylum, most_specific_rank, most_specific_taxa)


def make_strings_no_repeats(genome_taxa_dict: dict):
    labels = dict()
    seen = Counter()
    for genome, taxa_string in genome_taxa_dict.items():
        final_taxa_string = "%s_%s" % (taxa_string, str(seen[taxa_string]))
        seen[taxa_string] += 1
        labels[genome] = final_taxa_string
    return labels


def get_annotation_ids_by_row(data):
    functions = {i: j for i, j in ID_FUNCTION_DICT.items() if i in data.columns}
    missing = [i for i in ID_FUNCTION_DICT if i not in data.columns]
    logger.info(
        "Note: the fallowing id fields "
        f"were not in the annotations file and are not being used: {missing},"
        f" but these are {list(functions.keys())}"
    )
    out = data.apply(
        lambda x: {
            i
            for k, v in functions.items()
            if not pd.isna(x[k])
            for i in v(str(x[k]))
            if not pd.isna(i)
        },
        axis=1,
    )
    return out


def get_all_annotation_ids(data):
    annotation_series = data[DBSETS_COL]
    annotation_series.apply(list)
    out = Counter(chain(*annotation_series.values))
    return out


def rename_genomes_to_taxa(function_df, labels):
    function_df = function_df.copy()
    new_genome_column = [labels[i] for i in function_df["genome"]]
    function_df["genome"] = pd.Series(new_genome_column, index=function_df.index)
    return function_df


class DramUsageError(Exception):  # TODO maybe remove or make more specific
    "Raised when dram is not used corectly, usally it means you are missing a step"
    pass


def get_ordered_uniques(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x) or pd.isna(x))]


def main(annotations_tsv_path, groupby_column=DEFAULT_GROUPBY_COLUMN, genomes_per_product=DEFAULT_GENOMES_PER_PRODUCT,
         make_big_html=DEFAULT_MAKE_BIG_HTML):
    logger = logging.getLogger("dram2_log.viz")
    dram_config = {}
    output_dir = Path.cwd().resolve()
    annotations = pd.read_csv(annotations_tsv_path, sep="\t", index_col=0)
    print(annotations.columns)
    print(annotations.head())

    # db_kits_with_ids = [i for i in DB_KITS if i.selectable and i.can_get_ids]
    # print(f"{DB_KITS=}")
    # print(f"{db_kits_with_ids=}")
    # db_kits = [i(dram_config, logger) for i in db_kits_with_ids]
    # print(f"{db_kits=}")

    db_id_sets: pd.Series = get_annotation_ids_by_row(
        annotations
    )
    annotation_ids_by_row = annotations.copy()
    annotation_ids_by_row[DBSETS_COL] = db_id_sets

    module_steps_form = get_distillate_sheet(
        MODULE_STEPS_FORM_TAG, dram_config, logger
    )
    print(module_steps_form.columns)
    print(module_steps_form.head())
    etc_module_df = get_distillate_sheet(ETC_MODULE_DF_TAG, dram_config, logger)
    print(etc_module_df.columns)
    print(etc_module_df.head())
    function_heatmap_form = get_distillate_sheet(
        FUNCTION_HEATMAP_FORM_TAG, dram_config, logger
    )
    # make product
    if "bin_taxonomy" in annotations:
        genome_order = get_ordered_uniques(
            annotations.sort_values("bin_taxonomy")[groupby_column]
        )
        # if gtdb format then get phylum and most specific
        if all(
                [
                    i[:3] == "d__" and len(i.split(";")) == 7
                    for i in annotations["bin_taxonomy"].fillna("")
                ]
        ):
            taxa_str_parser = get_phylum_and_most_specific
        # else just throw in what is there
        else:

            def taxa_str_parser(x):
                return x

        labels = make_strings_no_repeats(
            {
                row[groupby_column]: taxa_str_parser(row["bin_taxonomy"])
                for _, row in annotations.iterrows()
            }
        )
    else:
        genome_order = get_ordered_uniques(
            annotations.sort_values(groupby_column)[groupby_column]
        )
        labels = None

    # make module coverage frame
    module_nets = {
        module: build_module_net(module_df)
        for module, module_df in module_steps_form.groupby("module")
        if module in HEATMAP_MODULES
    }

    module_coverage_df, etc_coverage_df, function_df = fill_product_dfs(
        annotations=annotations,
        module_nets=module_nets,
        etc_module_df=etc_module_df,
        function_heatmap_form=function_heatmap_form,
        annotation_ids_by_row=annotation_ids_by_row,
        groupby_column=groupby_column
    )
    product_df = make_product_df(module_coverage_df, etc_coverage_df, function_df)

    # module_coverage_df.to_csv("module_coverage_df.csv", index=False)
    # etc_coverage_df.to_csv("etc_coverage_df.csv", index=False)
    # function_df.to_csv("function_df.csv", index=False)

    #

    if labels is not None:
        function_df = rename_genomes_to_taxa(function_df, labels)

    if make_big_html or len(genome_order) < GENOMES_PRODUCT_LIMIT:
        # product = make_product_heatmap(
        plot = make_product_heatmap(
            module_coverage_df,
            etc_coverage_df,
            function_df,
            genome_order,
            labels,
        )
        plot.save(output_dir / "product.html", resources=INLINE_RESOURCES)
    product_df = make_product_df(module_coverage_df, etc_coverage_df, function_df)
    product_df.to_csv(output_dir / "product.tsv", sep="\t", index=False)
    logger.info("Completed visualization")


def make_product_heatmap(
        module_coverage_df,
        etc_coverage_df,
        function_df,
        genome_order,
        labels):
    charts = []
    titles = []
    module_charts = []
    module_charts.append(
        (heatmap(module_coverage_df, x_col="module_name", y_col="genome", c_col="step_coverage",
                 tooltip_cols=["genome", "module_name", "steps", "steps_present"],
                 title="Module",
                 ),
         "Module")
    )
    titles.append("Module")

    etc_tooltip_cols = ["genome", "module_name", "path_length", "path_length_coverage", "genes", "missing_genes"]
    etc_charts = []
    for i, (etc_complex, frame) in enumerate(etc_coverage_df.groupby("complex")):
        etc_charts.append(
            (heatmap(frame, x_col="module_name", y_col="genome", c_col="percent_coverage",
                     tooltip_cols=etc_tooltip_cols, y_axis_location=None,
                     title=etc_complex,
                     # margin=(0, 20, 0, 20)
                     ),
             etc_complex)
        )
        titles.append(etc_complex)
    add_colorbar(etc_charts[-1][0])

    func_tooltip_cols = ["genome", "category", "subcategory", ("Function IDs", "@function_ids"), "function_name",
                         "long_function_name", "gene_symbol"]
    function_df["present"] = function_df["present"].astype(str)
    function_charts = []
    for i, (group, frame) in enumerate(function_df.groupby("category", sort=False)):
        kw = dict()
        if i == len(function_df["category"].unique()) - 1:
            kw["rect_kw"] = dict(legend_field="present")
        function_charts.append(
            (heatmap(frame, x_col="function_name", y_col="genome", c_col="present",
                     tooltip_cols=func_tooltip_cols,
                     y_axis_location=None,
                     # rect_kw=dict(legend_field="present"),
                     title=group,
                     **kw
                     ),
             group)
        )
    l = function_charts[-1][0].legend[0]
    function_charts[-1][0].legend[:] = []
    function_charts[-1][0].add_layout(l, 'right')

    charts.append(format_chart_group(chart_group=function_charts))

    charts = [
        format_chart_group([bokeh_pane(p, title) for p, title in module_charts]),
        format_chart_group([bokeh_pane(p, title) for p, title in etc_charts], title="ETC Complexes"),
        format_chart_group([bokeh_pane(p, title) for p, title in function_charts]),
    ]

    plot = pn.Row(*charts)
    return plot


def bokeh_pane(p, title="", p_kw=None, title_kw=None):
    p_kw = p_kw or dict()
    title_kw = title_kw or dict()

    # return pn.pane.Bokeh(p, **p_kw)
    return pn.Column(
        # pn.pane.Markdown(f"## {title}", align="center"),
        pn.Row(pn.pane.Bokeh(p, **p_kw), align="center"),
    )


def format_chart_group(chart_group, title=""):
    if not isinstance(chart_group, list) or isinstance(chart_group, tuple):
        chart_group = [chart_group]
    return pn.Column(
        pn.pane.Markdown(f"## {title}", align="center"),
        pn.Row(*chart_group)
    )


def add_colorbar(p):
    color_bar = ColorBar(color_mapper=LinearColorMapper(palette=tuple(reversed(PALETTE_CONTINUOUS)), low=0, high=1),
                         height=HEATMAP_CELL_HEIGHT * 30)
    p.add_layout(color_bar, 'right')
    return p


def heatmap(df, x_col, y_col, c_col, tooltip_cols, title="", rect_kw=None, **fig_kwargs, ):
    rect_kw = rect_kw or {}
    df = df.sort_values(by=[y_col], ascending=False)

    tooltips = []
    for col in tooltip_cols:
        if isinstance(col, tuple):
            tooltips.append(col)
        else:
            tooltips.append((col.replace("_", " ").title(), f"@{col}"))

    half_ttl_ln = len(title) / 2
    if half_ttl_ln > len(
            df[x_col].unique()) and "min_border_left" not in fig_kwargs and "min_border_right" not in fig_kwargs:
        left_over = (half_ttl_ln - len(df[x_col].unique()))
        fig_kwargs["min_border_left"] = int(left_over / 1.4) * HEATMAP_CELL_WIDTH
        fig_kwargs["min_border_right"] = int(left_over / 1.4) * HEATMAP_CELL_WIDTH

    p = figure(
        frame_width=HEATMAP_CELL_WIDTH * len(df[x_col].unique()),
        frame_height=HEATMAP_CELL_WIDTH * len(df[y_col].unique()),

        x_range=sorted(list(df[x_col].unique())), y_range=list(df[y_col].unique()),
        tools="hover",
        toolbar_location=None,
        tooltips=tooltips,
        title=title,
        **fig_kwargs
    )

    if df[c_col].dtype == float:
        palette = tuple(reversed(PALETTE_CONTINUOUS))
        fill_color = linear_cmap(c_col, palette=palette, low=df[c_col].min(), high=df[c_col].max())
    else:
        factors = sorted(df[c_col].unique())
        max_factors = max(PALETTE_CATEGORICAL.keys())
        palette = PALETTE_CATEGORICAL[max(len(factors), 3)] if len(factors) <= max_factors else PALETTE_CONTINUOUS
        fill_color = factor_cmap(c_col, palette=tuple(reversed(palette)), factors=factors)

    p.rect(x=x_col, y=y_col,
           width=0.9, height=0.9,
           source=df,
           fill_alpha=0.9,
           color=fill_color,
           **rect_kw
           )

    p.title.align = "center"

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "14px"
    p.xaxis.major_label_orientation = pi / 2

    return p


if __name__ == "__main__":
    main(Path(__file__).resolve().parent / "data" / "annotations.tsv")
