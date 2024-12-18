CONDA_ENV_DEPS = file(params.main_environment).text.split('dependencies:')[1]

DNA_HEADER = """
===========================================
                   $workflow.manifest.name
(==(     )==)                 (==(     )==)
 `-.`. ,',-'                   `-.`. ,',-'
    _,-'"                         _,-'"
 ,-',' `.`-.                    ,-',' `.`-.
(==(     )==)                  (==(     )==)
 `-.`. ,',-'                    `-.`. ,',-'
    _,-'"                         _,-'"
 ,-',' `.`-.                   ,-',' `.`-.
(==(     )==)                 (==(     )==)
                   $workflow.manifest.version
===========================================
"""

DEPENDENCIES = """
    Software versions used:

    $CONDA_ENV_DEPS
"""

LICENSE = "MIT License. Micobial Ecosystems Lab, Colorado State University Fort Collins. 2024 (last updated 2024)"

REQUIRED_OPTIONS = """
REQUIRED $workflow.manifest.name profile options:
    -profile                STRING  <conda, conda_slurm, singularity, singularity_conda>
                                        Runs $workflow.manifest.name either using Conda (must be installed) or Singularity (must be installed).
                                        Runs $workflow.manifest.name with no scheduling or scheduling via SLURM.
                                        See SLURM options in full help menu.
"""
MAIN_HELP = """
Description:
    The purpose of $workflow.manifest.name is to provide FASTA annotation, across a vast array of databases, with expertly-currated distillation.
    $workflow.manifest.name can be used to call, annotate and distill annotations from input FASTA files.
    Call, annotate and distill can be run together or, each can be run idependently.

Bring up help menu:
    nextflow run $workflow.manifest.name --help (--h)

Bring up versions menu:
    nextflow run $workflow.manifest.name --version (--v)

Usage:
    nextflow run $workflow.manifest.name --rename --call --annotate --use_<database(s) --distill_topic <distillate(s)>

    Call genes using input fastas (use --rename to rename FASTA headers):
        nextflow run $workflow.manifest.name --call --rename --input_fasta <path/to/fasta/directory/>

    Annotate called genes using input fastas:
        nextflow run $workflow.manifest.name --annotate --input_genes <path/to/called/genes/directory>

    Distill using input annotations:
        nextflow run $workflow.manifest.name --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv>

    (Combined): Call, annotate and distill input fasta files:
        nextflow run $workflow.manifest.name --rename --call --annotate --use_<database(s) --distill_topic <distillate(s)>

    (Real) example: (on multiple lines for clarity)
    nextflow run $workflow.manifest.name --input_fasta ../test_data/
        --outdir $workflow.manifest.name-test-data-Feb012024/
        --call --rename
        --annotate --use_uniref --use_kegg --use_merops --use_viral --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur
        --add_annotations ../test-data/old-DRAM-annotations.tsv
        --distill_topic 'carbon transport energy' --distill_ecosystem 'eng_sys ag'
        --distill_custom bin/forms/distill_sheets/test.tsv -resume --slurm_node zenith
        --trnas ../test-data/trnas.tsv
        --rrnas ../test-data/rrnas.tsv
        --bin_quality ../test-data/checkM1-test-data.tsv
        --taxa ../test-data/gtdbtk.bac120.summary.tsv
        --generate_gff
        --generate_gbk
        --threads 5
        -with-report -with-trace -with-timeline

Main $workflow.manifest.name Operations:
    --call        : Call genes using prodigal
    --annotate    : Annotate called genes using downloaded databases
    --distill     : Distill the annotations into a multi-sheet distillate.xlsx
    --product     : Generate a product heatmap of the annotations
    --format_kegg : Format KEGG database for use in $workflow.manifest.name. Standalone operation, will exit after completion.
"""

CALL_OPTIONS = """
Call options:
    --call                          Call genes on the input FASTA files using Prodigal.

    --input_fasta           PATH    <path/to/fasta/directory/>
                                        Directory containing input fasta files.
                                        Default: <./input_fasta/*.fa*>

    --rename                OPTION  Rename FASTA headers based on file name.
                                        Example: sample1.fa --> (fasta header renamed to) > sample1......
                                        Why? $workflow.manifest.name output is focused on scaffolds/contigs with respect to each provided input sample.
                                            Thus, without renaming FASTA headers, the individual scaffolds/contigs will not be distinguashable.
                                            *If you have already renamed your FASTA headers, do not include '--call'.

    --prodigal_mode         STRING  <single|meta>
                                        Default: 'single'

    --prodigal_tras_table   NUMBER  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
                                        Specify a translation table to use (default: '1').

    --min_contig_len        NUMBER  <number in base pairs>
                                        Default: '2500'
"""

ANNOTATE_OPTIONS = """
Annotate options:
    --annotate                         Annotate called genes using downloaded databases

    --use_<db-name>         STRING   <camper|cant_hyd|dbcan|fegenie|kegg|kofam|merops|methyl|heme|pfam|sulfur|uniref]
                                        Specify databases to use. Can use more than one. Can be used in combination with --use_dbset.

    --use_dbset             STRING  <metabolism_kegg_set|metabolism_set|adjectives_kegg_set|adjectivs set>
                                        metabolism_kegg_set = kegg, dbcan, merops, pfam, heme
                                        metabolism_set      = kofam, dbcan, merops, pfam, heme
                                        adjectives_kegg_set = kegg, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                        adjectives_set      = kofam, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                        *Only one set can be used. Can be used in combination with --use_[db-name]

    --input_fasta           PATH    <path/to/fasta/directory/>
                                        Directory containing input fasta files.
                                        Default: './input_fasta/'
                                        Either '--input_fasta' or '--input_genes' may be used - not both.

    --input_genes           PATH    <path/to/called/genes/directory/>
                                        Directory containing called genes (.faa)

    --add_annotations       PATH    <path/to/old-annoations.tsv>
                                        Used to add in old annotations to the current run. (See example for format.)

    --generate_gff          OPTION  Will generate an output GFF for each sample based on the raw-annotations.tsv.

    --generate_gbk          OPTION  Will generate an output GBK for each sample based on the raw-annotations.tsv.
"""

DISTILL_OPTIONS = """
Distill options:
    --distill                          Distill the annotations into a multi-sheet distillate.xlsx

    --annotations           PATH     <path/to/annotations.tsv>
                                        Required if you are running distill without --call and --annotate.

    --rrnas                 PATH    <path/to/rRNA.tsv> (See example for format.)
                                        rRNA information will be included in distill output.

    --trnas                 PATH    <path/to/tRNA.tsv> (See example for format.)
                                        tRNA information will be included in distill output.

    --bin_quality           PATH    <path/to/bin-quality.tsv> (See example for format.)
                                        CheckM and CheckM2 compatible.

    --taxa                  PATH    <path/to/bin-taxonomy.tsv>
                                        Compatible with GTDB. (See example for format.)

    --distill_topic         STRING  <carbon|energy|misc|nitrogen|transport> OR <default = carbon, energy, misc, nitrogen, transport>
                                        If more than one topic included, they must be enclosed in single quotes

    --distill_ecosystem     STRING  <eng_sys|ag>
                                        If more than one ecosystem included, they must be enclosed in single quotes

    --distill_custom        STRING  <path/to/custom_distillate.tsv> (See example for format and options.)
                                        As of now, only one custom distillate may be included.
"""

PRODUCT_OPTIONS = """
Product options:
    --product                          Generate a product visualization of the annotations and save the output to the output directory.

    --annotations           PATH     <path/to/annotations.tsv>
                                        Required if you are running product without --call and --annotate.

    --groupby_column        STRING   Column to to group by in the annotations file for etc and function groupings
                                        Default: 'fasta'

    --module_steps_form     PATH     <path/to/module_steps_form.tsv>
                                        override default module steps form database TSV

    --etc_module_form       PATH     <path/to/etc_module_form.tsv>
                                        override default etc module form database TSV

    --function_heatmap_form PATH     <path/to/function_heatmap_form.tsv>
                                        override default function heatmap form database TSV
"""

TREE_OPTIONS = """
Tree options:
    --trees              OPTION  Will run Trees. (This option is not advised as $workflow.manifest.name Trees is in development.)
"""

FORMAT_KEGG_DB_OPTIONS = """
Format KEGG Database options:
    --format_kegg                      Format KEGG database for use in $workflow.manifest.name. Standalone operation, will exit after completion.

    --kegg_pep_loc        PATH    <path/to/kegg.pep>
                                        Path to and of the gene fasta files that are provided by the
                                        KEGG FTP server or a concatenated version of them
    --gene_ko_link_loc    PATH    <path/to/genes_ko_link>
                                        Path to the genes_ko_link file downloaded from KEGG

    --kegg_db            PATH    <path/to/kegg_db>
                                        Output path to the KEGG database directory. This is where the
                                        formatted KEGG database will be stored. Default: `${params.kegg_db}`
                                        relative to $workflow.manifest.name directory.
    --kegg_download_date STRING  <YYYY-MM-DD>
                                        The date the KEGG database was downloaded. If not provided, the
                                        current date will be used.
    --skip_gene_ko_link  OPTION         Skip the gene_ko_link file. If you are using an older version of KEGG that does 
                                        supply the gene_ko_link file you can use this option to skip the gene_ko_link 
                                        file. Otherwise, the gene_ko_link file is required.
"""

ADJECTIVES_OPTIONS = """
Adjectives options:
    --annotations_tsv_path PATH     Location of an annotations.tsv. You don't
                                    need to use this option if you are using the
                                    output_dir for dram with a project_config.
                                    If you use this option, you must also use
                                    the force flag to bypass the safeguards that
                                    prevent you from running distill with
                                    insufficient data.
    --adjectives_tsv_path PATH      Location of the output adjectives.tsv. if
                                    you leave this blank the adjectives.tsv file
                                    will be put in the output directory.
    -a, --adjectives TEXT           A list of adjectives, by name, to evaluate.
                                    This limits the number of adjectives that
                                    are evaluated, and is faster.
    -p, --plot_adjectives TEXT      A list of adjectives, by name, to plot. This
                                    limits the number of adjectives that are
                                    plotted and is probably needed for speed.
    -g, --plot_genomes TEXT
    --plot_path PATH                will become a folder of output plots, no
                                    path no plots.
    --strainer_tsv PATH             The path for a tsv that will pass to
                                    strainer to filter genes. The only option at
                                    this time is ‘pgtb’ for positive genes that
                                    are on true bugs.
    --strainer_type PATH            The type of process that should make the
                                    strainer file.
    --debug_ids_by_fasta_to_tsv PATH
                                    This is a tool to debug the list of IDs
                                    found by DRAM it is mostly for experts.
    --user_rules_tsv PATH           This is an optional path to a rules file
                                    with strict formatting. It will overwrite
                                    the original rules file that is stored with
                                    the script.
    --show_rules_path               Show the path to the default rules path.
    --list_name                     List the names for all adjectives_tsv that
                                    are available, you can pass these names to
                                    limit the adjectives that are evaluated
    --list_id                       List the names for all adjectives_tsv that
                                    are available, you can pass these names to
                                    limit the adjectives that are evaluated,
    --input_fasta           PATH    <path/to/fasta/directory/>
                                    Directory containing input fasta files.
                                    Default: <./input_fasta/*.fa*>
    """

GENERAL_OPTIONS = """
General options:
    --outdir                PATH    <path/to/output/directory>
                                        Default: './DRAM_output/'

    --threads               NUMBER  Number of threads to use for processing.
                                        Default: `${params.threads}`

    --slurm_node            string  <node_name>
                                        Example --slurm_queue c001

    --slurm_queue           string  <slurm partition name>
                                        Example:  --slurn_queue 'smith-hi,smith-low'
"""

CALL_DESCRIPTION = """
Call description: The purpose of $workflow.manifest.name --call is to call genes on input FASTA files.

Usage:

    Call genes using input fastas:
        nextflow run $workflow.manifest.name --call --input_fasta <path/to/fasta/directory/> --outdir <path/to/output/directory/> --threads <threads>
"""

ANNOTATE_DESCRIPTION = """
Annotate description: The purpose of $workflow.manifest.name '--annotate' is to annotate called genes on input (nucleotide) FASTA (fa*) files.

Usage:

    Annotate called genes using input called genes and the KOFAM database:
        nextflow run $workflow.manifest.name --annotate --input_genes <path/to/called/genes/directory> --use_kofam

    Annotate called genes using input fasta files and the KOFAM database:
        nextflow run $workflow.manifest.name --annotate --input_fasta <path/to/called/genes/directory> --use_kofam
"""

DISTILL_DESCRIPTION = """
Distill description:    The purpose of $workflow.manifest.name --distill is to distill down annotations based on curated distillation summary form(s).
                        User's may also provide a custom distillate via --distill_custom <path/to/file> (TSV forms).
                        Distill can be ran independent of --call and --annotate however, annotations must be provided (--annotations <path/to/annotations.tsv>).
                        Optional tRNA, rRNA and bin quality may also be provided.

Usage:
    nextflow run $workflow.manifest.name --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/> --threads <threads>
    *Important: if more than one topic or ecosystem is included, they must be enclosed in single quotes. Example: --distill_topic 'carbon transport'

Example:
    Call and Annotate genes using input fastas and KOFAM database. Distill using carbon topic and AG ecosystem:
        nextflow run $workflow.manifest.name --input_fasta <path/to/fasta/directory/> --outdir <path/to/output/directory/> --call --annotate --distill_topic carbon --distill_ecosystem ag --threads <threads> --use_kofam
"""

PRODUCT_DESCRIPTION = """
Product description: The purpose of $workflow.manifest.name --product is to generate a product visualization of the annotations
		     and save the output to the output directory.

Usage:
    nextflow run $workflow.manifest.name --product --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/>

Example:
    Create heatmap product visualization from annotations file and save to output directory:
        nextflow run $workflow.manifest.name --product --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/>
"""

ADJECTIVES_DESCRIPTION = """
Annotate description: The purpose of $workflow.manifest.name '--adjectives' is to evaluate genes and describe their features.

Usage:
    TO BE ADDED LATER

THIS IS A WORK IN PROGRESS AND IS NOT YET COMPLETE
"""

FORMAT_KEGG_DB_DESCRIPTION = """
Format KEGG DB description: The purpose of $workflow.manifest.name '--format_kegg' is to format the raw KEGG database (pep files)
			    for use in $workflow.manifest.name. Because use of the KEGG database requires a license, the user must download the raw KEGG database
			    themselves. Currently this only supports formatting a concatenated version of the KEGG pep files see
			    (https://github.com/WrightonLabCSU/DRAM/issues/305) for guidance on how to prepare your pep files for formatting.

Usage:
    nextflow run $workflow.manifest.name --format_kegg --kegg_pep_loc <path/to/kegg.pep> --kegg_db <path/to/save/output/kegg_db/to>
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Version menu
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def version() {
log.info"""
${DNA_HEADER}
${DEPENDENCIES}
"""
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Help menus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
log.info """
${DNA_HEADER}
${LICENSE}
${MAIN_HELP}
${REQUIRED_OPTIONS}
${CALL_OPTIONS}
${ANNOTATE_OPTIONS}
${DISTILL_OPTIONS}
${PRODUCT_OPTIONS}
${TREE_OPTIONS}
${FORMAT_KEGG_DB_OPTIONS}
${GENERAL_OPTIONS}
"""
}

def callHelpMessage() {
log.info """
${DNA_HEADER}
${LICENSE}
${CALL_DESCRIPTION}
${REQUIRED_OPTIONS}
${CALL_OPTIONS}
${GENERAL_OPTIONS}
"""
}

def annotateHelpMessage() {
log.info """
${DNA_HEADER}
${LICENSE}
${ANNOTATE_DESCRIPTION}
${REQUIRED_OPTIONS}
${ANNOTATE_OPTIONS}
${GENERAL_OPTIONS}
"""
}

def distillHelpMessage() {
log.info """
${DNA_HEADER}
${LICENSE}
${DISTILL_DESCRIPTION}
${REQUIRED_OPTIONS}
${DISTILL_OPTIONS}
${GENERAL_OPTIONS}
"""
}

def productHelpMessage() {
log.info """
${DNA_HEADER}
${LICENSE}
${PRODUCT_DESCRIPTION}
${REQUIRED_OPTIONS}
${PRODUCT_OPTIONS}
${GENERAL_OPTIONS}
"""
}

def adjectivesHelpMessage() {
log.info """
${DNA_HEADER}
${LICENSE}
${ADJECTIVES_DESCRIPTION}
${REQUIRED_OPTIONS}
${ADJECTIVES_OPTIONS}
${GENERAL_OPTIONS}
"""
}

def formatKeggHelpMessage() {
log.info """
${DNA_HEADER}
${LICENSE}
${FORMAT_KEGG_DB_DESCRIPTION}
${REQUIRED_OPTIONS}
${FORMAT_KEGG_DB_OPTIONS}
${GENERAL_OPTIONS}
"""
}