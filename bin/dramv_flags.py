#!/usr/bin/env python
"""Compute DRAM-v AMG flags, auxiliary_score, and is_transposon columns.

Phase 1+2b+2c (FASTA-default): no VirSorter input required. Adds:

  - is_transposon (bool): the gene's Pfam hits intersect TRANSPOSON_PFAMS.
  - amg_flags (str): concatenated single-letter flags in v1 order
    M / K / E / V / A / P / T / F / B, then a trailing N (not in v1).
    V fires when the gene has a VOG hit whose VOGdb FunctionalCategory
    is Xr (viral replication) or Xs (virion structure). The N flag is
    informational — it marks genes whose ids appear in amg_database.tsv
    rows where essential_viral_function=TRUE, per Martin et al. 2025
    (doi:10.1038/s41564-025-02095-4). It does NOT propagate to M and is
    NOT used by any downstream filter yet.
  - auxiliary_score (int 1-5): VirSorter-flank-based confidence per v1's
    calculate_auxiliary_scores. Lower = more confident viral context. In
    pure-FASTA mode (no virsorter_categories supplied) every gene scores
    5, matching v1's fallback. The B flag downgrades 1/2/3 to 4.

All semantics ported verbatim from DRAM v1
mag_annotator/annotate_vgfs.py (commit 6cd68f9).
"""

import re
from pathlib import Path

import click
import polars as pl

from utils.dramv_constants import (
    CELL_ENTRY_CAZYS,
    TRANSPOSON_PFAMS,
    VIRAL_PEPTIDASES_MEROPS,
)
from utils.logger import get_logger
from rule_parser.src.rules import ID_EXPR_DICT

logger = get_logger(filename=Path(__file__).stem)

DEFAULT_LENGTH_FROM_END = 5000

# polars str.extract_all returns *full* matches, not capture groups. The
# rule_parser's ID_EXPR_DICT entries for Pfam wrap the id in brackets
# (e.g. "[PF01609.1]") because they were written assuming capture-group
# semantics, which means TRANSPOSON_PFAMS ∩ <those ids> never fires.
# Override locally to extract the bare id. Fix upstream is left for a
# follow-up rule_parser change.
_PFAM_ID_RE = r"PF\d{5}"
_VOG_ID_RE = r"VOG\d+"
_ID_EXPR_DICT = dict(ID_EXPR_DICT)
_ID_EXPR_DICT["pfam_hits"] = pl.col("pfam_hits").str.extract_all(_PFAM_ID_RE)
_ID_EXPR_DICT["pfam_id"] = pl.col("pfam_id").str.extract_all(_PFAM_ID_RE)
# vogdb_id / vogdb_ids aren't in the upstream rule_parser ID_EXPR_DICT, so
# register them here. vogdb_id is the bare best-hit name (e.g. "VOG00177")
# and vogdb_ids is "; "-joined all hits per hmm_parser. Column names follow
# hmm_parser's `{db_name}_id` / `{db_name}_ids` convention with db_name=vogdb.
_ID_EXPR_DICT["vogdb_id"] = pl.col("vogdb_id").cast(pl.Utf8).str.extract_all(_VOG_ID_RE)
_ID_EXPR_DICT["vogdb_ids"] = pl.col("vogdb_ids").cast(pl.Utf8).str.extract_all(_VOG_ID_RE)

# VOGdb functional categories that count as "viral" for the V flag, per DRAM v1
# (mag_annotator/annotate_vgfs.py:get_metabolic_flags, commit 6cd68f9):
#   Xr — viral replication
#   Xs — virion structure
# Xh (host benefit), Xp (host integration), and Xu (unknown) are excluded.
_V_FLAG_CATEGORIES = frozenset({"Xr", "Xs"})

# VirSorter category strings recognised by the auxiliary_score algorithm
# (mag_annotator/annotate_vgfs.py, commit 6cd68f9). Categories 0/3 are phage
# and prophage hallmark genes; 1/4 are phage and prophage viral-like genes.
# Category 2 (uncategorized) and any None do not contribute.
_VIRSORTER_HALLMARK_CATEGORIES = frozenset({"0", "3"})
_VIRSORTER_VIRAL_LIKE_CATEGORIES = frozenset({"1", "4"})


def read_scaffold_lengths(fasta_path: str) -> dict[str, int]:
    """Return {scaffold_id: length_bp} from a (multi-)FASTA."""
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_len = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.rstrip())
    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


def _split_semicolon_ids(s: pl.Series) -> set[str]:
    out: set[str] = set()
    for cell in s.drop_nulls():
        for tok in str(cell).split(";"):
            tok = tok.strip()
            if tok:
                out.add(tok)
    return out


def build_amg_id_sets(
    amg_db_path: Path,
) -> tuple[set[str], set[str], set[str]]:
    """Return (all_amgs, verified_amgs, essential_amgs).

    all_amgs       — KO ∪ EC ∪ PFAM ids from every row of the AMG db.
    verified_amgs  — same union restricted to rows with verified=TRUE.
    essential_amgs — same union restricted to rows with
                     essential_viral_function=TRUE (Martin et al. 2025
                     paper-cautioned). The column is optional; if absent,
                     this set is empty.
    """
    db = pl.read_csv(amg_db_path, separator="\t", infer_schema_length=10_000)
    amgs = (
        _split_semicolon_ids(db["KO"])
        | _split_semicolon_ids(db["EC"])
        | _split_semicolon_ids(db["PFAM"])
    )

    def _ids_where(col: str) -> set[str]:
        if col not in db.columns:
            return set()
        sub = db.filter(
            pl.col(col).cast(pl.Utf8).str.strip_chars().str.to_uppercase() == "TRUE"
        )
        return (
            _split_semicolon_ids(sub["KO"])
            | _split_semicolon_ids(sub["EC"])
            | _split_semicolon_ids(sub["PFAM"])
        )

    return amgs, _ids_where("verified"), _ids_where("essential_viral_function")


def build_viral_vog_ids(vog_list_path: Path) -> set[str]:
    """Return the set of VOG ids whose FunctionalCategory is Xr or Xs.

    `vog_list_path` is VOGdb's vog_annotations_latest.tsv(.gz) — a 5-col TSV
    (GroupName, ProteinCount, SpeciesCount, FunctionalCategory,
    ConsensusFunctionalDescription). polars handles the .gz transparently.
    Categories are the v1 VOGdb codes (Xh/Xp/Xr/Xs/Xu); only Xr and Xs
    count as viral for the V flag. A category cell can hold multiple
    one-character codes concatenated (e.g. "XrXs"); we match on substring.
    """
    df = pl.read_csv(vog_list_path, separator="\t", infer_schema_length=10_000)
    name_col = "GroupName" if "GroupName" in df.columns else df.columns[0]
    cat_col = (
        "FunctionalCategory" if "FunctionalCategory" in df.columns else df.columns[3]
    )
    sub = df.filter(
        pl.any_horizontal(
            [pl.col(cat_col).cast(pl.Utf8).str.contains(c) for c in _V_FLAG_CATEGORIES]
        )
    )
    return {str(v) for v in sub[name_col].drop_nulls().to_list()}


def build_metabolic_genes(distill_sheets_dir: Path) -> set[str]:
    """Union of gene_id values across distill sheets where potential_amg=TRUE."""
    metabolic: set[str] = set()
    for f in sorted(distill_sheets_dir.glob("distill_*.tsv")):
        df = pl.read_csv(f, separator="\t", infer_schema_length=10_000)
        if "potential_amg" not in df.columns or "gene_id" not in df.columns:
            continue
        sub = df.filter(
            pl.col("potential_amg").cast(pl.Utf8).str.strip_chars().str.to_uppercase() == "TRUE"
        )
        for gid in sub["gene_id"].drop_nulls().to_list():
            gid = str(gid).strip()
            if gid:
                metabolic.add(gid)
    return metabolic


_GENOMAD_HALLMARK_SUFFIXES = frozenset({"VV", "Vv"})
_GENOMAD_VIRAL_LIKE_SUFFIXES = frozenset({"vV", "vv"})

# DRAM scaffold name on a CheckV-trimmed provirus carries an offset suffix:
#   <parent_assembly_contig>|provirus_<start>_<end>
# DRAM calls genes on the trimmed sub-sequence so gene start_position is
# 1-based within the provirus. To match a geNomad gene called on the parent
# assembly contig, translate via assembly_pos = start_offset + dram_pos - 1.
_PROVIRUS_RE = re.compile(r"^(.+)\|provirus_(\d+)_(\d+)$")
_GENOMAD_GENE_NUM_RE = re.compile(r"_\d+$")


def _is_truthy(v) -> bool:
    """geNomad's boolean columns serialise as 0/1 ints in modern releases and
    TRUE/FALSE in older ones. Handle both, plus polars-decoded Python bool."""
    if v is True or v == 1:
        return True
    if isinstance(v, str):
        s = v.strip().upper()
        return s == "TRUE" or s == "1"
    return False


def _parse_dram_scaffold(scaffold: str) -> tuple[str, int]:
    """Return (parent_assembly_contig, start_offset_1based).

    For provirus contigs (`<contig>|provirus_<start>_<end>`) the offset is
    parsed from the suffix. For whole-virus contigs the parent is the
    scaffold itself and the offset is 1, so translation is a no-op.
    """
    m = _PROVIRUS_RE.match(scaffold)
    if m:
        return m.group(1), int(m.group(2))
    return scaffold, 1


def _parse_genomad_parent_contig(gene_id: str) -> str:
    """Strip the trailing _<gene_number> from a geNomad gene id."""
    return _GENOMAD_GENE_NUM_RE.sub("", gene_id)


def parse_genomad_genes_tsv(path: Path) -> pl.DataFrame:
    """Adapter from geNomad's *_genes.tsv to VirSorter-style category codes.

    Two complementary signals (modern geNomad ≥1.5):

      - virus_hallmark (0/1 int)             → "0"  hallmark
      - marker name suffix .VV / .Vv         → "0"  hallmark
      - marker name suffix .vV / .vv         → "1"  viral-like

    Anything else (host marker, plasmid_hallmark, NA marker) is dropped.
    The phage / prophage distinction (v1's 0 vs 3, 1 vs 4) is collapsed —
    auxiliary_score treats {0,3} and {1,4} as equivalent, so the adapter
    emits only "0" or "1".

    Returns a DataFrame with columns (gene, contig, start, end, category) —
    one row per *kept* gene. `gene` is the full geNomad gene id; `contig`
    is the parent assembly contig (gene id with the trailing _<num>
    stripped). Empty rows are dropped, so no-signal samples produce an
    empty frame, not None.

    Joining with DRAM annotations is handled by
    `build_virsorter_categories_from_genomad`. It tries direct gene-id
    equality first — the common case when both tools name provirus
    contigs the same way and use provirus-local coords. It falls back to
    a position-overlap join on the parent assembly contig only when no
    direct match exists.
    """
    df = pl.read_csv(path, separator="\t", infer_schema_length=10_000,
                     null_values=["NA", ""])
    rows: list[dict] = []
    has_hallmark = "virus_hallmark" in df.columns
    has_marker = "marker" in df.columns
    for row in df.iter_rows(named=True):
        gene = row.get("gene")
        if not gene:
            continue
        gene = str(gene).strip()
        if not gene:
            continue
        category: str | None = None
        if has_hallmark and _is_truthy(row.get("virus_hallmark")):
            category = "0"
        elif has_marker:
            m = row.get("marker")
            if isinstance(m, str):
                m = m.strip()
                if m and m != "NA":
                    suffix = m.rsplit(".", 1)[-1] if "." in m else ""
                    if suffix in _GENOMAD_HALLMARK_SUFFIXES:
                        category = "0"
                    elif suffix in _GENOMAD_VIRAL_LIKE_SUFFIXES:
                        category = "1"
        if category is None:
            continue
        try:
            start = int(row.get("start"))
            end = int(row.get("end"))
        except (TypeError, ValueError):
            continue
        rows.append({
            "gene": gene,
            "contig": _parse_genomad_parent_contig(gene),
            "start": start,
            "end": end,
            "category": category,
        })
    return pl.DataFrame(
        rows,
        schema={"gene": pl.Utf8, "contig": pl.Utf8,
                "start": pl.Int64, "end": pl.Int64,
                "category": pl.Utf8},
    )


def build_virsorter_categories_from_genomad(
    genomad_df: pl.DataFrame,
    annotations: pl.DataFrame,
) -> dict[str, str]:
    """Map geNomad's gene-level categories onto DRAM's gene ids.

    Two-step join:

      1. **Direct gene-id equality.** geNomad and DRAM both name provirus
         contigs `<parent>|provirus_<start>_<end>` and number genes 1-based
         within the provirus, so a gene id like
         `k141_100971|provirus_1_18064_4` is identical in both tools.
         This is the common case and catches whole-virus contigs too
         (`k141_85503_1` etc).
      2. **Position-overlap on parent assembly contig** (fallback). Only
         consulted for DRAM genes that didn't get a direct match. Useful
         when DRAM and geNomad disagree on contig naming — e.g. CheckV
         further trimmed a geNomad-named provirus, or DRAM was run on a
         clustered catalog whose ids don't match per-sample geNomad. For
         each unmatched DRAM gene, parse its scaffold to extract the
         parent contig + start offset (offset 1 for whole-virus), then
         translate provirus coords back to parent coords and look for a
         geNomad gene on the same parent whose interval overlaps.
    """
    if genomad_df.is_empty():
        return {}

    # Direct gene-id index
    by_gene_id = dict(zip(
        genomad_df["gene"].to_list(),
        genomad_df["category"].to_list(),
    ))

    # Position-overlap fallback index
    by_contig: dict[str, list[tuple[int, int, str]]] = {}
    for row in genomad_df.iter_rows(named=True):
        by_contig.setdefault(row["contig"], []).append(
            (int(row["start"]), int(row["end"]), str(row["category"]))
        )

    out: dict[str, str] = {}
    needed = annotations.select(
        ["query_id", "scaffold", "start_position", "stop_position"]
    ).unique()
    for row in needed.iter_rows(named=True):
        qid = str(row["query_id"])
        # 1. Direct match
        cat = by_gene_id.get(qid)
        if cat is not None:
            out[qid] = cat
            continue
        # 2. Position-overlap fallback on parent assembly contig
        parent, offset = _parse_dram_scaffold(str(row["scaffold"]))
        hits = by_contig.get(parent)
        if not hits:
            continue
        orig_start = offset + int(row["start_position"]) - 1
        orig_end = offset + int(row["stop_position"]) - 1
        for g_start, g_end, c in hits:
            if g_start <= orig_end and g_end >= orig_start:
                out[qid] = c
                break
    return out


def _aux_score_for_scaffold(
    gene_ids: list[str],
    virsorter_categories: dict[str, str],
) -> dict[str, int]:
    """Verbatim port of v1 calculate_auxiliary_scores for one scaffold.

    `gene_ids` is the scaffold's genes in start-position order. v1 builds
    a (dram_gene, virsorter_gene, virsorter_category) tuple per gene and
    iterates with index `i`. In pure-FASTA mode every virsorter_gene is
    None, so the left/right context lists are empty and every gene falls
    through to the default score of 5. With VirSorter input, neighboring
    hallmark/viral-like calls drop the score per the v1 if-else chain.
    """
    n = len(gene_ids)
    scores: dict[str, int] = {}
    for i, gene_id in enumerate(gene_ids):
        if i == 0 or i == n - 1:
            scores[gene_id] = 5
            continue

        left_cats = {
            virsorter_categories[g]
            for g in gene_ids[:i] if g in virsorter_categories
        }
        right_cats = {
            virsorter_categories[g]
            for g in gene_ids[i + 1:] if g in virsorter_categories
        }
        hallmark_left = bool(left_cats & _VIRSORTER_HALLMARK_CATEGORIES)
        viral_like_left = bool(left_cats & _VIRSORTER_VIRAL_LIKE_CATEGORIES)
        hallmark_right = bool(right_cats & _VIRSORTER_HALLMARK_CATEGORIES)
        viral_like_right = bool(right_cats & _VIRSORTER_VIRAL_LIKE_CATEGORIES)

        own_cat = virsorter_categories.get(gene_id)

        if hallmark_left and hallmark_right:
            score = 1
        elif (hallmark_left and viral_like_right) or (viral_like_left and hallmark_right):
            score = 2
        elif viral_like_left and viral_like_right:
            score = 3
        elif hallmark_left or viral_like_left or hallmark_right or viral_like_right:
            score = 4
        elif own_cat in _VIRSORTER_HALLMARK_CATEGORIES or own_cat in _VIRSORTER_VIRAL_LIKE_CATEGORIES:
            score = 4
        else:
            score = 5
        scores[gene_id] = score
    return scores


def explode_gene_ids(annotations: pl.DataFrame) -> pl.DataFrame:
    """Add a `_gene_ids` list[str] column: union of ids across all besthit columns."""
    cols = [c for c in _ID_EXPR_DICT if c in annotations.columns]
    if not cols:
        return annotations.with_columns(
            pl.Series(
                "_gene_ids",
                [[] for _ in range(annotations.height)],
                dtype=pl.List(pl.Utf8),
            )
        )
    tmp_aliases = [f"_parsed_{c}" for c in cols]
    parsed = annotations.with_columns(
        [_ID_EXPR_DICT[col].alias(tmp) for col, tmp in zip(cols, tmp_aliases)]
    )
    parsed = parsed.with_columns(
        pl.concat_list([pl.col(t).fill_null([]) for t in tmp_aliases])
          .list.unique()
          .alias("_gene_ids")
    )
    return parsed.drop(tmp_aliases)


def compute_flags(
    annotations: pl.DataFrame,
    metabolic_genes: set[str],
    amgs: set[str],
    verified_amgs: set[str],
    scaffold_lengths: dict[str, int],
    length_from_end: int,
    essential_amgs: set[str] | None = None,
    viral_vog_ids: set[str] | None = None,
    virsorter_categories: dict[str, str] | None = None,
) -> pl.DataFrame:
    if essential_amgs is None:
        essential_amgs = set()
    if viral_vog_ids is None:
        viral_vog_ids = set()
    if virsorter_categories is None:
        virsorter_categories = {}
    df = explode_gene_ids(annotations)

    def _has_intersect(set_: set[str]) -> pl.Expr:
        if not set_:
            return pl.lit(False)
        return (
            pl.col("_gene_ids")
              .list.eval(pl.element().is_in(list(set_)))
              .list.any()
              .fill_null(False)
        )

    df = df.with_columns([
        _has_intersect(metabolic_genes).alias("_M"),
        _has_intersect(amgs).alias("_K"),
        _has_intersect(verified_amgs).alias("_E"),
        _has_intersect(viral_vog_ids).alias("_V"),
        _has_intersect(CELL_ENTRY_CAZYS).alias("_A"),
        _has_intersect(VIRAL_PEPTIDASES_MEROPS).alias("_P"),
        _has_intersect(TRANSPOSON_PFAMS).alias("is_transposon"),
        _has_intersect(essential_amgs).alias("_N"),
    ])

    # T flag: any gene on the same scaffold has is_transposon=True.
    df = df.with_columns(
        pl.col("is_transposon").any().over("scaffold").alias("_T")
    )

    # F flag: gene is within length_from_end bp of either contig end.
    if scaffold_lengths:
        length_df = pl.DataFrame({
            "scaffold": list(scaffold_lengths.keys()),
            "_scaffold_length": list(scaffold_lengths.values()),
        })
        df = df.join(length_df, on="scaffold", how="left")
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Int64).alias("_scaffold_length"))
    df = df.with_columns(
        (
            (pl.col("start_position") < length_from_end)
            | (pl.col("stop_position") > pl.col("_scaffold_length") - length_from_end)
        ).fill_null(False).alias("_F")
    )

    # B flag: any 3-consecutive-M-gene window on the same scaffold (sorted by start).
    # v1 sets B on all three genes of every (prev, self, next) triple where all three
    # carry M. Equivalently: gene X gets B iff at least one of the three rolling
    # 3-windows that include X (ending at X, X+1, or X+2) sums to 3 on _M.
    df = df.sort(["scaffold", "start_position"])
    m_int = pl.col("_M").cast(pl.Int32)
    roll = m_int.rolling_sum(window_size=3, min_samples=3).over("scaffold")
    df = df.with_columns([
        roll.alias("_roll_here"),
        roll.shift(-1).over("scaffold").alias("_roll_next"),
        roll.shift(-2).over("scaffold").alias("_roll_after"),
    ])
    df = df.with_columns(
        (
            (pl.col("_roll_here") == 3)
            | (pl.col("_roll_next") == 3)
            | (pl.col("_roll_after") == 3)
        ).fill_null(False).alias("_B")
    )

    # K forces M (per v1).
    df = df.with_columns((pl.col("_M") | pl.col("_K")).alias("_M"))

    # Build amg_flags string in v1 order (M K E V A P T F B), with non-v1 N
    # appended at the end.
    flag_str = (
        pl.when(pl.col("_M")).then(pl.lit("M")).otherwise(pl.lit(""))
        + pl.when(pl.col("_K")).then(pl.lit("K")).otherwise(pl.lit(""))
        + pl.when(pl.col("_E")).then(pl.lit("E")).otherwise(pl.lit(""))
        + pl.when(pl.col("_V")).then(pl.lit("V")).otherwise(pl.lit(""))
        + pl.when(pl.col("_A")).then(pl.lit("A")).otherwise(pl.lit(""))
        + pl.when(pl.col("_P")).then(pl.lit("P")).otherwise(pl.lit(""))
        + pl.when(pl.col("_T")).then(pl.lit("T")).otherwise(pl.lit(""))
        + pl.when(pl.col("_F")).then(pl.lit("F")).otherwise(pl.lit(""))
        + pl.when(pl.col("_B")).then(pl.lit("B")).otherwise(pl.lit(""))
        + pl.when(pl.col("_N")).then(pl.lit("N")).otherwise(pl.lit(""))
    ).alias("amg_flags")
    df = df.with_columns(flag_str)

    # auxiliary_score: per-scaffold v1 algorithm, then B-flag downgrade.
    scores: dict[str, int] = {}
    for scaffold_name, group in df.group_by("scaffold", maintain_order=True):
        gene_ids = group["query_id"].to_list()
        scores.update(_aux_score_for_scaffold(gene_ids, virsorter_categories))
    score_df = pl.DataFrame({
        "query_id": list(scores.keys()),
        "auxiliary_score": list(scores.values()),
    }, schema_overrides={"auxiliary_score": pl.Int64})
    df = df.join(score_df, on="query_id", how="left")
    # v1 downgrade: B-flagged genes with score < 4 are bumped up to 4.
    df = df.with_columns(
        pl.when(pl.col("_B") & (pl.col("auxiliary_score") < 4))
          .then(4)
          .otherwise(pl.col("auxiliary_score"))
          .alias("auxiliary_score")
    )

    drop_cols = [c for c in df.columns if c.startswith("_")]
    return df.drop(drop_cols)


@click.command()
@click.option("-i", "--input_file", required=True, type=click.Path(exists=True),
              help="Combined annotations TSV (output of COMBINE_ANNOTATIONS).")
@click.option("-o", "--output_file", required=True, type=click.Path(),
              help="Path to write annotations with amg_flags + is_transposon columns.")
@click.option("--catalog_fasta", required=True, type=click.Path(exists=True),
              help="Multi-FASTA of viral contigs the annotations were derived from.")
@click.option("--amg_db", default=None, type=click.Path(),
              help="AMG reference TSV. Defaults to the bundled assets/amg_database.tsv.")
@click.option("--distill_sheets_dir", default=None, type=click.Path(),
              help="Distill TSV directory. Defaults to bundled assets/forms/distill_sheets.")
@click.option("--vog_list", default=None, type=click.Path(),
              help="VOGdb vog_annotations_latest.tsv(.gz). Required for the V flag; "
                   "if omitted the V flag is never set.")
@click.option("--genomad_genes", multiple=True, type=click.Path(),
              help="geNomad *_genes.tsv. Drives auxiliary_score: virus_hallmark "
                   "→ '0', taxname starting 'Viruses' → '1'. Pass once per "
                   "sample (the flag may be repeated) — entries are merged into "
                   "a single {gene_id: category} dict. If omitted every gene "
                   "gets the v1 fallback score 5.")
@click.option("--length_from_end", default=DEFAULT_LENGTH_FROM_END, type=int,
              show_default=True,
              help="Window (bp) from contig ends used to set the F flag.")
def main(input_file, output_file, catalog_fasta, amg_db, distill_sheets_dir,
         vog_list, genomad_genes, length_from_end):
    """Add DRAM-v amg_flags and is_transposon columns to a combined annotations TSV."""
    here = Path(__file__).parent
    amg_db = Path(amg_db) if amg_db else here / "assets" / "amg_database.tsv"
    distill_sheets_dir = (
        Path(distill_sheets_dir) if distill_sheets_dir
        else here / "assets" / "forms" / "distill_sheets"
    )

    logger.info(f"Reading annotations: {input_file}")
    annotations = pl.read_csv(input_file, separator="\t", infer_schema_length=10_000)

    logger.info(f"Reading AMG database: {amg_db}")
    amgs, verified_amgs, essential_amgs = build_amg_id_sets(amg_db)
    logger.info(
        f"AMG ids: {len(amgs)} ({len(verified_amgs)} verified, "
        f"{len(essential_amgs)} essential viral function)"
    )

    logger.info(f"Building metabolic_genes set from {distill_sheets_dir}")
    metabolic_genes = build_metabolic_genes(distill_sheets_dir)
    logger.info(f"Metabolic genes: {len(metabolic_genes)}")

    viral_vog_ids: set[str] = set()
    if vog_list:
        logger.info(f"Building viral VOG id set (Xr/Xs) from {vog_list}")
        viral_vog_ids = build_viral_vog_ids(Path(vog_list))
        logger.info(f"Viral VOG ids: {len(viral_vog_ids)}")
    else:
        logger.info("No --vog_list provided; V flag will not be set.")

    virsorter_categories: dict[str, str] = {}
    if genomad_genes:
        parts: list[pl.DataFrame] = []
        for p in genomad_genes:
            logger.info(f"Parsing geNomad genes TSV for VirSorter category mapping: {p}")
            parts.append(parse_genomad_genes_tsv(Path(p)))
        non_empty = [df for df in parts if not df.is_empty()]
        genomad_df = pl.concat(non_empty) if non_empty else pl.DataFrame(
            schema={"contig": pl.Utf8, "start": pl.Int64, "end": pl.Int64,
                    "category": pl.Utf8})
        n_hallmark = int((genomad_df["category"] == "0").sum()) if not genomad_df.is_empty() else 0
        n_viral_like = int((genomad_df["category"] == "1").sum()) if not genomad_df.is_empty() else 0
        logger.info(
            f"geNomad-derived categories across {len(genomad_genes)} file(s): "
            f"{n_hallmark} hallmark + {n_viral_like} viral-like"
        )
        virsorter_categories = build_virsorter_categories_from_genomad(
            genomad_df, annotations
        )
        logger.info(
            f"DRAM genes with mapped category (after provirus position-join): "
            f"{len(virsorter_categories)} of {annotations.height}"
        )
    else:
        logger.info("No --genomad_genes provided; auxiliary_score will fall through to 5 for every gene.")

    logger.info(f"Reading scaffold lengths from {catalog_fasta}")
    scaffold_lengths = read_scaffold_lengths(catalog_fasta)
    logger.info(f"Scaffolds: {len(scaffold_lengths)}")

    annotated = compute_flags(
        annotations, metabolic_genes, amgs, verified_amgs,
        scaffold_lengths, length_from_end,
        essential_amgs=essential_amgs,
        viral_vog_ids=viral_vog_ids,
        virsorter_categories=virsorter_categories,
    )
    logger.info(f"Writing {output_file}")
    annotated.write_csv(output_file, separator="\t")


if __name__ == "__main__":
    main()
