#!/usr/bin/env python
"""Compute DRAM-v AMG flags and is_transposon column on a combined annotations TSV.

Phase 1 (FASTA-only): no VirSorter input required. Adds two columns:

  - is_transposon (bool): the gene's Pfam hits intersect TRANSPOSON_PFAMS.
  - amg_flags (str): concatenated single-letter flags in v1 order
    M / K / E / A / P / T / F / B, then a trailing N (not in v1).
    Semantics for M-B ported verbatim from DRAM v1
    mag_annotator/annotate_vgfs.py:get_metabolic_flags (commit 6cd68f9),
    minus the V (VOGdb) flag and the auxiliary_score block, both deferred
    to Phase 2. The N flag is informational — it marks genes whose ids
    appear in amg_database.tsv rows where essential_viral_function=TRUE,
    per Martin et al. 2025 (doi:10.1038/s41564-025-02095-4). It does NOT
    propagate to M and is NOT used by any downstream filter yet.
"""

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
# vog_id / vog_ids aren't yet in the upstream rule_parser ID_EXPR_DICT, so
# register them here. vog_id is the bare best-hit name (e.g. "VOG00177") and
# vog_ids is "; "-joined all hits per hmm_parser. extract_all on /VOG\d+/
# tolerates both shapes.
_ID_EXPR_DICT["vog_id"] = pl.col("vog_id").cast(pl.Utf8).str.extract_all(_VOG_ID_RE)
_ID_EXPR_DICT["vog_ids"] = pl.col("vog_ids").cast(pl.Utf8).str.extract_all(_VOG_ID_RE)


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
) -> pl.DataFrame:
    if essential_amgs is None:
        essential_amgs = set()
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

    # Build amg_flags string in v1 order, with non-v1 N appended at the end.
    flag_str = (
        pl.when(pl.col("_M")).then(pl.lit("M")).otherwise(pl.lit(""))
        + pl.when(pl.col("_K")).then(pl.lit("K")).otherwise(pl.lit(""))
        + pl.when(pl.col("_E")).then(pl.lit("E")).otherwise(pl.lit(""))
        + pl.when(pl.col("_A")).then(pl.lit("A")).otherwise(pl.lit(""))
        + pl.when(pl.col("_P")).then(pl.lit("P")).otherwise(pl.lit(""))
        + pl.when(pl.col("_T")).then(pl.lit("T")).otherwise(pl.lit(""))
        + pl.when(pl.col("_F")).then(pl.lit("F")).otherwise(pl.lit(""))
        + pl.when(pl.col("_B")).then(pl.lit("B")).otherwise(pl.lit(""))
        + pl.when(pl.col("_N")).then(pl.lit("N")).otherwise(pl.lit(""))
    ).alias("amg_flags")
    df = df.with_columns(flag_str)

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
@click.option("--length_from_end", default=DEFAULT_LENGTH_FROM_END, type=int,
              show_default=True,
              help="Window (bp) from contig ends used to set the F flag.")
def main(input_file, output_file, catalog_fasta, amg_db, distill_sheets_dir, length_from_end):
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

    logger.info(f"Reading scaffold lengths from {catalog_fasta}")
    scaffold_lengths = read_scaffold_lengths(catalog_fasta)
    logger.info(f"Scaffolds: {len(scaffold_lengths)}")

    annotated = compute_flags(
        annotations, metabolic_genes, amgs, verified_amgs,
        scaffold_lengths, length_from_end,
        essential_amgs=essential_amgs,
    )
    logger.info(f"Writing {output_file}")
    annotated.write_csv(output_file, separator="\t")


if __name__ == "__main__":
    main()
