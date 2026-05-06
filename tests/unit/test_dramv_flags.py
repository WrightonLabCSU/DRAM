"""Unit tests for bin/dramv_flags.py compute_flags + helpers.

Run with:  pytest tests/unit/test_dramv_flags.py
"""

import polars as pl
import pytest

from dramv_flags import compute_flags, read_scaffold_lengths


def _ann(rows):
    """Build a small annotations DataFrame with the columns dramv_flags reads."""
    return pl.DataFrame(rows, schema_overrides={
        "kofam_id": pl.Utf8,
        "pfam_hits": pl.Utf8,
        "dbcan_id": pl.Utf8,
    })


def _flags_by_query(df: pl.DataFrame) -> dict[str, str]:
    return {row["query_id"]: row["amg_flags"] for row in df.iter_rows(named=True)}


def test_individual_flags_fire():
    """One scaffold, five genes, each carrying a distinct flag-trigger."""
    ann = _ann([
        # K00001 metabolic, near-start (F), part of M-triple (B)
        {"query_id": "g1", "scaffold": "s1", "start_position": 100, "stop_position": 500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
        # K00002 metabolic, mid-scaffold (no F), part of M-triple (B)
        {"query_id": "g2", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K00002", "pfam_hits": None, "dbcan_id": None},
        # K00003 metabolic, part of M-triple (B), near-end (F)
        {"query_id": "g3", "scaffold": "s1", "start_position": 7000, "stop_position": 7500,
         "kofam_id": "K00003", "pfam_hits": None, "dbcan_id": None},
        # GH18 → A flag (CELL_ENTRY_CAZYS); no M
        {"query_id": "g4", "scaffold": "s1", "start_position": 8000, "stop_position": 8500,
         "kofam_id": None, "pfam_hits": None, "dbcan_id": "GH18"},
        # K99999 in amgs (not metabolic) → K forces M; PF01609 → is_transposon → T fires for all on s1
        {"query_id": "g5", "scaffold": "s1", "start_position": 9000, "stop_position": 9500,
         "kofam_id": "K99999", "pfam_hits": "[PF01609.1]", "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001", "K00002", "K00003"},
        amgs={"K99999"},
        verified_amgs=set(),
        scaffold_lengths={"s1": 10_000},
        length_from_end=5_000,
    )
    flags = _flags_by_query(out)

    assert flags["g1"] == "MTFB"
    assert flags["g2"] == "MTFB"   # 6000–6500 falls within F window (10000-5000=5000)
    assert flags["g3"] == "MTFB"
    assert flags["g4"] == "ATF"    # A only; T scaffold-wide; F because near end
    assert flags["g5"] == "MKTF"   # K forces M; T from own pfam; F near end; no B (g4 broke triple)

    # is_transposon must only be True for the row carrying the trigger pfam.
    transposons = {row["query_id"]: row["is_transposon"] for row in out.iter_rows(named=True)}
    assert transposons == {"g1": False, "g2": False, "g3": False, "g4": False, "g5": True}


def test_b_flag_does_not_cross_scaffold_boundaries():
    """Two scaffolds, B must only fire within the same scaffold's M-triple."""
    rows = [
        # s1: 4 M genes in a row → all four get B
        {"query_id": "s1g1", "scaffold": "s1", "start_position": 100,  "stop_position": 500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "s1g2", "scaffold": "s1", "start_position": 2000, "stop_position": 2500,
         "kofam_id": "K00002", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "s1g3", "scaffold": "s1", "start_position": 4000, "stop_position": 4500,
         "kofam_id": "K00003", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "s1g4", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K00004", "pfam_hits": None, "dbcan_id": None},
        # s2: a single M gene, must not inherit B from s1's triple
        {"query_id": "s2g1", "scaffold": "s2", "start_position": 100, "stop_position": 500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
    ]
    out = compute_flags(
        _ann(rows),
        metabolic_genes={"K00001", "K00002", "K00003", "K00004"},
        amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 8_000, "s2": 3_000},
        length_from_end=5_000,
    )
    flags = _flags_by_query(out)
    for g in ("s1g1", "s1g2", "s1g3", "s1g4"):
        assert "B" in flags[g], f"{g}: expected B, got {flags[g]!r}"
    assert "B" not in flags["s2g1"]


def test_sub_three_gene_scaffolds_get_no_b():
    """v1 only set B inside iter through middle indices; scaffolds with fewer
    than 3 genes can never form an M-triple."""
    rows = [
        # s1: 2 M genes — no B possible
        {"query_id": "s1g1", "scaffold": "s1", "start_position": 100,  "stop_position": 500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "s1g2", "scaffold": "s1", "start_position": 2000, "stop_position": 2500,
         "kofam_id": "K00002", "pfam_hits": None, "dbcan_id": None},
        # s2: 1 M gene alone — no B
        {"query_id": "s2g1", "scaffold": "s2", "start_position": 100, "stop_position": 500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
    ]
    out = compute_flags(
        _ann(rows),
        metabolic_genes={"K00001", "K00002"},
        amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 8_000, "s2": 3_000},
        length_from_end=5_000,
    )
    flags = _flags_by_query(out)
    for g in ("s1g1", "s1g2", "s2g1"):
        assert "B" not in flags[g], f"{g}: should not have B, got {flags[g]!r}"


def test_k_flag_forces_m():
    """A row whose only AMG-relevant id is in `amgs` (not in metabolic_genes)
    should still get M because v1's K rule force-sets M."""
    ann = _ann([
        {"query_id": "g1", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K99999", "pfam_hits": None, "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes=set(),     # deliberately empty
        amgs={"K99999"},
        verified_amgs=set(),
        scaffold_lengths={"s1": 10_000},
        length_from_end=5_000,
    )
    flags = _flags_by_query(out)
    assert flags["g1"].startswith("MK"), f"expected M before K (v1 order), got {flags['g1']!r}"


def test_e_flag_for_verified_amgs():
    """A verified AMG mid-scaffold (no F) should produce exactly MKE."""
    ann = _ann([
        {"query_id": "g1", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001"},
        amgs={"K00001"},
        verified_amgs={"K00001"},
        scaffold_lengths={"s1": 14_000},  # F-free window is 5000..9000
        length_from_end=5_000,
    )
    assert _flags_by_query(out)["g1"] == "MKE"


def test_f_flag_window():
    """F triggers iff start < length_from_end OR stop > scaffold_length - length_from_end."""
    ann = _ann([
        {"query_id": "near_start", "scaffold": "s1", "start_position": 100,  "stop_position": 500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "middle",     "scaffold": "s1", "start_position": 5500, "stop_position": 5800,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "near_end",   "scaffold": "s1", "start_position": 9500, "stop_position": 9900,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001"},
        amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 11_000},  # F window covers <5000 or >6000
        length_from_end=5_000,
    )
    flags = _flags_by_query(out)
    assert "F" in flags["near_start"]
    assert "F" not in flags["middle"]
    assert "F" in flags["near_end"]


def test_n_flag_fires_for_essential_viral_function():
    """A row whose only AMG-relevant id sits in essential_amgs gets N at the
    end of amg_flags. N is informational — it must NOT force M."""
    ann = _ann([
        # PF06508 = QueC, paper-flagged essential viral function.
        {"query_id": "g1", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": None, "pfam_hits": "[PF06508.1]", "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes=set(), amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 14_000},  # F-free window 5000..9000
        length_from_end=5_000,
        essential_amgs={"PF06508"},
    )
    flags = _flags_by_query(out)
    # Just N, with no M (no force from N), no K, no E, no F.
    assert flags["g1"] == "N", f"expected exactly 'N', got {flags['g1']!r}"


def test_n_flag_combines_with_k_and_e():
    """Same id can be in amgs, verified_amgs, and essential_amgs simultaneously
    (e.g. mazG: verified=TRUE and essential_viral_function=TRUE in the bundled
    db). Resulting flag string must include K, E, and trailing N, with M from
    the v1 K-forces-M rule."""
    ann = _ann([
        {"query_id": "g1", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K04765", "pfam_hits": None, "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes=set(),
        amgs={"K04765"},
        verified_amgs={"K04765"},
        scaffold_lengths={"s1": 14_000},
        length_from_end=5_000,
        essential_amgs={"K04765"},
    )
    assert _flags_by_query(out)["g1"] == "MKEN"


def test_n_flag_default_is_off_when_no_essential_set_passed():
    """Backward-compat: the existing call signature (no essential_amgs kwarg)
    still works, and produces no N in the output."""
    ann = _ann([
        {"query_id": "g1", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001"},
        amgs={"K00001"},
        verified_amgs={"K00001"},
        scaffold_lengths={"s1": 14_000},
        length_from_end=5_000,
        # essential_amgs intentionally omitted
    )
    assert _flags_by_query(out)["g1"] == "MKE"


def test_build_amg_id_sets_returns_essential(tmp_path):
    """build_amg_id_sets parses the new essential_viral_function column."""
    from dramv_flags import build_amg_id_sets
    db = tmp_path / "amg.tsv"
    db.write_text(
        "KO\tEC\tPFAM\tgene\tmodule\tmetabolism\treference\tverified\tessential_viral_function\n"
        "K00001\t\t\tg1\t\t\tref\tTRUE\tFALSE\n"
        "K00002\t\tPF99999\tg2\t\t\tref\tFALSE\tTRUE\n"
        "K00003\t\t\tg3\t\t\tref\tFALSE\tFALSE\n"
    )
    amgs, verified, essential = build_amg_id_sets(db)
    assert amgs == {"K00001", "K00002", "K00003", "PF99999"}
    assert verified == {"K00001"}
    assert essential == {"K00002", "PF99999"}


def test_build_amg_id_sets_tolerates_missing_essential_column(tmp_path):
    """An older amg_database.tsv (no essential_viral_function column) must still
    load; essential set comes back empty."""
    from dramv_flags import build_amg_id_sets
    db = tmp_path / "amg.tsv"
    db.write_text(
        "KO\tEC\tPFAM\tgene\tmodule\tmetabolism\treference\tverified\n"
        "K00001\t\t\tg1\t\t\tref\tTRUE\n"
    )
    amgs, verified, essential = build_amg_id_sets(db)
    assert amgs == {"K00001"}
    assert verified == {"K00001"}
    assert essential == set()


def test_read_scaffold_lengths(tmp_path):
    """Scaffold lengths are read correctly across multi-line records, including
    trailing newlines and IDs with embedded spaces."""
    p = tmp_path / "tiny.fa"
    p.write_text(">contig_a some description\nACGT\nACGT\n>contig_b\nNNN\n")
    assert read_scaffold_lengths(str(p)) == {"contig_a": 8, "contig_b": 3}
