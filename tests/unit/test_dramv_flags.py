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


def test_v_flag_fires_for_xr_xs_vog_hits():
    """V flag fires when the gene's vogdb_ids contain at least one VOG whose
    FunctionalCategory is Xr or Xs. Inserted between E and A in v1 order."""
    ann = _ann([
        # Xr hit only → V at v1 position (between E and A)
        {"query_id": "g_xr", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": None, "pfam_hits": None, "dbcan_id": None,
         "vogdb_id": "VOG00001", "vogdb_ids": "VOG00001"},
        # Xs in a multi-id list (joined by "; ") still fires V
        {"query_id": "g_xs", "scaffold": "s1", "start_position": 6700, "stop_position": 7000,
         "kofam_id": None, "pfam_hits": None, "dbcan_id": None,
         "vogdb_id": "VOG00099", "vogdb_ids": "VOG00099; VOG00002"},
        # Xu (unknown) — must NOT fire V
        {"query_id": "g_xu", "scaffold": "s1", "start_position": 7100, "stop_position": 7400,
         "kofam_id": None, "pfam_hits": None, "dbcan_id": None,
         "vogdb_id": "VOG00099", "vogdb_ids": "VOG00099"},
    ])
    out = compute_flags(
        ann,
        metabolic_genes=set(), amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 14_000}, length_from_end=5_000,
        viral_vog_ids={"VOG00001", "VOG00002"},
    )
    flags = _flags_by_query(out)
    assert flags["g_xr"] == "V"
    assert flags["g_xs"] == "V"
    assert flags["g_xu"] == ""


def test_v_flag_combines_with_other_flags_in_v1_order():
    """V slots between E and A. A gene with M+K+E+V+A should produce 'MKEVA'."""
    ann = _ann([
        # K00001 metabolic + verified AMG + GH18 (A) + VOG hit Xr
        {"query_id": "g1", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": "GH18",
         "vogdb_id": "VOG00001", "vogdb_ids": "VOG00001"},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001"},
        amgs={"K00001"}, verified_amgs={"K00001"},
        scaffold_lengths={"s1": 14_000}, length_from_end=5_000,
        viral_vog_ids={"VOG00001"},
    )
    assert _flags_by_query(out)["g1"] == "MKEVA"


def test_v_flag_default_off_when_no_viral_vog_set_passed():
    """Backward compat: omitting viral_vog_ids never produces V even if a row
    has a vogdb_ids hit."""
    ann = _ann([
        {"query_id": "g1", "scaffold": "s1", "start_position": 6000, "stop_position": 6500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None,
         "vogdb_id": "VOG00001", "vogdb_ids": "VOG00001"},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001"},
        amgs={"K00001"}, verified_amgs={"K00001"},
        scaffold_lengths={"s1": 14_000}, length_from_end=5_000,
        # viral_vog_ids intentionally omitted
    )
    assert "V" not in _flags_by_query(out)["g1"]


def test_build_viral_vog_ids(tmp_path):
    """build_viral_vog_ids picks up GroupNames whose FunctionalCategory contains
    Xr or Xs, and rejects Xh/Xp/Xu."""
    from dramv_flags import build_viral_vog_ids
    p = tmp_path / "vog_annot.tsv"
    p.write_text(
        "GroupName\tProteinCount\tSpeciesCount\tFunctionalCategory\tConsensusFunctionalDescription\n"
        "VOG00001\t10\t5\tXr\tviral replication thing\n"
        "VOG00002\t10\t5\tXs\tvirion structural\n"
        "VOG00003\t10\t5\tXrXs\tcombined replication+structure\n"
        "VOG00004\t10\t5\tXh\thost benefit\n"
        "VOG00005\t10\t5\tXp\tintegration\n"
        "VOG00006\t10\t5\tXu\tunknown\n"
        "VOG00007\t10\t5\t\tno category\n"
    )
    assert build_viral_vog_ids(p) == {"VOG00001", "VOG00002", "VOG00003"}


def test_auxiliary_score_pure_fasta_mode_is_all_fives():
    """Without virsorter_categories, every gene gets the v1 fallback score 5."""
    ann = _ann([
        {"query_id": "g1", "scaffold": "s1", "start_position": 100,  "stop_position": 500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "g2", "scaffold": "s1", "start_position": 1000, "stop_position": 1500,
         "kofam_id": None,    "pfam_hits": None, "dbcan_id": None},
        {"query_id": "g3", "scaffold": "s1", "start_position": 2000, "stop_position": 2500,
         "kofam_id": None,    "pfam_hits": None, "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001"},
        amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 10_000}, length_from_end=5_000,
    )
    by_id = {row["query_id"]: row["auxiliary_score"] for row in out.iter_rows(named=True)}
    assert by_id == {"g1": 5, "g2": 5, "g3": 5}


def test_auxiliary_score_v1_if_chain_with_virsorter_categories():
    """Reproduce each branch of v1's calculate_auxiliary_scores using a 5-gene
    scaffold and supplied virsorter_categories. Indices 0 and 4 are ends → 5."""
    ann = _ann([
        {"query_id": f"g{i}", "scaffold": "s1", "start_position": i * 1000,
         "stop_position": i * 1000 + 500,
         "kofam_id": None, "pfam_hits": None, "dbcan_id": None}
        for i in range(5)
    ])
    out = compute_flags(
        ann,
        metabolic_genes=set(), amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 100_000},  # F-window stays well off
        length_from_end=5_000,
        virsorter_categories={"g0": "0", "g1": "0", "g3": "0", "g4": "0"},  # hallmark on both flanks of g2
    )
    by_id = {row["query_id"]: row["auxiliary_score"] for row in out.iter_rows(named=True)}
    assert by_id["g0"] == 5    # first
    assert by_id["g4"] == 5    # last
    assert by_id["g2"] == 1    # hallmark left + hallmark right


def test_auxiliary_score_branches():
    """Each non-end-gene branch of the v1 if-chain on its own 3-gene scaffold."""
    cases = [
        # (scaffold, vs_categories_for_neighbors, own_cat, expected_score)
        ("s_hh", {0: "0",  2: "0"},  None, 1),  # hallmark left + hallmark right
        ("s_hv", {0: "0",  2: "1"},  None, 2),  # hallmark left + viral_like right
        ("s_vv", {0: "1",  2: "1"},  None, 3),  # viral_like both sides
        ("s_hl", {0: "0"},           None, 4),  # hallmark left only, nothing right
        ("s_sh", {},                 "0",  4),  # own=hallmark, no neighbor cats
        ("s_no", {},                 None, 5),  # nothing anywhere
    ]
    rows = []
    vs_all: dict[str, str] = {}
    for sc, neighbor_cats, own_cat, _ in cases:
        for i in range(3):
            qid = f"q_{sc}_{i}"
            rows.append({
                "query_id": qid, "scaffold": sc,
                "start_position": i * 1000, "stop_position": i * 1000 + 500,
                "kofam_id": None, "pfam_hits": None, "dbcan_id": None,
            })
            if i in neighbor_cats:
                vs_all[qid] = neighbor_cats[i]
            if i == 1 and own_cat is not None:
                vs_all[qid] = own_cat
    out = compute_flags(
        _ann(rows),
        metabolic_genes=set(), amgs=set(), verified_amgs=set(),
        scaffold_lengths={s: 100_000 for s, _, _, _ in cases},
        length_from_end=5_000,
        virsorter_categories=vs_all,
    )
    by_id = {row["query_id"]: row["auxiliary_score"] for row in out.iter_rows(named=True)}
    for sc, _, _, expected in cases:
        assert by_id[f"q_{sc}_1"] == expected, \
            f"{sc}: expected {expected}, got {by_id[f'q_{sc}_1']}"


def test_b_flag_downgrades_low_auxiliary_score_to_four():
    """v1: a gene with B in amg_flags AND auxiliary_score < 4 is bumped to 4."""
    # 3 metabolic genes back-to-back → all get B. Place hallmark VirSorter
    # neighbors so g2 (the middle) would otherwise score 1.
    ann = _ann([
        {"query_id": "g0", "scaffold": "s1", "start_position": 1000,  "stop_position": 1500,
         "kofam_id": None,    "pfam_hits": None, "dbcan_id": None},
        {"query_id": "g1", "scaffold": "s1", "start_position": 2000,  "stop_position": 2500,
         "kofam_id": "K00001", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "g2", "scaffold": "s1", "start_position": 3000,  "stop_position": 3500,
         "kofam_id": "K00002", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "g3", "scaffold": "s1", "start_position": 4000,  "stop_position": 4500,
         "kofam_id": "K00003", "pfam_hits": None, "dbcan_id": None},
        {"query_id": "g4", "scaffold": "s1", "start_position": 5000,  "stop_position": 5500,
         "kofam_id": None,    "pfam_hits": None, "dbcan_id": None},
    ])
    out = compute_flags(
        ann,
        metabolic_genes={"K00001", "K00002", "K00003"},
        amgs=set(), verified_amgs=set(),
        scaffold_lengths={"s1": 100_000}, length_from_end=5_000,
        virsorter_categories={"g0": "0", "g4": "0"},  # hallmark on both flanks
    )
    by_id = {row["query_id"]: (row["amg_flags"], row["auxiliary_score"])
             for row in out.iter_rows(named=True)}
    # g2 has B (centre of M-triple) AND would have scored 1 (hallmark both
    # flanks); v1 bumps it to 4.
    assert "B" in by_id["g2"][0]
    assert by_id["g2"][1] == 4
    # End genes still 5, no downgrade applies.
    assert by_id["g0"][1] == 5


def test_read_scaffold_lengths(tmp_path):
    """Scaffold lengths are read correctly across multi-line records, including
    trailing newlines and IDs with embedded spaces."""
    p = tmp_path / "tiny.fa"
    p.write_text(">contig_a some description\nACGT\nACGT\n>contig_b\nNNN\n")
    assert read_scaffold_lengths(str(p)) == {"contig_a": 8, "contig_b": 3}
