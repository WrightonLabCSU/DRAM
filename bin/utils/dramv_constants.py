"""DRAM-v reference sets used to derive AMG flags.

Ported verbatim from DRAM v1 (commit 6cd68f9 in the upstream master line,
mag_annotator/annotate_vgfs.py). Kept here so the FASTA-only DRAM-v
Phase 1 path and any future VirSorter-aware Phase 2 path share the same
source of truth.
"""

# Pfam ids used to mark genes as transposons (T flag).
TRANSPOSON_PFAMS = {
    "PF01609", "PF00872", "PF01610", "PF01527", "PF02371", "PF01710",
    "PF01385", "PF01548", "PF01526", "PF01797", "PF02899", "PF05717",
    "PF07592", "PF03050", "PF04754", "PF04986", "PF03400",
}

# CAZy families that mediate viral attachment / host cell entry (A flag).
CELL_ENTRY_CAZYS = {
    "CBM50", "GH102", "GH103", "GH104", "GH108", "GH18", "GH19", "GH22",
    "GH23", "GH24", "GH25", "GH73", "PL9", "CBM12", "CBM14", "CBM18",
    "CBM19", "AA15", "AA10", "GH32", "GH58", "PL16",
}

# MEROPS peptidase families typical of viruses (P flag).
VIRAL_PEPTIDASES_MEROPS = {
    "A02H", "A02G", "A02F", "A02E", "A02D", "A02C", "A02B", "A02A",
    "A03B", "A03A", "A11B", "A11A", "A22B", "A22A", "A33",
    "C01B", "C01A", "C02B", "C02A", "C03H", "C03G", "C03F", "C03E",
    "C03D", "C03C", "C03B", "C03A", "C04", "C05", "C06", "C07", "C08",
    "C09", "C104", "C105", "C107", "C108", "C113", "C14B", "C14A",
    "C16B", "C16A", "C18", "C19", "C21", "C23", "C24", "C26", "C27",
    "C28", "C30", "C31", "C32", "C33", "C36", "C37", "C39", "C40",
    "C42", "C44", "C46", "C48", "C49", "C51", "C53", "C57", "C59",
    "C60B", "C60A", "C62", "C63", "C71", "C74", "C76", "C85B", "C85A",
    "C87", "C89", "C97", "C99",
    "G02",
    "I02", "I04", "I08", "I24", "I25C", "I25B", "I25A", "I29", "I32",
    "I36", "I43", "I50B", "I50A", "I51", "I63", "I75", "I87", "I91",
    "M03C", "M03B", "M03A", "M10C", "M10B", "M10A", "M12C", "M12B",
    "M12A", "M13", "M14C", "M14D", "M14B", "M14A", "M15D", "M15C",
    "M15B", "M15A", "M16C", "M16B", "M16A", "M20E", "M20A", "M20D",
    "M20F", "M20B", "M20C", "M23B", "M23A", "M27", "M38", "M41", "M42",
    "M43B", "M43A", "M44", "M48C", "M48B", "M48A", "M56", "M60",
    "M67C", "M67B", "M67A", "M78", "M79", "M86",
    "N01", "N02", "N04", "N05", "N07", "N08", "N09", "N10E", "N10D",
    "N10C", "N10B", "N11",
    "S01F", "S01E", "S01D", "S01C", "S01B", "S01A", "S03", "S06", "S07",
    "S08C", "S08B", "S08A", "S09D", "S09C", "S09B", "S09A", "S11", "S12",
    "S14", "S16", "S21", "S24", "S26C", "S26B", "S26A", "S28", "S29",
    "S30", "S31", "S32", "S33", "S49C", "S49B", "S49A", "S50", "S53",
    "S54", "S62", "S65", "S69", "S73", "S74", "S75", "S77", "S78",
    "S80", "S81",
    "T01B", "T01A", "T03",
    "U32",
    "U40",
}
