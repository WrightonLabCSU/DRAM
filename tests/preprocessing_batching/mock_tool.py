#!/usr/bin/env python3
"""Small stand-ins that validate generated preprocessing commands, not biology."""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

args = sys.argv[1:]
tool = Path(sys.argv[0]).name


def option(name):
    return args[args.index(name) + 1]


def sample(path):
    with open(path) as handle:
        return handle.readline().strip()


def fail_if_requested(name):
    if os.environ.get("FAIL_SAMPLE") == name:
        sys.exit(42)


if tool == "python":
    source = Path(args[0]).read_text()
    assignments = re.findall(r"^batch_inputs = (\[.*\])$", source, re.MULTILINE)
    assert assignments, source
    batch_inputs = json.loads(assignments[-1])
    is_trna = "tRNAscan-SE" in source
    assert ('--thread", "1"' in source) if is_trna else ("threads=1" in source)
    for name, fasta in batch_inputs:
        assert Path(fasta).exists()
        fail_if_requested(name)
        suffix = "_processed_trnas.tsv" if is_trna else "_processed_rrnas.tsv"
        Path(name + suffix).write_text("fasta\tquery_id\n" + name + "\t" + name + "\n")
elif tool == "sbatch":
    script = Path(args[-1]).resolve()
    contents = script.read_text()
    array = re.search(r"^#SBATCH --array[= ](\d+)-(\d+)", contents, re.MULTILINE)
    indices = range(int(array[1]), int(array[2]) + 1) if array else [None]
    job_id = str(os.getpid())
    with Path(os.environ["MOCK_SLURM_LOG"]).open("a") as handle:
        handle.write(("array" if array else "single") + "\t" + str(script) + "\n")
    for index in indices:
        env = dict(os.environ, SLURM_JOB_ID=job_id)
        if index is not None:
            env.update(SLURM_ARRAY_JOB_ID=job_id, SLURM_ARRAY_TASK_ID=str(index))
        subprocess.run(
            ["bash", str(script)],
            env=env,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    print("Submitted batch job " + job_id)
elif tool in ["squeue", "scancel"]:
    pass
elif tool == "reformat.sh":
    values = dict(arg.split("=", 1) for arg in args)
    name = sample(values["in"])
    fail_if_requested(name)
    Path(values["out"]).write_text(name + "\n")
elif tool == "prodigal":
    name = sample(option("-i"))
    fail_if_requested(name)
    Path(option("-o")).write_text(f"{name}\tprodigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n")
    Path(option("-d")).write_text(name + "\n")
    Path(option("-a")).write_text(name + "\n")
elif tool == "parse_faa.sh":
    Path(args[1]).write_text(sample(args[0]) + "\n")
elif tool == "gff_replace_id_with_scaffold_gene_number.sh":
    Path(args[1]).write_text(Path(args[0]).read_text())
elif tool == "quast.py":
    output = Path(option("-o"))
    output.mkdir(exist_ok=True)
    (output / "report.tsv").write_text("quast\n")
    Path("quast.log").write_text("quast\n")
elif tool == "process_quast.py":
    Path("collected_quast.tsv").write_text("quast\n")
elif tool == "tRNAscan-SE":
    assert option("--thread") == "1"
    name = sample(args[-1])
    fail_if_requested(name)
    Path(option("-o")).write_text(
        "mock header\n"
        "Sequence\tName\tBegin\tEnd\tType\tCodon\tScore\tNote\n"
        "separator\n"
        f"{name}\t{name}\t1\t10\tAla\tTGC\t50\t\n"
    )
elif tool == "barrnap":
    assert option("--threads") == "1"
    name = sample(args[-1])
    fail_if_requested(name)
    print(f"{name}\tbarrnap\trRNA\t1\t10\t1e-5\t+\t.\tName=16S")
else:
    raise AssertionError((tool, args))
