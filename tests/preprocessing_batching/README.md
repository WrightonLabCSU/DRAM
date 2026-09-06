# Preprocessing batching integration checks

These checks execute the real `CALL_GENES`, `TRNA_SCAN`, and `RRNA_SCAN` process
scripts with small mock executables. They cover default one-input tasks, batching,
per-input byte limits, output unpacking, resource caps, resume behavior, member
failure, and batching within Slurm job arrays.

```bash
NXF_VER=25.04.7 python3 -m unittest discover -s tests/preprocessing_batching -v
NXF_VER=25.10.0 python3 -m unittest discover -s tests/preprocessing_batching -v
```
