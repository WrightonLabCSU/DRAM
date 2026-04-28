"""pytest configuration: put bin/ on sys.path so test files can import modules
that the Nextflow processes ship as bin/<script>.py."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))
