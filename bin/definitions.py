from pathlib import Path
import os

ROOT = Path(__file__).resolve().parent.parent
DRAM_DIR = ROOT / "bin"

IS_WINDOWS = os.name == "nt"
