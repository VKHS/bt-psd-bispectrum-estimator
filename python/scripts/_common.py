"""Shared helpers for the scripts in this folder."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
RESULTS = ROOT / "results"
DATA = ROOT / "data"
FIGS = ROOT / "docs" / "figures"
VECTORS = ROOT / "sim" / "vectors"
