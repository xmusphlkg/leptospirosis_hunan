#!/usr/bin/env python3
"""Compatibility entry point; all current figures are rendered by ggplot2 in R."""
from pathlib import Path
import argparse
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "manuscript/figures")
    parser.add_argument("--input-dir", type=Path,
                        default=ROOT / "output/hunan/current/network_redesign_2026-11-13")
    args = parser.parse_args()
    expected = ROOT / "output/hunan/current/network_redesign_2026-11-13"
    if args.input_dir.resolve() != expected:
        parser.error("The R figure suite uses the audited 2026-11-13 analysis directory.")
    runner = shutil.which("Rscript")
    if not runner:
        parser.error("Rscript and the ggplot2 dependencies are required.")
    subprocess.run([runner, str(ROOT / "script/hunan/plot_sentinel_network.R"),
                    str(args.output_dir.resolve())], check=True, cwd=ROOT)


if __name__ == "__main__":
    main()
