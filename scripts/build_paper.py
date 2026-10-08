#!/usr/bin/env python3
# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Compile the final manuscript without modifying paper/; fail on unresolved refs and overfull boxes."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import subprocess


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = root/"paper"/"trace_nonpreserving_bound.tex"
    out = root/"build"/"paper"
    out.mkdir(parents=True, exist_ok=True)
    if shutil.which("pdflatex") is None:
        parser.error("pdfLaTeX is required")
    shutil.copy2(source, out/source.name)
    for iteration in (1, 2):
        completed = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source.name],
            cwd=out, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (out/f"pass-{iteration}.txt").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode:
            print(completed.stdout)
            return completed.returncode
    text = (out/"trace_nonpreserving_bound.log").read_text(errors="replace")
    problems = [line for line in text.splitlines() if "Overfull" in line or
                ("Warning" in line and any(t in line.lower() for t in ["undefined", "rerun", "multiply defined"]))]
    if problems:
        print("\n".join(problems))
        return 1
    pdf = out/"trace_nonpreserving_bound.pdf"
    if shutil.which("pdftotext"):
        content = subprocess.check_output(["pdftotext", str(pdf), "-"], text=True)
        for expected in ["By Lemma 4", "Therefore Lemma 5", "Corollary 11", "By Proposition 9"]:
            if expected not in " ".join(content.split()):
                print("Missing expected cross-reference:", expected)
                return 1
    print(pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
