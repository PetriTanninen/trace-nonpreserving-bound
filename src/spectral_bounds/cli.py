# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
"""Command-line interface; failures return nonzero status."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .experiments import reproduce
from .symbolic import exact_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Reproducibility diagnostics, not a proof certificate.")
    sub = parser.add_subparsers(dest="command", required=True)
    rep = sub.add_parser("reproduce", help="Run exact and seeded numerical checks")
    rep.add_argument("--output", type=Path, default=Path("results/reproduced"))
    rep.add_argument("--seed", type=int, default=20261007)
    rep.add_argument("--samples", type=int, default=40)
    sym = sub.add_parser("symbolic", help="Run exact algebra checks")
    sym.add_argument("--max-dimension", type=int, default=6)
    args = parser.parse_args(argv)
    if args.command == "reproduce":
        if args.samples < 1 or args.seed < 0:
            parser.error("samples must be positive and seed nonnegative")
        result = reproduce(args.output, seed=args.seed, samples=args.samples)
    else:
        if args.max_dimension < 2:
            parser.error("max-dimension must be at least 2")
        result = exact_report(args.max_dimension)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0
