# Copyright (c) 2026 Petri Tanninen, Sentient Machine Corporation
# SPDX-License-Identifier: MIT
import json

from spectral_bounds.cli import main
from spectral_bounds.experiments import reproduce
from spectral_bounds.symbolic import exact_report


def test_symbolic():
    r = exact_report(5)
    assert r["all_checks_passed"]
    assert r["identities"]["replacement_rhs_minus_trace"] == "-a + d"


def test_reproduce(tmp_path):
    r = reproduce(tmp_path,seed=42,samples=2)
    assert r["all_internal_checks_passed"]
    assert r["counts"]["random_cp_maps"] == 10
    assert r["counts"]["outside_cp_cases"] == 8
    assert (tmp_path/"epsilon_limits.csv").exists()
    assert json.loads((tmp_path/"validation_summary.json").read_text())["seed"] == 42


def test_cli(capsys):
    assert main(["symbolic","--max-dimension","2"]) == 0
    assert json.loads(capsys.readouterr().out)["all_checks_passed"]
