"""Smoke test: ogni script di `src/` (e `src/sim/`) si importa senza errori e risponde a `--help`.

Nessuna banda statistica: esito binario (codice di uscita 0 e testo "usage").
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from conftest import ROOT, SRC, load_script

PY_SCRIPTS = sorted(p for p in SRC.glob("*.py")) + sorted((SRC / "sim").glob("*.py"))


@pytest.mark.parametrize("path", PY_SCRIPTS, ids=lambda p: p.name)
def test_import(path):
    load_script(path.stem)


@pytest.mark.parametrize(
    "path", [p for p in PY_SCRIPTS if p.name != "scpipe.py"], ids=lambda p: p.name
)
def test_help(path):
    r = subprocess.run(
        [sys.executable, str(path), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert r.returncode == 0, r.stderr[-2000:]
    assert "usage" in r.stdout.lower()


def test_download_tcr_help():
    r = subprocess.run(
        ["bash", str(SRC / "download_tcr.sh"), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert r.returncode == 0, r.stderr
    assert "uso" in r.stdout.lower() or "usage" in r.stdout.lower()
