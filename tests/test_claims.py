"""Runs the claim ledger in ``scripts/verify_claims.py`` under pytest.

Each entry there checks one numbered identity from the manuscript against an
independent evaluation of its definition.  This test fails if any residual
exceeds its budget, so a regression in the code that produced a figure cannot
silently invalidate a stated result.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import verify_claims as vc  # noqa: E402


def test_every_claim_in_the_ledger_passes():
    vc.RESULTS.clear()
    vc.check_szego()
    vc.check_unsharp()
    vc.check_gksl()
    vc.check_winding()
    vc.check_kg()
    vc.check_universal_bound()
    assert vc.RESULTS, "the ledger produced no checks"
    failures = [
        (tag, residual, tol, statement)
        for tag, statement, residual, _note, tol in vc.RESULTS
        if tol == tol and residual >= tol          # tol != nan
    ]
    assert not failures, failures


def test_the_ledger_covers_the_headline_identities():
    vc.RESULTS.clear()
    vc.check_szego()
    vc.check_unsharp()
    vc.check_gksl()
    vc.check_winding()
    vc.check_kg()
    vc.check_universal_bound()
    tags = " ".join(t for t, *_ in vc.RESULTS)
    for needed in ("szego", "static", "lindblad", "winding", "kg", "proper", "universal"):
        assert needed in tags, needed
