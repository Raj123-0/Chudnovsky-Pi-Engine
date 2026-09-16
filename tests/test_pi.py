"""Unit tests for the Chudnovsky Pi Engine.

The module lives in a file with spaces in its name, so it is loaded via
importlib rather than a normal import.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parent.parent / "Chudnovsky Pi Engine.py"
_spec = importlib.util.spec_from_file_location("pi_engine", MODULE_PATH)
pi_engine = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pi_engine)

# Published decimal expansion of pi (OEIS A000796).
KNOWN_DIGITS = "31415926535897932384626433832795028841971693993751"


class TestBinarySplitting:
    """The Chudnovsky binary-splitting machinery."""

    def test_single_terms_match_series(self):
        """Leaf values must match the Chudnovsky term structure."""
        from fractions import Fraction

        A, B = pi_engine.A, pi_engine.B
        C3_OVER_24 = pi_engine.C3_OVER_24

        # term ratio for k >= 1: P(k)/Q(k); k=0 contributes T = A
        for k in (1, 2, 5, 10):
            P, Q, T = pi_engine.bs_chudnovsky(k, k + 1)
            expected_P = (2 * k - 1) * (6 * k - 5) * (6 * k - 1)
            assert P == expected_P
            assert Q == k ** 3 * C3_OVER_24
            expected_T = expected_P * (A + B * k)
            if k % 2 == 1:
                expected_T = -expected_T
            assert T == expected_T

    def test_bs_equals_iterative_sum(self):
        """bs_chudnovsky(0, N) must equal the term-by-term series construction."""
        A, B = pi_engine.A, pi_engine.B
        C3_OVER_24 = int(pi_engine.C3_OVER_24)
        N = 30

        _, Q, T = pi_engine.bs_chudnovsky(0, N)

        # Iteratively rebuild the same exact (Q, T) the binary splitting produces.
        # State after k terms: Q = Q(0,k), T = T(0,k), P = P(0,k).
        P, Qacc, Tacc = 1, 1, int(A)  # k = 0 state
        for k in range(1, N):
            p = (2 * k - 1) * (6 * k - 5) * (6 * k - 1)
            q = k ** 3 * C3_OVER_24
            t = p * int(A + B * k) * (-1 if k % 2 else 1)
            Tacc = Tacc * q + P * t
            P = P * p
            Qacc = Qacc * q

        assert int(T) == Tacc
        assert int(Q) == Qacc


class TestComputePi:
    """End-to-end digit generation."""

    def test_known_digits(self):
        """First 50 digits match the published expansion of pi."""
        assert pi_engine.compute_pi(len(KNOWN_DIGITS)) == KNOWN_DIGITS

    def test_prefix_property(self):
        """Digits computed at lower precision are a prefix of higher precision."""
        short = pi_engine.compute_pi(15)
        long = pi_engine.compute_pi(150)
        assert long.startswith(short)

    def test_requested_length_returned(self):
        for n in (1, 2, 10, 77):
            assert len(pi_engine.compute_pi(n)) == n

    def test_rejects_nonpositive_digits(self):
        with pytest.raises(ValueError):
            pi_engine.compute_pi(0)
        with pytest.raises(ValueError):
            pi_engine.compute_pi(-5)


class TestOeisFiles:
    """OEIS output file generation."""

    def test_files_written_correctly(self, tmp_path, monkeypatch):
        """Both output files appear with the expected names, content, and b-file format."""
        monkeypatch.chdir(tmp_path)
        digits = pi_engine.compute_pi(20)
        pi_engine.save_oeis_files("Pi", digits, 20)

        raw = tmp_path / "Pi_20_digits.txt"
        bfile = tmp_path / "b_file_Pi_20.txt"

        assert raw.read_text() == digits
        lines = bfile.read_text().splitlines()
        assert len(lines) == 20
        assert lines[0] == f"1 {digits[0]}"
        assert lines[19] == f"20 {digits[19]}"

    def test_extra_digits_truncated(self, tmp_path, monkeypatch):
        """A digit string longer than target_digits is truncated."""
        monkeypatch.chdir(tmp_path)
        pi_engine.save_oeis_files("Test", "31415926789", 5)
        assert (tmp_path / "Test_5_digits.txt").read_text() == "31415"


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
