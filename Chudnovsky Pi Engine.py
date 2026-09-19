#!/usr/bin/env python3
"""
Chudnovsky Pi Engine (OEIS Edition)
====================================
Computes pi to exactly N significant digits using the Chudnovsky
hypergeometric series with binary splitting in exact big-integer
arithmetic (gmpy2), and strict truncation for OEIS database submission.

Mathematical Formula
-------------------
    1/pi = 12 * sum_{k>=0} (-1)^k (6k)! (13591409 + 545140134 k)
           / ((3k)! (k!)^3 640320^(3k + 3/2))

The series is evaluated by binary splitting: each interval [a, b)
collapses to an exact integer triple (P, Q, T), and the final value is
produced by one integer square root and one division:

    pi = Q * 426880 * isqrt(10005 * 10^(2s)) / T   (scaled by 10^s)

Note (bug fix): the previous revision divided by (T * scale) instead of T,
which produced the single digit "3" for every requested precision — the
committed Pi_1000_digits.txt contained exactly that. The formula is
corrected here and the output files regenerated.

Note (parallelism): the previous 12-process chunked pool was removed —
benchmarked 6x SLOWER than sequential binary splitting even at 100k digits,
because pickling the giant intermediate integers across workers dominates
any parallel gain at these problem sizes.
"""

from __future__ import annotations

import argparse
import functools
import math
import sys
import time

from gmpy2 import mpz, isqrt



# str(pi_scaled) can exceed CPython's default 4300-digit int->str limit.
sys.set_int_max_str_digits(0)

C = 640320
C3_OVER_24 = mpz(C) ** 3 // 24  # 10939058860032000
A = mpz(13591409)
B = mpz(545140134)
DIGITS_PER_TERM = 14.18164742394  # log10(640320^3 / 24)


@functools.lru_cache(maxsize=None)
def bs_chudnovsky(a: int, b: int) -> tuple[mpz, mpz, mpz]:
    """Binary splitting of the Chudnovsky series over the interval [a, b).

    Returns the exact integer triple (P, Q, T) with
        P(a,b) = prod of term-ratio numerators,
        Q(a,b) = prod of term-ratio denominators,
        T(a,b) = the partial sum scaled by Q(a,b).
    """
    if b - a == 1:
        if a == 0:
            return mpz(1), mpz(1), A
        P = mpz(2 * a - 1) * mpz(6 * a - 5) * mpz(6 * a - 1)
        Q = mpz(a) ** 3 * C3_OVER_24
        T = P * (A + B * a)
        if a % 2 == 1:
            T = -T
        return P, Q, T
    m = (a + b) // 2
    P1, Q1, T1 = bs_chudnovsky(a, m)
    P2, Q2, T2 = bs_chudnovsky(m, b)
    return P1 * P2, Q1 * Q2, T1 * Q2 + P1 * T2


def compute_pi(target_digits: int) -> str:
    """Compute the first `target_digits` digits of pi.

    Args:
        target_digits: number of decimal digits to produce (>= 1).

    Returns:
        Digit string of length `target_digits` (leading '3' included,
        decimal point removed — OEIS b-file convention).

    Raises:
        ValueError: if target_digits is not a positive integer.
    """
    if target_digits < 1:
        raise ValueError("target_digits must be a positive integer")

    # +50 guard digits prevent rounding drift at the truncation boundary.
    dps_working = target_digits + 50
    terms = int(math.ceil(dps_working / DIGITS_PER_TERM)) + 2

    _, Q, T = bs_chudnovsky(0, terms)

    # pi scaled by 10^(dps_working + 10): one integer sqrt, one division.
    scale = mpz(10) ** (dps_working + 10)
    sqrt_10005 = isqrt(mpz(10005) * scale * scale)
    pi_scaled = (Q * mpz(426880) * sqrt_10005) // T

    return str(pi_scaled)[:target_digits]


def save_oeis_files(constant_name: str, digits_str: str, target_digits: int) -> None:
    """Write the raw digit string and an OEIS b-file to the current directory.

    Args:
        constant_name: name used in output filenames.
        digits_str: decimal digit string (leading integer-part digit included).
        target_digits: expected number of digits (extra digits are truncated).
    """
    clean_digits = digits_str.replace(".", "")[:target_digits]

    raw_filename = f"{constant_name}_{target_digits}_digits.txt"
    with open(raw_filename, "w", encoding="utf-8") as f:
        f.write(clean_digits)
    print(f"Saved raw digit output to {raw_filename}")

    b_filename = f"b_file_{constant_name}_{target_digits}.txt"
    with open(b_filename, "w", encoding="utf-8") as f:
        # writelines consumes the generator lazily — O(1) memory even for
        # millions of digits.
        f.writelines(f"{idx} {digit}\n" for idx, digit in enumerate(clean_digits, start=1))
    print(f"Saved OEIS b-file output to {b_filename}")


def main() -> None:
    """Entry point: parse arguments, compute, and save OEIS output files."""
    parser = argparse.ArgumentParser(description="Pi OEIS Calculator")
    parser.add_argument(
        "-n", "--digits", type=int, default=1000,
        help="Target digits (default: 1000)")
    args = parser.parse_args()

    if args.digits < 1:
        parser.error("--digits must be a positive integer")

    t0 = time.time()
    digits = compute_pi(args.digits)
    t1 = time.time()

    save_oeis_files("Pi", digits, args.digits)
    print(f"Execution finished in {t1 - t0:.4f} seconds.")


if __name__ == "__main__":
    main()
