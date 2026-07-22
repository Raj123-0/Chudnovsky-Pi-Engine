#!/usr/bin/env python3
"""
Chudnovsky Pi Engine (HPC OEIS Edition)
=======================================
Calculates Pi (pi) to exactly [N] significant digits using the Chudnovsky 
hypergeometric series, parallel binary-splitting across exactly 12 worker processes, 
C-accelerated gmpy2 big-integer arithmetic, and strict truncation for OEIS database submission.

Mathematical Formula:
  1 / pi = 12 * sum_{k=0}^{infinity} (-1)^k * (6k)! * (13591409 + 545140134 * k) / ((3k)! * (k!)^3 * 640320^(3k + 3/2))

Chunking Strategy & Parallelization:
  The summation domain [0, K) is partitioned into 12 contiguous mathematical intervals.
  Each chunk is evaluated independently via a 12-worker multiprocessing.Pool, returning 
  compact (P, Q, T) integer triples to minimize IPC serialization overhead.
"""

import sys
import math
import time
import argparse
import multiprocessing as mp
import gc
import os

os.environ['MPMATH_GMPY2'] = '1'
import gmpy2
from gmpy2 import mpz, isqrt

sys.set_int_max_str_digits(0)

# Constants
NUM_WORKERS = 12
C = 640320
C3_OVER_24 = mpz(C)**3 // 24 # 10939058860032000
A = mpz(13591409)
B = mpz(545140134)

def bs_chudnovsky_range(a, b):
    """Binary splitting over interval [a, b)."""
    if b - a == 1:
        if a == 0:
            P = mpz(1)
            Q = mpz(1)
            T = A
        else:
            P = mpz(2*a - 1) * mpz(6*a - 5) * mpz(6*a - 1)
            Q = mpz(a)**3 * C3_OVER_24
            T = P * (A + B * a)
            if a % 2 == 1:
                T = -T
        return P, Q, T

    m = (a + b) // 2
    P1, Q1, T1 = bs_chudnovsky_range(a, m)
    P2, Q2, T2 = bs_chudnovsky_range(m, b)

    P = P1 * P2
    Q = Q1 * Q2
    T = T1 * Q2 + P1 * T2
    return P, Q, T

def worker_chunk(args):
    """Worker task evaluating a specific mathematical chunk."""
    a, b = args
    return bs_chudnovsky_range(a, b)

def save_oeis_files(constant_name, digits_str, target_digits):
    """Saves raw digit string and OEIS b-file format."""
    clean_digits = digits_str.replace(".", "")[:target_digits]
    
    # 1. Raw digits text file
    raw_filename = f"{constant_name}_{target_digits}_digits.txt"
    with open(raw_filename, "w", encoding="utf-8") as f:
        f.write(clean_digits)
    print(f"Saved raw digit output to {raw_filename}")

    # 2. OEIS b-file format
    b_filename = f"b_file_{constant_name}_{target_digits}.txt"
    with open(b_filename, "w", encoding="utf-8") as f:
        for idx, digit in enumerate(clean_digits, start=1):
            f.write(f"{idx} {digit}\n")
    print(f"Saved OEIS b-file output to {b_filename}")

def compute_pi_hpc(target_digits):
    # Safety margin of +50 digits to prevent rounding drift
    dps_working = target_digits + 50
    terms = int(math.ceil(dps_working / 14.18164742394)) + 2

    # Break terms into exactly 12 contiguous chunks
    chunk_size = math.ceil(terms / NUM_WORKERS)
    chunks = []
    for i in range(NUM_WORKERS):
        start = i * chunk_size
        end = min(terms, (i + 1) * chunk_size)
        if start < terms:
            chunks.append((start, end))

    # Parallel chunk execution on 12 cores
    with mp.Pool(processes=NUM_WORKERS) as pool:
        results = pool.map(worker_chunk, chunks)

    # Main process aggregation
    P, Q, T = results[0]
    for P_next, Q_next, T_next in results[1:]:
        P = P * P_next
        Q = Q * Q_next
        T = T * Q_next + P * T_next

    del results
    gc.collect()

    # Square root and final high-precision division
    prec_bits = int(dps_working * 3.3219280948873626) + 200
    gmpy2.get_context().precision = prec_bits

    scale = mpz(10)**(dps_working + 10)
    sqrt_10005 = isqrt(mpz(10005) * scale * scale)
    
    numerator = Q * mpz(426880) * sqrt_10005
    pi_int = numerator // (T * scale)

    pi_raw = str(pi_int)
    # Strictly truncate at N digits
    pi_truncated = pi_raw[:target_digits]
    
    del pi_int, numerator, sqrt_10005
    gc.collect()

    save_oeis_files("Pi", pi_truncated, target_digits)
    return pi_truncated

def main():
    parser = argparse.ArgumentParser(description="HPC Pi OEIS Calculator")
    parser.add_argument("-n", "--digits", type=int, default=1000, help="Target digits (default: 1000)")
    args = parser.parse_args()

    t0 = time.time()
    digits = compute_pi_hpc(args.digits)
    t1 = time.time()

    print(f"Execution finished in {t1 - t0:.4f} seconds using {NUM_WORKERS} cores.")

if __name__ == "__main__":
    main()
