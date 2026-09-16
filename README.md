[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

===============================================================================
PROJECT: Chudnovsky Pi Engine
===============================================================================

OVERVIEW:
Calculates Pi (pi) to extreme arbitrary precision (N digits) using the world-record
Chudnovsky Ramanujan-type hypergeometric series formula with parallel binary splitting 
and C-accelerated gmpy2 big-integer arithmetic.

ALGORITHM & MATHEMATICS:
- Chudnovsky Formula:
    1 / pi = 12 * sum_{k=0}^{infinity} (-1)^k * (6k)! * (13591409 + 545140134 * k) / ((3k)! * (k!)^3 * 640320^(3k + 3/2))
- Convergence Rate: ~14.1816 decimal digits per series term.
- Multi-Core Parallel Binary Splitting: Divides range [0, K) into a balanced binary tree 
  evaluating P, Q, and T integer products across all available CPU cores.

## Usage

```bash
python "Chudnovsky Pi Engine.py" --help
```
