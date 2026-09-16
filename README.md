[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

===============================================================================
PROJECT: Chudnovsky Pi Engine
===============================================================================

OVERVIEW:
Calculates Pi (pi) to extreme arbitrary precision (N digits) using the world-record
Chudnovsky Ramanujan-type hypergeometric series formula with exact-integer binary 
splitting and C-accelerated gmpy2 big-integer arithmetic.

ALGORITHM & MATHEMATICS:
- Chudnovsky Formula:
    1 / pi = 12 * sum_{k=0}^{infinity} (-1)^k * (6k)! * (13591409 + 545140134 * k) / ((3k)! * (k!)^3 * 640320^(3k + 3/2))
- Convergence Rate: ~14.1816 decimal digits per series term.
- Binary Splitting: the series is summed in exact big-integer arithmetic
  and collapses to a single (P, Q, T) triple with one final division.
  (A 12-process chunked pool was removed: benchmarked 6x SLOWER than
  sequential at 100k digits due to giant-integer serialization overhead.)

NOTE: the previous revision contained a scaling bug that made the engine
output the single digit "3" for every requested precision; the committed
Pi_1000_digits.txt contained exactly that. Fixed here, output regenerated
and verified against the published decimal expansion of pi (OEIS A000796).

## Usage

```bash
python "Chudnovsky Pi Engine.py" --help
```

TESTS:
    pytest tests/
