"""Exact expected costs for cost(q)=q**k; successful final query is charged.

Usage: python scripts/query.py 100 9
The supplied JSONL is used when the requested pair is present. Add --compute
to recompute using arbitrary-precision Python integer arithmetic.
"""
import argparse
import json
from fractions import Fraction
from math import comb
from pathlib import Path
from common import DEFAULT_DATA, read_rows


def solve(j, k):
    if j < 1 or not 1 <= k <= 9:
        raise ValueError('Require j >= 1 and 1 <= k <= 9')
    size = j+2
    d = [[0]*size for _ in range(size)]
    b = [[0]*size for _ in range(size)]
    a = [[0]*size for _ in range(size)]
    root = [[0]*size for _ in range(size)]
    power = [q**k for q in range(size)]
    for length in range(1,j+1):
        for left in range(1,j-length+2):
            right = left+length-1
            value, count, q = min(
                (length*power[q]+d[left][q-1]+d[q+1][right],
                 length+a[left][q-1]+a[q+1][right], q)
                for q in range(left,right+1))
            d[left][right], a[left][right], root[left][right] = value,count,q
            q = (left+right)//2
            b[left][right] = length*power[q]+b[left][q-1]+b[q+1][right]
    def moments(optimal):
        answer=[0]*(k+1)
        pending=[(1,j)]
        while pending:
            left,right=pending.pop()
            if left>right:
                continue
            q=root[left][right] if optimal else (left+right)//2
            for r in range(k+1):
                answer[r]+=(right-left+1)*(j-q)**r
            pending.extend(((left,q-1),(q+1,right)))
        return answer
    m, bm = moments(True),moments(False)
    pc=[(-1)**r*comb(k,r)*m[r] for r in range(k+1)]
    bc=[(-1)**r*comb(k,r)*bm[r] for r in range(k+1)]
    c,g=d[1][j],b[1][j]
    return dict(j=j,k=k,optimal_sum=str(c),binary_sum=str(g),root=root[1][j],
                optimal_expectation=str(Fraction(c,j)),binary_expectation=str(Fraction(g,j)),
                saving=str(Fraction(g-c,j)),saving_rate=str(Fraction(g-c,g)),
                moments=list(map(str,m)),binary_moments=list(map(str,bm)),
                polynomial_coefficients=list(map(str,pc)),
                binary_polynomial_coefficients=list(map(str,bc)),
                saving_polynomial_coefficients=[str(y-x) for x,y in zip(pc,bc)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("j", type=int)
    parser.add_argument("k", type=int)
    parser.add_argument("--compute", action="store_true",
                        help="Recompute with Python integer DP instead of looking up data.")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    if args.j < 1 or not 1 <= args.k <= 9:
        parser.error("Require j >= 1 and 1 <= k <= 9")
    if args.compute:
        record = solve(args.j, args.k)
    else:
        record = next((row for row in read_rows(args.data)
                       if (row["j"], row["k"]) == (args.j, args.k)), None)
        if record is None:
            parser.error("Pair not in dataset; use --compute to explicitly recompute.")
    print(json.dumps(record, ensure_ascii=False, indent=2))
