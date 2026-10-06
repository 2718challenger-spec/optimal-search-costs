"""Verify exact data independently of graph generation.

Checks all stored identities and the midpoint baseline; exhaustively enumerates
trees for j<=8. Optionally compares a fresh C++ output with the reference data.
"""
import argparse
import json
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path
from common import DEFAULT_DATA, read_rows
from query import solve


@lru_cache(None)
def all_shapes(n):
    """Inorder tuple of subtree sizes for each possible ordered binary tree."""
    if n == 0:
        return ((),)
    return tuple(left + (n,) + right
                 for count in range(n)
                 for left in all_shapes(count)
                 for right in all_shapes(n-count-1))


def check(condition, message):
    if not condition:
        raise ValueError(message)


def verify(data, recomputed=None):
    by = {}
    for row in read_rows(data):
        key = row['j'], row['k']
        check(key not in by, f'Duplicate row {key}')
        by[key] = row
    n_max = max(j for j, k in by)
    check(set(by) == {(j,k) for j in range(1,n_max+1) for k in range(1,10)},
          'Incomplete data grid')

    # Construct midpoint trees independently, then compute all their moments.
    for j in range(1, n_max+1):
        bm = [0]*10
        pending = [(1,j)]
        while pending:
            a,b = pending.pop()
            if a>b:
                continue
            q = (a+b)//2
            value = 1
            for r in range(10):
                bm[r] += (b-a+1)*value
                value *= j-q
            pending.extend(((a,q-1),(q+1,b)))
        for k in range(1,10):
            row = by[j,k]
            c,g = int(row['optimal_sum']), int(row['binary_sum'])
            m = list(map(int,row['moments']))
            check(len(m)==k+1, f'Wrong moment length at {(j,k)}')
            pc = [(-1)**r*comb(k,r)*m[r] for r in range(k+1)]
            bc = [(-1)**r*comb(k,r)*bm[r] for r in range(k+1)]
            check(sum(pc[r]*j**(k-r) for r in range(k+1))==c,
                  f'Optimal moment identity failed at {(j,k)}')
            check(sum(bc[r]*j**(k-r) for r in range(k+1))==g,
                  f'Midpoint baseline failed at {(j,k)}')
            check(0<c<=g and 1<=row['root']<=j, f'Invalid row {(j,k)}')
            expected = {
                'optimal_expectation': str(Fraction(c,j)),
                'binary_expectation': str(Fraction(g,j)),
                'saving': str(Fraction(g-c,j)),
                'saving_rate': str(Fraction(g-c,g)),
                'polynomial_coefficients': list(map(str,pc)),
                'binary_moments': list(map(str,bm[:k+1])),
                'binary_polynomial_coefficients': list(map(str,bc)),
                'saving_polynomial_coefficients': list(map(str,[b-a for a,b in zip(pc,bc)])),
            }
            for field,value in expected.items():
                if field in row:
                    check(row[field]==value, f'{field} failed at {(j,k)}')
            if k>=2 and j>=3:
                mid=(j+1)//2
                child=mid//2
                gain=child*mid**k-(j-mid+1)*child**k
                check(0<gain<=g-c, f'Rotation bound failed at {(j,k)}')
            if k==1 and (j+1)%32==0:
                check(((j+1)//32)**2<=g-c, f'Skeleton bound failed at {(j,k)}')

    for j in range(1,min(8,n_max)+1):
        for k in range(1,10):
            best = min((sum(size*q**k for q,size in enumerate(shape,1)),sum(shape))
                       for shape in all_shapes(j))
            row=by[j,k]
            check(best==(int(row['optimal_sum']),int(row['moments'][0])),
                  f'Exhaustive enumeration failed at {(j,k)}')

    compared = 0
    if recomputed:
        seen=set()
        for row in read_rows(recomputed):
            key=row['j'],row['k']
            check(key in by and key not in seen, f'Invalid recomputed row {key}')
            seen.add(key)
            for field in ['optimal_sum','binary_sum','root','moments']:
                check(row[field]==by[key][field], f'C++ comparison failed: {key} {field}')
            if row['j']<=8:
                independent=solve(*key)
                for field in ['optimal_sum','binary_sum','root','moments']:
                    check(row[field]==independent[field], f'Python DP failed: {key} {field}')
            compared+=1
        check(compared>0, 'Empty recomputed file')
        fresh_max=max(j for j,k in seen)
        check(seen=={(j,k) for j in range(1,fresh_max+1) for k in range(1,10)},
              'Incomplete recomputed grid')
    return {'verified_rows':len(by), 'maximum_j':n_max,
            'exhaustive_maximum_j':min(8,n_max),
            'recomputed_rows_compared':compared, 'status':'passed'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=DEFAULT_DATA)
    parser.add_argument('--recomputed',type=Path)
    args=parser.parse_args()
    print(json.dumps(verify(args.data,args.recomputed),indent=2))
