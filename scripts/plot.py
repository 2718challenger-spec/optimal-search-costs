"""Reproduce the saving-rate figure and normalized-rate figure from exact data.

Data processing and plotting are separated from the C++ interval DP.
Only plotting converts reduced fractions to floating point.
"""
import argparse
import math
from fractions import Fraction
from pathlib import Path
from common import DEFAULT_DATA, ROOT, read_rows

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
parser.add_argument("--output-dir", type=Path, default=ROOT / "figures")
args = parser.parse_args()
by = {}
for row in read_rows(args.data):
    c, g = int(row["optimal_sum"]), int(row["binary_sum"])
    row["saving_rate"] = str(Fraction(g-c, g))
    key = row["j"], row["k"]
    if key in by:
        raise ValueError(f"Duplicate data row: {key}")
    by[key] = row
n_max = max(j for j, k in by)
if set(by) != {(j,k) for j in range(1,n_max+1) for k in range(1,10)}:
    raise ValueError("Expected a complete grid of j=1,...,N and k=1,...,9")
output = args.output_dir
output.mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.labelcolor':'#253248','text.color':'#253248',
                     'axes.titleweight':'bold','grid.color':'#DCE2EA','grid.linewidth':0.6})
colors = ['#145B8D','#257A70','#A75B2B','#69569B','#327884','#9D5A75','#587E3E','#966A2D','#344A84']
fig, axes = plt.subplots(3,3,figsize=(14,10),layout='constrained')
for k, ax in enumerate(axes.flat,1):
    j = list(range(1, n_max + 1))
    y = [100*float(Fraction(by[n,k]['saving_rate'])) for n in j]
    ax.plot(j,y,lw=1.0,color=colors[k-1])
    dj = [2**h-1 for h in range(2, (n_max + 1).bit_length()) if 2**h-1 <= n_max]
    dy = [100*float(Fraction(by[n,k]['saving_rate'])) for n in dj]
    ax.scatter(dj,dy,s=16,color='#E2A324',edgecolor='white',linewidth=.35,zorder=3)
    ax.set_xscale('log',base=2)
    ax.set_xticks([1,4,16,64,256,1024],labels=['1','4','16','64','256','1024'])
    ax.set_title(f'k = {k}',loc='left')
    ax.set_ylim(bottom=0);ax.grid(True,alpha=.8)
    ax.set_xlabel('Interval size j');ax.set_ylabel('Saving rate (%)')
fig.suptitle('Exact optimal search vs. lower-midpoint binary search',fontsize=18,fontweight='bold')
fig.supxlabel(f'All integer j = 1,...,{n_max}. Gold markers: j = 2^h - 1. Exact fractions are rounded only for plotting.',fontsize=11)
fig.savefig(output/'saving_rates.png',dpi=180,bbox_inches='tight')
fig.savefig(output/'saving_rates.svg',bbox_inches='tight')
plt.close(fig)

if n_max < 32:
    raise SystemExit(0)

fig, axes = plt.subplots(3,3,figsize=(14,10),layout='constrained')
for k, ax in enumerate(axes.flat,1):
    j = list(range(32, n_max + 1))
    y = [float(Fraction(by[n,k]['saving_rate']))*math.log2(n) for n in j]
    ax.plot(j,y,lw=1.0,color=colors[k-1])
    ax.set_xscale('log',base=2)
    ax.set_xticks([32,64,128,256,512,1024,2048],labels=['32','64','128','256','512','1024','2048'])
    ax.tick_params(axis='x',labelsize=8)
    ax.set_title(f'k = {k}',loc='left')
    ax.set_ylim(bottom=0);ax.grid(True,alpha=.8)
    ax.set_xlabel('Interval size j');ax.set_ylabel(r'$R_{j,k}\,\log_2 j$')
fig.suptitle('Normalized saving rate: dyadic variation persists in the computed range',fontsize=17,fontweight='bold')
fig.supxlabel('Numerical evidence only for the shape; the proven order is R(j,k) = Theta(1 / log j).',fontsize=11)
fig.savefig(output/'normalized_rates.png',dpi=180,bbox_inches='tight')
fig.savefig(output/'normalized_rates.svg',bbox_inches='tight')
plt.close(fig)
