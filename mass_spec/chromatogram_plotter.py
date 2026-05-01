import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
 
PALETTE = ["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000"]
 
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica Neue", "Arial"],
    "font.size": 6, "axes.linewidth": 0.8, "lines.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 3, "ytick.major.size": 3,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.facecolor": "white", "figure.facecolor": "white",
    "axes.edgecolor": "#555", "xtick.color": "#555",
    "ytick.color": "#555", "text.color": "#555",
})
 
 
def load(path):
    data = np.genfromtxt(path, delimiter=",", invalid_raise=False, encoding="utf-8-sig")
    data = data[~np.isnan(data).any(axis=1)]   # drop header/bad rows
    return data[:, 1], data[:, 2]              # Agilent: col1=time, col2=counts
 
 
def break_marks(ax_l, ax_r, d=0.015):
    for ax, x in [(ax_l, 1), (ax_r, 0)]:
        ax.plot([x, x], [-d, d], color="#555", lw=0.8,
                transform=ax.transAxes, clip_on=False)
 
 
def plot(datasets, break_range=None, width_ratios=None, stack=False,
         fill=True, output=None, figwidth=3.5):
 
    n = len(datasets)
    colors = [PALETTE[i % len(PALETTE)] for i in range(n)]
 
    if stack:
        max_y = max(np.nanmax(np.abs(y)) for _, y in datasets.values())
        offsets = [i * max_y * 1.2 for i in range(n)]
    else:
        offsets = [0] * n
 
    if break_range:
        x_lo, x_hi = break_range
        all_x = np.concatenate([x for x, _ in datasets.values()])
        if width_ratios is None:
            span = np.ptp(all_x)
            width_ratios = [max((x_lo - all_x.min()) / span, 0.15),
                            max((all_x.max() - x_hi) / span, 0.15)]
        fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(figwidth, 2.2),
            gridspec_kw={"width_ratios": width_ratios, "wspace": 0.04})
    else:
        fig, ax = plt.subplots(figsize=(figwidth, 2.2))
 
    handles = []
    for i, (label, (x, y)) in enumerate(datasets.items()):
        col, yp = colors[i], y + offsets[i]
 
        def draw(a, _yp=yp, _off=offsets[i], _col=col, _label=label):
            ln, = a.plot(x, _yp, color=_col, lw=0.9 if n > 4 else 1.1, label=_label)
            if fill:
                a.fill_between(x, _off, _yp, color=_col, alpha=0.10)
            return ln
 
        if break_range:
            handles.append(draw(ax_l)); draw(ax_r)
        else:
            handles.append(draw(ax))
 
    if break_range:
        ax_l.set_xlim(right=x_lo);  ax_r.set_xlim(left=x_hi)
        all_y = np.concatenate([y for _, y in datasets.values()])
        ylim = (all_y.min() - 0.05*np.ptp(all_y), all_y.max() + 0.10*np.ptp(all_y))
        ax_l.set_ylim(ylim);  ax_r.set_ylim(ylim)
        ax_l.set_ylabel("au");  ax_r.set_yticks([])
        ax_l.spines["right"].set_visible(False)
        ax_r.spines["left"].set_visible(False)
        fig.text(0.5, -0.04, "time (min)", ha="center", fontsize=6, color="#555")
        break_marks(ax_l, ax_r)
        if n > 1:
            ax_r.legend(handles=handles, frameon=False, handlelength=1.2)
    else:
        ax.set_xlabel("time (min)"); ax.set_ylabel("au")
        if n > 1:
            ax.legend(handles=handles, frameon=False, handlelength=1.2)
 
    fig.tight_layout()
    fig.savefig(output) if output else plt.show()
 
 
def main():
    p = argparse.ArgumentParser()
    p.add_argument("files", nargs="*")
    p.add_argument("--break-axis", nargs=2, type=float, metavar=("LO", "HI"))
    p.add_argument("--width-ratios", nargs=2, type=float)
    p.add_argument("--stack", action="store_true")
    p.add_argument("--no-fill", action="store_true")
    p.add_argument("--figure-width", type=float, default=3.5)
    p.add_argument("-o", "--output", default=None)
    args = p.parse_args()
 
    datasets = {Path(f).stem: load(f) for f in args.files}
    plot(datasets,
         break_range=tuple(args.break_axis) if args.break_axis else None,
         width_ratios=args.width_ratios,
         stack=args.stack,
         fill=not args.no_fill,
         output=args.output,
         figwidth=args.figure_width)
 
if __name__ == "__main__":
    main()