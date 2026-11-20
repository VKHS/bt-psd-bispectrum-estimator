"""Helpers for multi-panel figure sheets: cards with a tag, title and subtitle, and the paper-reported marker style.

Model data (this repository, synthetic signals) use filled markers and the series palette; values reported in the
paper are always drawn as hollow black diamonds, so the two can never be confused.
"""
from matplotlib.colors import LinearSegmentedColormap

from btpsd.viz import INK, INK2, SEQ_BLUE

PAPER = dict(marker="D", mfc="#ffffff", mec=INK, mew=1.1, color=INK, ls="none", ms=4.4)
SEQ_ORANGE = ["#fde3d8", "#fbc3ac", "#f6a07d", "#ef7d4f", "#d85c2f", "#a8431f", "#6e2a12"]


def cmap_blue():
    return LinearSegmentedColormap.from_list("seqblue", ["#ffffff"] + SEQ_BLUE)


def cmap_div():
    """Diverging blue - neutral - orange map for signed data (neutral midpoint, no hue at zero)."""
    return LinearSegmentedColormap.from_list("div", [SEQ_BLUE[5], SEQ_BLUE[3], SEQ_BLUE[1], "#f4f4f2", SEQ_ORANGE[1], SEQ_ORANGE[3], SEQ_ORANGE[5]])


def cmap_orange():
    return LinearSegmentedColormap.from_list("seqorange", ["#ffffff"] + SEQ_ORANGE)


def card(d, x, y, w, h, tag, title, sub=None, plot=True, pad=(8.8, 7.4, 3.0, 9.2), fs=8.2):
    """White card with a tag/title and an optional subtitle; returns the plot-area axes when ``plot`` is true."""
    d.card(x, y, w, h)
    d.panel(x + 1.8, y + h - 3.3, tag, title, fs=fs)
    if sub:
        d.text(x + 1.8, y + h - 6.6, sub, fs=5.6, ha="left", color=INK2)
    if plot:
        left, bottom, right, top = pad
        return d.inset(x + left, y + bottom, w - left - right, h - bottom - top)
    return None


def labels(ax, xl=None, yl=None, fs=5.8):
    if xl:
        ax.set_xlabel(xl, fontsize=fs, labelpad=1.2)
    if yl:
        ax.set_ylabel(yl, fontsize=fs, labelpad=1.2)


def legend(ax, loc="upper right", fs=5.2, **kw):
    opts = dict(handlelength=1.5, handletextpad=0.5, borderpad=0.2, labelspacing=0.28)
    opts.update(kw)
    return ax.legend(fontsize=fs, frameon=False, loc=loc, **opts)
