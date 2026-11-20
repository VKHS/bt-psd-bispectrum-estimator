"""Plot/diagram helpers: palette and block-diagram primitives (matplotlib only).

Colours follow a validated categorical palette (light surface); text always uses ink
tokens, never the series colour.  Diagrams are drawn in 'drawing units' of 0.086 inch.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle  # noqa: E402

PAL = {
    "blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
    "magenta": "#e87ba4", "green": "#008300", "violet": "#6250d6", "red": "#e34948",
}
SERIES = [PAL[k] for k in ("blue", "orange", "aqua", "yellow", "magenta", "green", "violet", "red")]
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8984", "#fcfcfb", "#e4e3df"
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

#: block role -> palette colour (identity only; the label carries the meaning)
KIND = {
    "memory": "blue", "compute": "aqua", "control": "violet", "rom": "orange",
    "accum": "yellow", "io": None, "free": "magenta",
}


def tint(hex_color: str, a: float = 0.84) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    mix = lambda c: int(round(c + (255 - c) * a))  # noqa: E731
    return "#%02x%02x%02x" % (mix(r), mix(g), mix(b))


def style_axes(ax, grid=True):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8, width=0.8)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.title.set_color(INK)
    if grid:
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)


def savefig(fig, path_stem: Path, dpi: int = 200):
    path_stem = Path(path_stem)
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_stem.with_suffix(".svg"), facecolor=SURF, metadata={"Date": None})
    fig.savefig(path_stem.with_suffix(".png"), dpi=dpi, facecolor=SURF)
    plt.close(fig)


plt.rcParams["svg.hashsalt"] = "bt-psd-bispectrum-estimator"  # reproducible ids: regenerating a figure changes the SVG only if the figure changed
UNIT = 0.086  # inches per drawing unit


@dataclass
class Box:
    x: float
    y: float
    w: float
    h: float

    @property
    def l(self): return (self.x, self.y + self.h / 2)  # noqa: E743
    @property
    def r(self): return (self.x + self.w, self.y + self.h / 2)
    @property
    def t(self): return (self.x + self.w / 2, self.y + self.h)
    @property
    def b(self): return (self.x + self.w / 2, self.y)
    @property
    def c(self): return (self.x + self.w / 2, self.y + self.h / 2)

    def at(self, fx, fy):
        return (self.x + fx * self.w, self.y + fy * self.h)


class View:
    """Translated drawing context: coordinates given to its methods are relative to ``(dx, dy)``.

    Lets a block diagram that was laid out from the origin be moved up a sheet that also carries data panels.
    Returned ``Box`` objects are in the local (shifted) frame, so their ports can be fed back to the same view.
    """

    def __init__(self, d: "Diagram", dx: float = 0.0, dy: float = 0.0):
        self.d, self.dx, self.dy = d, dx, dy

    @property
    def ax(self):
        return self.d.ax

    def _l(self, fn, x, y, *a, **k):
        return fn(x + self.dx, y + self.dy, *a, **k)

    def box(self, x, y, *a, **k):
        b = self._l(self.d.box, x, y, *a, **k)
        return Box(b.x - self.dx, b.y - self.dy, b.w, b.h)

    def text(self, x, y, *a, **k): return self._l(self.d.text, x, y, *a, **k)
    def panel(self, x, y, *a, **k): return self._l(self.d.panel, x, y, *a, **k)
    def junction(self, x, y, *a, **k): return self._l(self.d.junction, x, y, *a, **k)
    def adder(self, x, y, *a, **k): return self._l(self.d.adder, x, y, *a, **k)
    def rect(self, x, y, *a, **k): return self._l(self.d.rect, x, y, *a, **k)
    def card(self, x, y, *a, **k): return self._l(self.d.card, x, y, *a, **k)
    def inset(self, x, y, *a, **k): return self._l(self.d.inset, x, y, *a, **k)
    def legend(self, x, y, *a, **k): return self._l(self.d.legend, x, y, *a, **k)

    def arrow(self, pts, *a, **k):
        return self.d.arrow([(x + self.dx, y + self.dy) for x, y in pts], *a, **k)

    def view(self, dx, dy):
        return View(self.d, self.dx + dx, self.dy + dy)


class Diagram:
    """Block-diagram canvas: 1 drawing unit = ``unit`` inches (default 0.086).

    Boxes have a coloured header band, a white body and a soft shadow; ``inset`` places ordinary
    matplotlib axes at diagram coordinates so data panels and block diagrams share one figure sheet.
    """

    def __init__(self, w: float, h: float, title: str | None = None, unit: float = UNIT, subtitle: str | None = None):
        self.w, self.h, self.u = w, h, unit
        self.fig, self.ax = plt.subplots(figsize=(w * unit, h * unit))
        self.ax.set_xlim(0, w)
        self.ax.set_ylim(0, h)
        self.ax.set_aspect("equal")
        self.ax.axis("off")
        self.fig.subplots_adjust(0, 0, 1, 1)
        self.fig.patch.set_facecolor(SURF)
        if title:
            self.ax.text(1.5, h - 2.0, title, fontsize=11.5, fontweight="bold", color=INK, va="top", ha="left")
        if subtitle:
            self.ax.text(1.5, h - 5.6, subtitle, fontsize=7.6, color=INK2, va="top", ha="left")

    # ---- primitives -----------------------------------------------------
    def _lh(self, fs):
        return fs * 1.14 / (72.0 * self.u)

    def box(self, x, y, w, h, text="", kind="compute", fs=7.2, step=None, bold=True, round_=1.0, lw=1.0, dashed=False,
            header=True):
        col = KIND.get(kind)
        edge = PAL[col] if col else INK2
        band = tint(PAL[col], 0.78) if col else "#e7e6e1"
        light = tint(PAL[col], 0.93) if col else "#f8f7f5"
        ls = "--" if dashed else "-"
        style = f"round,pad=0,rounding_size={round_}"
        lines = text.split("\n") if text else []
        banded = header and bold and len(lines) > 1
        self.ax.add_patch(FancyBboxPatch((x + 0.30, y - 0.30), w, h, boxstyle=style, fc="#000000", ec="none",
                                         alpha=0.08, zorder=1))
        body = FancyBboxPatch((x, y), w, h, boxstyle=style, fc=light if banded else band, ec="none", zorder=2)
        self.ax.add_patch(body)
        if banded:
            hh = self._lh(fs) + 2.0
            hh = min(hh, h * 0.5)
            r = Rectangle((x, y + h - hh), w, hh, fc=band, ec="none", zorder=2.1)
            self.ax.add_patch(r)
            r.set_clip_path(body)
            self.ax.plot([x, x + w], [y + h - hh] * 2, color=edge, lw=0.6, alpha=0.55, zorder=2.15, solid_capstyle="butt")
        self.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, fc="none", ec=edge, lw=lw, ls=ls, zorder=2.2))
        if lines:
            if banded:
                lh = self._lh(fs - 0.4)
                n = len(lines) - 1
                self.ax.text(x + w / 2, y + h - hh / 2, lines[0], fontsize=fs, fontweight="bold", color=INK, ha="center",
                             va="center", zorder=3)
                yc = y + (h - hh) / 2
                for i, ln in enumerate(lines[1:]):
                    self.ax.text(x + w / 2, yc + ((n - 1) / 2 - i) * lh, ln, fontsize=fs - 0.4, color=INK2, ha="center",
                                 va="center", zorder=3)
            else:
                self.ax.text(x + w / 2, y + h / 2, text, fontsize=fs, color=INK, ha="center", va="center",
                             fontweight="bold" if bold and len(lines) == 1 else "normal", linespacing=1.15, zorder=3)
        if step is not None:
            self.ax.add_patch(Circle((x + 0.2, y + h - 0.2), 1.35, fc=INK, ec=SURF, lw=1.0, zorder=5))
            self.ax.text(x + 0.2, y + h - 0.2, str(step), color=SURF, fontsize=6.5, ha="center", va="center",
                         fontweight="bold", zorder=6)
        return Box(x, y, w, h)

    def text(self, x, y, s, fs=7, color=INK2, ha="center", va="center", bold=False, italic=False, rot=0, mono=False):
        self.ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va, fontweight="bold" if bold else "normal",
                     fontstyle="italic" if italic else "normal", rotation=rot, linespacing=1.15, zorder=6,
                     family="DejaVu Sans Mono" if mono else None)

    def panel(self, x, y, tag, title="", fs=8.6):
        """Panel label such as '(a)' followed by a bold title."""
        self.ax.text(x, y, tag, fontsize=fs, fontweight="bold", color=INK, ha="left", va="center", zorder=6)
        if title:
            self.ax.text(x + 0.82 * len(tag) * fs / 8.6, y, title, fontsize=fs - 0.6, fontweight="bold", color=INK,
                         ha="left", va="center", zorder=6)

    def arrow(self, pts, label=None, color=INK2, lw=1.05, ls="-", lpos=0.5, loff=(0, 1.3), fs=6.6, head=True,
              ha="center", bus=None):
        """Polyline arrow through ``pts``; ``bus`` draws a width slash and label on the first segment."""
        xs, ys = zip(*pts)
        if head and len(pts) >= 2:
            if len(pts) > 2:
                self.ax.plot(xs[:-1], ys[:-1], color=color, lw=lw, ls=ls, solid_capstyle="butt", zorder=2.5)
            self.ax.plot([xs[-2], xs[-1]], [ys[-2], ys[-1]], color=color, lw=lw, ls=ls, solid_capstyle="butt", zorder=2.5)
            self.ax.annotate("", xy=pts[-1], xytext=pts[-2],
                             arrowprops=dict(arrowstyle="-|>,head_length=0.38,head_width=0.2", color=color, lw=lw,
                                             shrinkA=0, shrinkB=0, ls="-"), zorder=2.6)
        else:
            self.ax.plot(xs, ys, color=color, lw=lw, ls=ls, zorder=2.5)
        if bus:
            (x0, y0), (x1, y1) = pts[0], pts[1]
            mx, my = (x0 + x1) / 2, (y0 + y1) / 2
            dx, dy = x1 - x0, y1 - y0
            n = (dx * dx + dy * dy) ** 0.5 or 1.0
            ux, uy = dx / n, dy / n
            self.ax.plot([mx - 0.45 * uy - 0.25 * ux, mx + 0.45 * uy + 0.25 * ux],
                         [my + 0.45 * ux - 0.25 * uy, my - 0.45 * ux + 0.25 * uy], color=color, lw=0.9, zorder=2.7)
            off = (0.0, 1.35) if abs(dx) >= abs(dy) else (1.1, 0.0)
            self.ax.text(mx + off[0], my + off[1], bus, fontsize=5.6, color=INK2, ha="center" if abs(dx) >= abs(dy) else "left",
                         va="center", family="DejaVu Sans Mono", zorder=6)
        if label:
            i = min(int(lpos * (len(pts) - 1)), len(pts) - 2)
            fx = lpos * (len(pts) - 1) - i
            x = pts[i][0] + fx * (pts[i + 1][0] - pts[i][0]) + loff[0]
            y = pts[i][1] + fx * (pts[i + 1][1] - pts[i][1]) + loff[1]
            self.text(x, y, label, fs=fs, ha=ha)

    def junction(self, x, y, r=0.35):
        self.ax.add_patch(Circle((x, y), r, fc=INK2, ec=INK2, zorder=4))

    def adder(self, x, y, sign="+", r=1.15, color=INK):
        self.ax.add_patch(Circle((x, y), r, fc="#ffffff", ec=color, lw=1.0, zorder=4))
        self.ax.text(x, y - 0.05, sign, fontsize=7.5, color=color, ha="center", va="center", fontweight="bold", zorder=5)

    def rect(self, x, y, w, h, fc=None, ec=MUTED, lw=0.8, ls="--", label=None, fs=7):
        self.ax.add_patch(Rectangle((x, y), w, h, fc=fc or "none", ec=ec, lw=lw, ls=ls, zorder=0.5))
        if label:
            self.text(x + 1.0, y + h - 1.0, label, fs=fs, ha="left", va="top", color=INK2, italic=True)

    def card(self, x, y, w, h, ec=GRID, fc="#ffffff", lw=0.9, round_=1.2):
        """Plain white card used as a panel background."""
        self.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={round_}", fc=fc, ec=ec,
                                         lw=lw, zorder=0.6))

    def inset(self, x, y, w, h, title=None, fs=6.2, grid=True):
        """Matplotlib axes placed at diagram coordinates (plot area only; tick labels spill outside)."""
        ax = self.fig.add_axes([x / self.w, y / self.h, w / self.w, h / self.h], zorder=3)
        style_axes(ax, grid=grid)
        ax.tick_params(labelsize=fs, length=2, pad=1.5, width=0.6)
        ax.patch.set_facecolor("#ffffff")
        for s in ("left", "bottom"):
            ax.spines[s].set_linewidth(0.7)
        if title:
            ax.set_title(title, fontsize=fs + 0.8, color=INK, loc="left", pad=3.5, fontweight="bold")
        return ax

    def legend(self, x, y, items, fs=6.8, dx=24):
        """Horizontal legend strip of (kind, label)."""
        for i, (kind, lab) in enumerate(items):
            col = KIND.get(kind)
            self.ax.add_patch(FancyBboxPatch((x + i * dx, y - 0.9), 2.6, 1.8, boxstyle="round,pad=0,rounding_size=0.4",
                                             fc=tint(PAL[col], 0.78) if col else "#e7e6e1", ec=PAL[col] if col else INK2, lw=0.9))
            self.text(x + i * dx + 3.4, y, lab, fs=fs, ha="left")

    def view(self, dx: float = 0.0, dy: float = 0.0) -> "View":
        return View(self, dx, dy)

    def save(self, stem):
        savefig(self.fig, stem)
