"""
Centralised plot style for the mt_dsnow project.

Usage
-----
    from plot_style import apply_style, C, ALPHA, LS, HATCH, SUBSET_COLOR, FIG

    apply_style()          # call once per script / at the top of a notebook

    ax.plot(t, swe, color=C.DSNOW,     label="ΔSnow")
    ax.plot(t, swe, color=C.HS2SWE,    label="HS2SWE")
    ax.plot(t, swe, color=C.OBS,       label="CRNS median")
    ax.fill_between(t, lo, hi, color=C.OBS, alpha=ALPHA.BAND)

    # SNOWPACK subsets — decreasing hue encodes subset scope
    ax.bar(x, y, color=SUBSET_COLOR.ALL,        label="All stations")
    ax.bar(x, y, color=SUBSET_COLOR.RAIN_GAUGE, label="Rain-gauge subset")
    ax.bar(x, y, color=SUBSET_COLOR.BELOW_2000, label="< 2000 m subset")

    # Optimizer variants — DE hatched, Nelder-Mead solid
    ax.bar(x, y, color=col, hatch=HATCH.DE, label="Differential Evolution")
    ax.bar(x, y, color=col, hatch=HATCH.NM, label="Nelder-Mead")

    # Variable linestyles — SWE and rho_bulk are primary (most legible)
    ax.plot(t, swe,      color=c, linestyle=LS.SWE)
    ax.plot(t, rho_bulk, color=c, linestyle=LS.RHO_BULK)
    ax.plot(t, hs,       color=c, linestyle=LS.HS)
    ax.plot(t, hnw,      color=c, linestyle=LS.HNW)

Color palette
-------------
All colours are taken from the Okabe-Ito palette or close derivatives —
verified distinguishable under protanopia, deuteranopia, and tritanopia.

    OBS        near-black     authoritative observations  (CRNS, manual pits)
    WINKLER    vermilion      Winkler dataset
    MAGNUSSON  steel-blue     Magnusson / Mag25 dataset
    SNOWPACK   rose-purple    SNOWPACK dataset (full hue = all-stations reference)
    DSNOW      amber-orange   ΔSnow model output
    HS2SWE     bluish-green   HS2SWE model output
    HS         mid-grey       Snow height / hs (auxiliary)

SNOWPACK subset colours  (decreasing hue / lightness)
------------------------------------------------------
    SUBSET_COLOR.ALL         #8E2A84   full hue    — all stations
    SUBSET_COLOR.RAIN_GAUGE  #B06AA9   ~70 % hue   — rain-gauge subset
    SUBSET_COLOR.BELOW_2000  #D2AACE   ~40 % hue   — altitude < 2000 m subset

Optimizer hatch variants
------------------------
    HATCH.NM     ""      solid (no decoration)  — Nelder-Mead
    HATCH.DE     ////    forward diagonal       — Differential Evolution

Variable linestyles  (LS)
--------------------------
Priority pair (most legible — thick, high-contrast):
    LS.SWE       "-"     solid          — SWE
    LS.RHO_BULK  "-."    dash-dot       — bulk snow density

Secondary pair (thinner, less prominent):
    LS.HS        "--"    dashed         — snow height
    LS.HNW       ":"     dotted         — new-snow water equivalent
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# ── Colours ──────────────────────────────────────────────────────────────────

class C:
    OBS       = "#222222"   # near-black   — observed reference
    WINKLER   = "#D55E00"   # vermilion    — Winkler dataset
    MAGNUSSON = "#0072B2"   # steel-blue   — Mag25 dataset
    SNOWPACK  = "#8E2A84"   # deep plum    — SNOWPACK dataset
    DSNOW     = "#E69F00"   # amber        — ΔSnow model
    HS2SWE    = "#009E73"   # blue-green   — HS2SWE model
    HS        = "#999999"   # mid-grey     — snow height (hs), auxiliary
    NM        = "#35A6DF"  #
    DE        = "#882255"   # dark grey    — Differential Evolution optimizer



    # Convenience list for cycling through models / datasets
    CYCLE = [DSNOW, HS2SWE, MAGNUSSON, SNOWPACK, WINKLER, OBS, HS]


class ALPHA:
    BAND   = 0.18   # uncertainty / percentile band
    FILL   = 0.30   # lighter area fills
    MARKER = 0.70   # semi-transparent scatter points


class FIG:
    """Standard figure sizes and save settings for the mt_dsnow project.

    All figures are saved as PNG at 350 dpi — no PDF output.

    Usage
    -----
        fig, ax  = plt.subplots(1, 1, figsize=FIG.SINGLE)
        fig, axs = plt.subplots(3, 1, figsize=FIG.TALL3)
        fig, axs = plt.subplots(3, 3, figsize=FIG.GRID3x3)
        fig, axs = plt.subplots(1, 4, figsize=FIG.BOX4)
        fig.savefig(path, **FIG.SAVE)
    """
    SINGLE  = (10,  5)   # single panel or simple boxplot
    TALL3   = (10,  9)   # 3-panel tall (e.g. seasonal climatology column)
    GRID3x3 = (12,  9)   # 3×3 panel grid
    WIDE2x2 = (15, 10)   # 2×2 boxplot / comparison grid
    BOX4    = (12,  5)   # 1×4 side-by-side boxplots
    SCAT1   = (6,  6)       # single scatter plots

    DPI  = 350
    SAVE = dict(dpi=DPI, bbox_inches="tight", format="png")


class LS:
    """Linestyles for the four primary snow variables.

    SWE and rho_bulk are the priority pair — solid and dash-dot are the two
    most legible patterns.  HS and HNW use dashed / dotted and are intentionally
    less prominent so they don't compete with the key variables.
    """
    SWE      = "-"    # solid       — SWE           (primary)
    RHO_BULK = "-."   # dash-dot    — bulk density   (primary)
    HS       = "--"   # dashed      — snow height    (secondary)
    HNW      = ":"    # dotted      — new-snow HNW   (secondary)


class SUBSET_COLOR:
    """Face colours for the three SNOWPACK calibration subsets.

    Decreasing lightness encodes subset scope; all three are tints of
    C.SNOWPACK (#8E2A84) blended toward white at 100 %, 70 %, and 40 %.
    """
    ALL        = "#8E2A84"   # full hue  — all stations
    RAIN_GAUGE = "#B06AA9"   # ~70 % hue — rain-gauge subset
    BELOW_2000 = "#D2AACE"   # ~40 % hue — altitude < 2000 m subset

    CYCLE = [ALL, RAIN_GAUGE, BELOW_2000]


class HATCH:
    """Hatch patterns for calibration optimizer variants.

    Nelder-Mead is rendered as a solid fill; Differential Evolution gets a
    forward-diagonal hatch so both optimizers remain distinguishable when
    subset colours are already encoding the station scope.
    """
    NM = ""       # solid — Nelder-Mead
    DE = "////"   # forward diagonal — Differential Evolution

    CYCLE = [NM, DE]


# ── rcParams style ────────────────────────────────────────────────────────────

_STYLE: dict = {
    # Figure — slightly smaller default, tighter outer padding
    "figure.figsize":        (9, 4.5),
    "figure.dpi":            100,
    "figure.facecolor":      "white",
    "savefig.dpi":           350,
    "savefig.bbox":          "tight",
    "savefig.pad_inches":    0.08,
    "savefig.facecolor":     "white",

    # Font — softer sizes, lighter weights, muted colours
    "font.family":           "sans-serif",
    "font.size":             10,
    "axes.titlesize":        12,
    "axes.titleweight":      "normal",
    "axes.titlecolor":       "#333333",
    "figure.titlesize":      16,
    "figure.titleweight":    "semibold",
    "axes.labelsize":        10,
    "axes.labelcolor":       "#444444",
    "xtick.labelsize":       9,
    "ytick.labelsize":       9,
    "legend.fontsize":       9,
    "legend.title_fontsize": 9,

    # Axes — lighter borders
    "axes.facecolor":        "white",
    "axes.edgecolor":        "#BBBBBB",
    "axes.linewidth":        1.0,
    "axes.spines.top":       True,
    "axes.spines.right":     True,
    "axes.prop_cycle":       mpl.cycler(color=C.CYCLE),

    # Subplot spacing — tighter defaults for multi-panel figures
    "figure.subplot.wspace": 0.22,
    "figure.subplot.hspace": 0.30,

    # Boxplot defaults
    "boxplot.medianprops.color": "black",
    "boxplot.medianprops.linewidth": 1.4,

    # Grid — subtle, only horizontal by default
    "axes.grid":             True,
    "axes.grid.axis":        "y",
    "grid.color":            "#EBEBEB",
    "grid.linewidth":        0.6,
    "grid.linestyle":        "-",

    # Lines & markers
    "lines.linewidth":       1.6,
    "lines.markersize":      5,

    # Legend
    "legend.framealpha":     0.85,
    "legend.edgecolor":      "#DDDDDD",
    "legend.frameon":        True,

    # Ticks — muted colour
    "xtick.direction":       "out",
    "ytick.direction":       "out",
    "xtick.color":           "#777777",
    "ytick.color":           "#777777",
    "xtick.major.size":      3.5,
    "ytick.major.size":      3.5,

    # Date formatting
    "date.autoformatter.month": "%b",
    "date.autoformatter.year":  "%Y",
}


def apply_style() -> None:
    """Apply the project-wide rcParams. Call once per script or notebook."""
    mpl.rcParams.update(_STYLE)


# ── Quick reference ────────────────────────────────────────────────────────────

def print_palette() -> None:
    """Print a colour swatch to the terminal / notebook for reference."""
    names  = ["OBS",  "WINKLER",  "MAGNUSSON",  "DSNOW",  "HS2SWE",  "HS"]
    colors = [C.OBS,  C.WINKLER,  C.MAGNUSSON,  C.DSNOW,  C.HS2SWE,  C.HS]

    fig, ax = plt.subplots(figsize=(len(names) * 1.4, 1.2))
    for i, (name, col) in enumerate(zip(names, colors)):
        ax.add_patch(plt.Rectangle((i, 0), 0.9, 0.8, color=col))
        ax.text(i + 0.45, -0.15, name, ha="center", va="top", fontsize=9)
        ax.text(i + 0.45, 0.9,   col,  ha="center", va="bottom", fontsize=7,
                color="#555555")
    ax.set_xlim(0, len(names))
    ax.set_ylim(-0.3, 1.1)
    ax.axis("off")
    ax.set_title("mt_dsnow colour palette", fontsize=10, pad=4)
    plt.tight_layout()
    plt.show()


def print_subset_colors() -> None:
    """Print a colour swatch for the three SNOWPACK subset tints."""
    labels = ["All stations",       "Rain-gauge subset",  "< 2000 m subset"]
    colors = [SUBSET_COLOR.ALL,     SUBSET_COLOR.RAIN_GAUGE, SUBSET_COLOR.BELOW_2000]

    fig, ax = plt.subplots(figsize=(len(labels) * 2.0, 1.4))
    for i, (label, col) in enumerate(zip(labels, colors)):
        ax.add_patch(plt.Rectangle((i * 1.9, 0), 1.7, 0.9, color=col))
        ax.text(i * 1.9 + 0.85, -0.15, label, ha="center", va="top", fontsize=9)
        ax.text(i * 1.9 + 0.85, 1.0,   col,   ha="center", va="bottom",
                fontsize=7, color="#555555")
    ax.set_xlim(0, len(labels) * 1.9)
    ax.set_ylim(-0.4, 1.2)
    ax.axis("off")
    ax.set_title("SNOWPACK subset colours  (decreasing hue)", fontsize=10, pad=4)
    plt.tight_layout()
    plt.show()


def print_hatches() -> None:
    """Print a hatch swatch for the two optimizer variants (DE vs NM)."""
    labels  = ["Nelder-Mead (NM)",  "Differential Evolution (DE)"]
    hatches = [HATCH.NM,            HATCH.DE]

    fig, ax = plt.subplots(figsize=(len(labels) * 2.5, 1.4))
    for i, (label, hatch) in enumerate(zip(labels, hatches)):
        ax.add_patch(plt.Rectangle(
            (i * 2.2, 0), 2.0, 0.9,
            facecolor=C.SNOWPACK, edgecolor="white",
            hatch=hatch, linewidth=0.0,
        ))
        ax.text(i * 2.2 + 1.0, -0.15, label, ha="center", va="top", fontsize=9)
        hatch_str = repr(hatch) if hatch else '""  (solid)'
        ax.text(i * 2.2 + 1.0, 1.0, hatch_str, ha="center", va="bottom",
                fontsize=7, color="#555555")
    ax.set_xlim(0, len(labels) * 2.2)
    ax.set_ylim(-0.4, 1.2)
    ax.axis("off")
    ax.set_title("Optimizer hatch variants  (DE hatched, NM solid)", fontsize=10, pad=4)
    plt.tight_layout()
    plt.show()


def dataset_color(dataset: str, subset: str = "ALL") -> str:
    """Return the canonical face colour for a dataset.

    Parameters
    ----------
    dataset : str
        Case-insensitive name: 'SNOWPACK', 'Win21', 'Magnusson'/'Mag25',
        'dSnow', 'HS2SWE'.
    subset : str, default 'ALL'
        SNOWPACK subset — 'ALL', 'RAIN_GAUGE', or 'BELOW_2000'.
        Ignored for all other datasets.

    Returns
    -------
    str  — hex colour

    Examples
    --------
    >>> dataset_color('SNOWPACK')               # full hue
    '#8E2A84'
    >>> dataset_color('SNOWPACK', 'RAIN_GAUGE') # ~70 % hue
    '#B06AA9'
    >>> dataset_color('Win21')
    '#D55E00'
    """
    key = dataset.strip().upper().replace('-', '_').replace(' ', '_')
    if 'SNOWPACK' in key:
        return getattr(SUBSET_COLOR, subset.upper(), SUBSET_COLOR.ALL)
    mapping = {
        'WIN21':      C.WINKLER,
        'WINKLER':    C.WINKLER,
        'MAGNUSSON':  C.MAGNUSSON,
        'MAG25':      C.MAGNUSSON,
        'DSNOW':      C.DSNOW,
        'DELTANOW':   C.DSNOW,
        'HS2SWE':     C.HS2SWE,
        'HS':         C.HS,
        'OBS':        C.OBS,
    }
    return mapping.get(key, '#888888')


def build_palette(datasets, subset_map=None) -> dict:
    """Build a ``{dataset_name: colour}`` dict for seaborn / matplotlib.

    Parameters
    ----------
    datasets : list of str
        Dataset names as they appear in a DataFrame column.
    subset_map : dict, optional
        ``{dataset_name: subset_key}`` overrides for SNOWPACK entries.
        Keys are dataset names; values are 'ALL', 'RAIN_GAUGE', or
        'BELOW_2000'.  Entries not listed default to 'ALL'.

    Examples
    --------
    >>> build_palette(['SNOWPACK', 'Win21'])
    {'SNOWPACK': '#8E2A84', 'Win21': '#D55E00'}

    >>> build_palette(
    ...     ['SNOWPACK_all', 'SNOWPACK_rg', 'Win21'],
    ...     subset_map={'SNOWPACK_rg': 'RAIN_GAUGE'},
    ... )
    {'SNOWPACK_all': '#8E2A84', 'SNOWPACK_rg': '#B06AA9', 'Win21': '#D55E00'}
    """
    subset_map = subset_map or {}
    return {ds: dataset_color(ds, subset_map.get(ds, 'ALL')) for ds in datasets}


def get_dataset_mean_legend_handles():
    """Legend handles for Win21, Mag25, SNOWPACK, and mean marker.

    Use this in multi-panel comparison plots when you want a consistent
    legend without hatch encoding.
    """
    return [
        Patch(facecolor=C.WINKLER, edgecolor="#1A1A1A", label="Win21"),
        Patch(facecolor=C.MAGNUSSON, edgecolor="#1A1A1A", label="Mag25"),
        Patch(facecolor=C.SNOWPACK, edgecolor="#1A1A1A", label="SNOWPACK"),
        Line2D(
            [0],
            [0],
            marker="D",
            color="none",
            markerfacecolor="#DDDDDD",
            markeredgecolor="#555555",
            markeredgewidth=0.8,
            markersize=6,
            label="Mean marker",
        ),
    ]


def add_subplot_labels(axs, start="a", x=0.01, y=0.99):
    """Add subplot labels (a, b, c, ...) to one or more axes.

    Parameters
    ----------
    axs : matplotlib Axes or array-like of Axes
        Single axis or collection (e.g. from plt.subplots).
    start : str, default "a"
        Starting lowercase letter for labels.
    x, y : float
        Axes-relative text position (0..1). Defaults place label at top-left.
    """
    axes = axs.ravel() if hasattr(axs, "ravel") else [axs]
    start_ord = ord(start)
    for i, ax in enumerate(axes):
        ax.text(
            x,
            y,
            f"{chr(start_ord + i)})",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            fontweight="normal",
            color="#777777",
        )
