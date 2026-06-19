import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import pandas as pd
import xarray as xr


def _filter_season(df, full_season=False):
    """Filter DataFrame to snow season (Nov 1 – Apr 30) unless full_season=True."""

    df = _set_datime_index(df)

    if full_season:
        return df.copy()

    else:
        mask = (df.index.month >= 11) | (df.index.month <= 4)
        return df.loc[mask].copy()

def _set_datime_index(df):
    """Ensure DataFrame has a DatetimeIndex."""

    df.index = pd.to_datetime(df.time)
    return df


def _calculate_metrics(x, y):
    """
    Calculate validation metrics.
    x = observed
    y = modeled
    """

    residuals = y - x

    rmse = np.sqrt(np.mean(residuals**2))
    bias = np.mean(residuals)
    pbias = np.sum(residuals) / np.sum(x) if np.sum(x) != 0 else np.nan

    ss_res = np.sum((x - y)**2)
    ss_tot = np.sum((x - np.mean(x))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

    return {
        "RMSE": rmse,
        "Bias": bias,
        "Rel_BIAS": pbias,
        "R2": r2
    }


def make_model_cmap(color, name="_model_cmap"):
    """Build a near-white → *color* LinearSegmentedColormap for density plots.

    Parameters
    ----------
    color : str
        Any matplotlib-compatible colour string, typically a project hex
        constant such as ``C.HS2SWE`` (``"#009E73"``) or ``C.DSNOW``
        (``"#E69F00"``).
    name : str, optional
        Internal name registered with matplotlib (default ``"_model_cmap"``).

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        Ramp from ``#F0F0F0`` (sparse bins) → *color* (dense bins).

    Example
    -------
    >>> from plot_style import C
    >>> cmap = val_helper.make_model_cmap(C.HS2SWE)
    """
    return LinearSegmentedColormap.from_list(name, ["#F0F0F0", color], N=256)


def _resolve_cmap(cmap):
    """Return a Colormap object from a name, hex colour, or existing Colormap.

    Behaviour
    ---------
    - Named matplotlib colormap (e.g. ``"viridis"``): returned directly.
    - Hex or CSS colour string (e.g. ``"#009E73"``, ``"steelblue"``):
      a two-stop LinearSegmentedColormap from near-white → that colour is
      created (near-white = sparse bins, deep colour = dense bins).
    - Colormap object: passed through unchanged.

    This lets notebooks do ``cmap = C.HS2SWE`` (a hex string) and get a
    single-hue ramp that matches the project colour scheme.
    """
    if not isinstance(cmap, str):
        return cmap  # already a Colormap object

    # Named matplotlib colormap → use directly
    try:
        return plt.get_cmap(cmap)
    except (ValueError, KeyError):
        pass

    # Treat as a plain colour → build a near-white → colour ramp
    if mcolors.is_color_like(cmap):
        return LinearSegmentedColormap.from_list("_model_cmap", ["#F0F0F0", cmap], N=256)

    # Last resort: let matplotlib decide (will raise a clear error if invalid)
    return cmap


def _nice_vmax(n):
    """Round len(x)/10 up to the nearest power of 10 for clean colorbar ticks."""
    raw = max(1, n / 10)
    return float(10 ** int(np.ceil(np.log10(max(raw, 1.01)))))


def _cb_ticks(vmax):
    """Return power-of-10 tick values from 1 up to vmax (inclusive)."""
    n_decades = int(round(np.log10(vmax)))
    return [10**i for i in range(n_decades + 1)]


def _plot_validation(x, y, stats, model_name, lim, xlabel, ylabel, cmap="viridis"):
    """Generic single-panel validation density scatter."""

    vmax = _nice_vmax(len(x))

    fig, ax = plt.subplots(figsize=(8, 7))

    h, xedges, yedges, img = ax.hist2d(
        x, y,
        bins=50,
        range=[lim, lim],
        norm=LogNorm(vmin=1, vmax=vmax),
        cmap=_resolve_cmap(cmap),
    )

    # Equal aspect first so make_axes_locatable sees the final square size
    ax.set_aspect("equal")

    # Colorbar locked to the square data area, not the full bounding box
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cb = fig.colorbar(img, cax=cax)
    cb.set_label("Number of observations")

    ticks = _cb_ticks(vmax)
    cb.set_ticks(ticks)
    cb.set_ticklabels([str(t) for t in ticks])

    ax.plot(lim, lim, "--", color="gray", linewidth=1.3)

    ticks_xy = np.linspace(lim[0], lim[1], 5)
    ax.set_xticks(ticks_xy)
    ax.set_yticks(ticks_xy)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(model_name)

    textstr = (
        f"$R^2$: {stats['R2']:.2f}\n"
        f"Bias: {stats['Bias']:.2f}\n"
        f"RMSE: {stats['RMSE']:.1f}\n"
        f"Rel_BIAS: {stats['Rel_BIAS']:.1%}\n"
    )

    ax.text(
        0.03, 0.97, textstr,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.grid(False)
    fig.tight_layout()


def _plot_validation_ax(ax, x, y, stats, title, lim, xlabel, ylabel, fig, cmap="viridis"):
    """Draw a density validation scatter into an existing Axes."""

    vmax = _nice_vmax(len(x))

    h, xedges, yedges, img = ax.hist2d(
        x, y,
        bins=50,
        range=[lim, lim],
        norm=LogNorm(vmin=1, vmax=vmax),
        cmap=_resolve_cmap(cmap),
    )

    # Equal aspect first so make_axes_locatable sees the final square size
    ax.set_aspect("equal")

    # Colorbar locked to the square data area, not the full bounding box
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cb = fig.colorbar(img, cax=cax)
    cb.set_label("Number of observations")

    ticks = _cb_ticks(vmax)
    cb.set_ticks(ticks)
    cb.set_ticklabels([str(t) for t in ticks])

    ax.plot(lim, lim, "--", color="gray", linewidth=1.3)

    ticks_xy = np.linspace(lim[0], lim[1], 5)
    ax.set_xticks(ticks_xy)
    ax.set_yticks(ticks_xy)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    textstr = (
        f"$R^2$: {stats['R2']:.2f}\n"
        f"Bias: {stats['Bias']:.2f}\n"
        f"RMSE: {stats['RMSE']:.1f}\n"
        f"Rel_BIAS: {stats['Rel_BIAS']:.1%}\n"
        f"N: {stats['N']}"
    )

    ax.text(
        0.03, 0.97, textstr,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.grid(False)


def validate_hnw_mag25(df,
                       model_name,
                       obs_col="HNW_obs",
                       mod_col="HNW_mod",
                       save_dir=None,
                       filename="hnw_validation.png",
                       full_season=False,
                       drop_weisfluh_joch=True,
                       ax=None,
                       cmap="viridis"):
    if drop_weisfluh_joch:
        df = df[df["station"] != "Weisfluh_Joch"].copy()

    df = _filter_season(df, full_season)

    df_valid = df.dropna(subset=[obs_col, mod_col])

    # Caution: great difference between >= 0 and > 0
    df_valid = df_valid[df_valid[obs_col] >= 0]

    df_valid = df_valid[
        np.isfinite(df_valid[obs_col]) &
        np.isfinite(df_valid[mod_col])
    ]

    x = df_valid[obs_col].values
    y = df_valid[mod_col].values

    stats = _calculate_metrics(x, y)
    stats["N"] = len(df_valid)

    print(stats)

    _plot_validation(
        x,
        y,
        stats,
        model_name,
        lim=[0, 100],
        xlabel="Observed HNW (mm)",
        ylabel="Modeled HNW (mm)",
        cmap=cmap,
    )

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, dpi=350, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()

    return stats


def validate_swe_mag25(df,
                       model_name,
                       obs_col="SWE_obs",
                       mod_col="SWE_mod",
                       save_dir=None,
                       filename="swe_validation.png",
                       full_season=False,
                       drop_weisfluh_joch=True,
                       cmap="viridis"):

    if drop_weisfluh_joch:
        df = df[df["station"] != "Weisfluh_Joch"].copy()

    df = _filter_season(df, full_season)

    df_valid = df.dropna(subset=[obs_col, mod_col])

    df_valid = df_valid[df_valid[obs_col] >= 0]

    print(f"Number of valid observations after filtering: {len(df_valid)}")

    df_valid = df_valid[
        np.isfinite(df_valid[obs_col]) &
        np.isfinite(df_valid[mod_col])
    ]

    x = df_valid[obs_col].values
    y = df_valid[mod_col].values

    stats = _calculate_metrics(x, y)
    stats["N"] = len(df_valid)

    print(stats)

    _plot_validation(
        x,
        y,
        stats,
        model_name,
        lim=[0, 1000],
        xlabel="Observed SWE (mm)",
        ylabel="Modeled SWE (mm)",
        cmap=cmap,
    )

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, dpi=350, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()

    return stats


def validate_hnw_swe_combined(hnw_df, swe_df, model_name,
                               params=None,
                               hnw_obs_col="HNW_obs", hnw_mod_col="HNW_mod",
                               swe_obs_col="SWE_obs", swe_mod_col="SWE_mod",
                               full_season=False, drop_weisfluh_joch=True,
                               save_dir=None, filename="hnw_swe_validation_combined.png",
                               cmap="viridis"):
    """
    Plot HNW and SWE validation side-by-side in one figure.

    Parameters
    ----------
    hnw_df : pd.DataFrame  — must contain hnw_obs_col, hnw_mod_col, station, time
    swe_df : pd.DataFrame  — must contain swe_obs_col, swe_mod_col, station, time
    model_name : str       — shown as figure suptitle
    params : dict, optional
        Model parameters to annotate in the figure.
    full_season : bool     — if False, restrict to Nov–Apr
    drop_weisfluh_joch : bool
    save_dir : str or Path, optional
    filename : str
    cmap : str or Colormap
        Controls the 2-D histogram colour scheme for both panels.
        - Named matplotlib colormap (e.g. ``"viridis"``, ``"plasma"``): used directly.
        - Hex / CSS colour string (e.g. ``C.HS2SWE``, ``C.DSNOW``): a single-hue
          ramp from near-white (few obs) to that colour (many obs) is generated
          automatically, tying each plot to its model brand colour.

    Returns
    -------
    dict with keys "HNW" and "SWE", each a stats dict
    (RMSE, Bias, Rel_BIAS, R2, N).
    """

    # ── HNW preparation ───────────────────────────────────────────────────────
    hnw = hnw_df.copy()
    if drop_weisfluh_joch:
        hnw = hnw[hnw["station"] != "Weisfluh_Joch"]
    hnw = _filter_season(hnw, full_season)
    hnw = hnw.dropna(subset=[hnw_obs_col, hnw_mod_col])
    hnw = hnw[hnw[hnw_obs_col] >= 0]
    hnw = hnw[np.isfinite(hnw[hnw_obs_col]) & np.isfinite(hnw[hnw_mod_col])]

    x_hnw = hnw[hnw_obs_col].values
    y_hnw = hnw[hnw_mod_col].values
    stats_hnw = _calculate_metrics(x_hnw, y_hnw)
    stats_hnw["N"] = len(hnw)

    # ── SWE preparation ───────────────────────────────────────────────────────
    swe = swe_df.copy()
    if drop_weisfluh_joch:
        swe = swe[swe["station"] != "Weisfluh_Joch"]
    swe = _filter_season(swe, full_season=True)   # SWE uses full year
    swe = swe.dropna(subset=[swe_obs_col, swe_mod_col])
    swe = swe[swe[swe_obs_col] >= 0]
    swe = swe[np.isfinite(swe[swe_obs_col]) & np.isfinite(swe[swe_mod_col])]

    x_swe = swe[swe_obs_col].values
    y_swe = swe[swe_mod_col].values
    stats_swe = _calculate_metrics(x_swe, y_swe)
    stats_swe["N"] = len(swe)

    print("HNW stats:", stats_hnw)
    print("SWE stats:", stats_swe)

    # ── Build params annotation string ────────────────────────────────────────
    def _fmt_param(v):
        """Format a scalar as z.zz×10^x, stripping leading zeros from exponent."""
        if v == 0:
            return "0"
        s = f"{v:.2e}"                          # e.g. "1.11e+02"
        mantissa, exp_part = s.split("e")
        sign = exp_part[0]                       # '+' or '-'
        exp_int = int(exp_part[1:])              # strip leading zeros
        exp_str = f"{sign}{exp_int}" if sign == "-" else str(exp_int)
        return f"{mantissa}×10^{exp_str}"

    param_str = ""
    if params is not None:
        parts = [f"{k}={_fmt_param(v)}" for k, v in params.items()]
        n_per_line = 5
        lines = [" | ".join(parts[i:i + n_per_line]) for i in range(0, len(parts), n_per_line)]
        param_str = "\n".join(lines)

    # ── Figure ────────────────────────────────────────────────────────────────
    # (14, 7): each panel ~6" wide after margins/colorbar, equal-aspect clamps
    # the square data area to ~6"×6" — height matches, no wasted whitespace.
    # constrained_layout is OFF; make_axes_locatable is incompatible with it.
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    _plot_validation_ax(
        ax=axes[0], x=x_hnw, y=y_hnw, stats=stats_hnw,
        title="HNW",
        lim=[0, 110],
        xlabel="Observed HNW (mm)",
        ylabel="Modeled HNW (mm)",
        fig=fig,
        cmap=cmap,
    )

    _plot_validation_ax(
        ax=axes[1], x=x_swe, y=y_swe, stats=stats_swe,
        title="SWE",
        lim=[0, 1000],
        xlabel="Observed SWE (mm)",
        ylabel="Modeled SWE (mm)",
        fig=fig,
        cmap=cmap,
    )

    fig.suptitle(model_name)
    fig.tight_layout()

    # Param string sits below the subplots, clear of the x-axis labels
    if param_str:
        fig.text(0.5, -0.02, param_str,
                 ha="center", va="top",
                 fontsize=8, color="#444444")

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        fig.savefig(save_path, dpi=350, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()

    return {"HNW": stats_hnw, "SWE": stats_swe}


# ─────────────────────────────────────────────────────────────────────────────
# Metrics-only variants (no plots)
# ─────────────────────────────────────────────────────────────────────────────

def compute_metrics(df,
                    obs_col,
                    mod_col,
                    full_season=False,
                    drop_weisfluh_joch=True,
                    obs_min=0):
    """
    Calculate RMSE, Bias, Rel_BIAS, R2, and N for one variable.

    Applies the same filtering pipeline as the plotting functions
    (snow-season mask, non-negative observed values, finite pairs)
    but returns only the stats dict — no plot.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``obs_col``, ``mod_col``, ``station``, and ``time``.
    obs_col : str
        Column name for observed values.
    mod_col : str
        Column name for modelled values.
    full_season : bool, default False
        If False, restrict to Nov 1 – Apr 30.
    drop_weisfluh_joch : bool, default True
        Exclude rows where station == "Weisfluh_Joch".
    obs_min : float, default 0
        Drop rows where observed < obs_min (mirrors the >= 0 filter
        used in the plotting functions; pass None to skip).

    Returns
    -------
    dict with keys: RMSE, Bias, Rel_BIAS, R2, N
    """
    d = df.copy()

    if drop_weisfluh_joch and "station" in d.columns:
        d = d[d["station"] != "Weisfluh_Joch"]

    d = _filter_season(d, full_season)
    d = d.dropna(subset=[obs_col, mod_col])

    if obs_min is not None:
        d = d[d[obs_col] >= obs_min]

    d = d[np.isfinite(d[obs_col]) & np.isfinite(d[mod_col])]

    x = d[obs_col].values
    y = d[mod_col].values

    stats = _calculate_metrics(x, y)
    stats["N"] = len(d)

    return stats


def compute_metrics_hnw_swe(hnw_df, swe_df,
                             hnw_obs_col="HNW_obs", hnw_mod_col="HNW_mod",
                             swe_obs_col="SWE_obs", swe_mod_col="SWE_mod",
                             full_season=False,
                             drop_weisfluh_joch=True):
    """
    Calculate metrics for HNW and SWE simultaneously — no plots.

    Mirrors the filtering in ``validate_hnw_swe_combined``:
    HNW uses the snow-season mask (Nov–Apr unless full_season=True);
    SWE always uses the full year.

    Parameters
    ----------
    hnw_df, swe_df : pd.DataFrame
    hnw_obs_col, hnw_mod_col : str
    swe_obs_col, swe_mod_col : str
    full_season : bool
        Applied to HNW only; SWE always uses the full year.
    drop_weisfluh_joch : bool

    Returns
    -------
    dict with keys "HNW" and "SWE", each a dict:
        {RMSE, Bias, Rel_BIAS, R2, N}
    """
    stats_hnw = compute_metrics(
        hnw_df,
        obs_col=hnw_obs_col,
        mod_col=hnw_mod_col,
        full_season=full_season,
        drop_weisfluh_joch=drop_weisfluh_joch,
    )

    stats_swe = compute_metrics(
        swe_df,
        obs_col=swe_obs_col,
        mod_col=swe_mod_col,
        full_season=True,          # SWE: always full year
        drop_weisfluh_joch=drop_weisfluh_joch,
    )

    return {"HNW": stats_hnw, "SWE": stats_swe}
