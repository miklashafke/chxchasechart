#!/usr/bin/env python3
"""Quantify a cycloheximide (CHX) chase from ImageJ measurement tables.

Takes two ImageJ Results tables -- one for the protein of interest, one for the
PGK1 loading control -- optionally subtracts a background measurement, normalizes
the chosen measurement (default ``Mean``) to PGK1, expresses every time point as
a percentage of the first time point, and plots the decay curves.

The two tables must have the same number of rows, one row per measured band,
and the rows must be in the order the experiment was laid out (see ``--order``;
the default is condition -> time point -> replicate, condition varying slowest).

``--background`` takes a two-row table: row 1 is the background for the raw
measurement, row 2 the background for the PGK1 control. It is read with the same
measurement column and subtracted from every band of the matching table.

Example
-------
The dataset in ``example/`` is three conditions x four time points x three
biological replicates, with a measured PGK1 control and a background table.
Blue is WT (half-life 30 min), orange a completely stabilized mutant, and green
an intermediate mutant (half-life 60 min)::

    python chx_chase.py example/raw.csv example/pgk1_control.csv 3 -r 3 \
        --background example/background.tsv

Simplest possible call, taking the defaults (4 time points, 1 replicate,
0/30/60/90 min, no background)::

    python chx_chase.py raw.csv pgk1_control.csv 4

Both headed ImageJ exports (`` ,Area,Mean,Min,Max``) and raw headerless
copy-paste dumps are accepted, comma- or tab-separated.

Measure the bands in ImageJ and export the Results table; this script takes it
from there.
"""

from __future__ import annotations

import argparse
import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Categorical palette, fixed slot order (validated for line charts, light surface).
PALETTE = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#4a3aa7",  # violet
    "#e34948",  # red
]
LINESTYLES = ["-", "--", ":"]  # secondary encoding once the palette wraps

# Marks are specified in millimeters, as figure panels usually are.
POINTS_PER_MM = 72 / 25.4
LINE_WIDTH_MM = 0.5
MARKER_SIZE_MM = 1.5  # dot diameter, white outline included
MARKER_EDGE_MM = 0.2  # white ring separating the dot from the line behind it

SURFACE = "#ffffff"
INK = "#000000"
INK_SECONDARY = "#000000"
INK_MUTED = "#000000"
GRID = "#e6e6e6"
AXIS = "#000000"

# ImageJ's default measurement set, for tables pasted without a header row.
IMAGEJ_DEFAULT_COLUMNS = ["Area", "Mean", "Min", "Max"]

ORDER_NAMES = {"c": "condition", "t": "timepoint", "r": "replicate"}


class InputError(Exception):
    """Something is wrong with the files or the arguments the user supplied."""


# --------------------------------------------------------------------------- #
# reading ImageJ tables
# --------------------------------------------------------------------------- #

def _looks_numeric(text: object) -> bool:
    """True if ``text`` is a number, ignoring pandas' ``.1`` duplicate suffix."""
    stripped = re.sub(r"\.\d+$", "", str(text).strip())
    try:
        float(stripped)
    except ValueError:
        return False
    return True


def read_imagej_table(path: Path) -> pd.DataFrame:
    """Read an ImageJ Results table, with or without a header row."""
    if not path.is_file():
        raise InputError(f"no such file: {path}")
    if path.suffix.lower() in (".tif", ".tiff"):
        raise InputError(
            f"{path.name} is a blot image, not a measurement table. Measure the bands "
            "in ImageJ (Analyze > Measure) and export the Results table."
        )

    try:
        table = pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:  # pandas raises a zoo of parser errors
        raise InputError(f"could not parse {path.name}: {exc}") from exc

    # A headerless table makes pandas promote the first data row to column names.
    if all(_looks_numeric(name) for name in table.columns):
        table = pd.read_csv(path, sep=None, engine="python", header=None)
        if table.shape[1] >= len(IMAGEJ_DEFAULT_COLUMNS):
            # ImageJ measures Area, Mean, Min, Max by default and puts them last;
            # anything before them is a row number, label or slice index.
            leading = ["Index", "Label", "Slice"]
            extra = table.shape[1] - len(IMAGEJ_DEFAULT_COLUMNS)
            names = [leading[i] if i < len(leading) else f"Col{i + 1}" for i in range(extra)]
            table.columns = names + IMAGEJ_DEFAULT_COLUMNS
            print(
                f"note: {path.name} has no header; read its last four columns as "
                f"{', '.join(IMAGEJ_DEFAULT_COLUMNS)}. Export with a header row if they differ.",
                file=sys.stderr,
            )
        else:
            table.columns = ["Index"] + [f"Col{i}" for i in range(1, table.shape[1])]
            print(
                f"note: {path.name} has no header and only {table.shape[1]} columns; "
                f"columns named {list(table.columns)} -- pass --column to pick one.",
                file=sys.stderr,
            )

    table.columns = [str(name).strip() for name in table.columns]
    if table.columns[0] in ("", "Unnamed: 0"):
        table = table.rename(columns={table.columns[0]: "Index"})
    # ImageJ sometimes leaves a trailing separator, giving an all-empty column.
    table = table.dropna(axis="columns", how="all")

    if table.empty:
        raise InputError(f"{path.name} contains no measurements")
    return table


def pick_measurement(table: pd.DataFrame, column: str, path: Path) -> np.ndarray:
    """Pull one measurement column out of a table, case-insensitively."""
    lookup = {name.lower(): name for name in table.columns}
    key = column.strip().lower()
    if key not in lookup:
        raise InputError(
            f"{path.name} has no column {column!r}; available: {', '.join(table.columns)}"
        )

    values = pd.to_numeric(table[lookup[key]], errors="coerce")
    if values.isna().any():
        bad = (values.index[values.isna()] + 1).tolist()
        raise InputError(f"{path.name}, column {column!r}: non-numeric value(s) in row(s) {bad}")
    return values.to_numpy(dtype=float)


# --------------------------------------------------------------------------- #
# layout and normalization
# --------------------------------------------------------------------------- #

def build_layout(n_conditions: int, n_timepoints: int, n_replicates: int, order: str) -> pd.DataFrame:
    """Map row position in the table onto (condition, timepoint, replicate)."""
    sizes = {"c": n_conditions, "t": n_timepoints, "r": n_replicates}
    combos = itertools.product(*(range(sizes[axis]) for axis in order))
    layout = pd.DataFrame(combos, columns=[ORDER_NAMES[axis] for axis in order])
    return layout[["condition", "timepoint", "replicate"]]


def subtract_background(
    signal: np.ndarray, control: np.ndarray, path: Path, column: str
) -> tuple[np.ndarray, np.ndarray, tuple[float, float]]:
    """Subtract a two-row background table: row 1 from raw, row 2 from PGK1."""
    table = read_imagej_table(path)
    if len(table) != 2:
        raise InputError(
            f"{path.name} has {len(table)} rows; the background table needs exactly 2 "
            "(row 1 = raw, row 2 = PGK1 control)"
        )

    background = pick_measurement(table, column, path)
    signal = signal - background[0]
    control = control - background[1]

    if (control <= 0).any():
        rows = (np.flatnonzero(control <= 0) + 1).tolist()
        raise InputError(
            f"PGK1 control is <= 0 after subtracting background {background[1]:g} "
            f"in row(s) {rows}; check the background ROI"
        )
    if (signal <= 0).any():
        rows = (np.flatnonzero(signal <= 0) + 1).tolist()
        print(
            f"warning: raw signal is <= 0 after subtracting background {background[0]:g} "
            f"in row(s) {rows}; those points are at or below background.",
            file=sys.stderr,
        )
    return signal, control, (float(background[0]), float(background[1]))


def normalize(
    layout: pd.DataFrame, signal: np.ndarray, control: np.ndarray
) -> pd.DataFrame:
    """Normalize to PGK1, then to the first time point of each series."""
    tidy = layout.copy()
    tidy["signal"] = signal
    tidy["control"] = control

    if (tidy["control"] == 0).any():
        rows = (tidy.index[tidy["control"] == 0] + 1).tolist()
        raise InputError(f"PGK1 control is zero in row(s) {rows}; cannot normalize")
    tidy["ratio"] = tidy["signal"] / tidy["control"]

    first = (
        tidy.loc[tidy["timepoint"] == 0, ["condition", "replicate", "ratio"]]
        .rename(columns={"ratio": "ratio_t0"})
    )
    tidy = tidy.merge(first, on=["condition", "replicate"], how="left")

    if (tidy["ratio_t0"] == 0).any():
        bad = tidy.loc[tidy["ratio_t0"] == 0, ["condition", "replicate"]].drop_duplicates()
        pairs = ", ".join(
            f"condition {c + 1}/replicate {r + 1}" for c, r in bad.itertuples(index=False)
        )
        raise InputError(f"first time point is zero for {pairs}; cannot normalize to it")

    tidy["percent"] = 100 * tidy["ratio"] / tidy["ratio_t0"]
    return tidy.drop(columns="ratio_t0")


def summarize(tidy: pd.DataFrame) -> pd.DataFrame:
    """Mean and (when there is more than one replicate) SD per condition/time point."""
    summary = tidy.groupby(["condition", "timepoint"], as_index=False).agg(
        mean=("percent", "mean"),
        sd=("percent", "std"),  # ddof=1 -> NaN for a single replicate
        n=("percent", "size"),
    )
    return summary.sort_values(["condition", "timepoint"], ignore_index=True)


# --------------------------------------------------------------------------- #
# plotting
# --------------------------------------------------------------------------- #

def plot_chase(
    summary: pd.DataFrame,
    condition_names: list[str],
    times: np.ndarray,
    time_labels: list[str],
    xlabel: str,
    ylabel: str,
    title: str,
    show_errorbars: bool,
    show_legend: bool,
    figsize_cm: tuple[float, float],
    legend_loc: str,
    yticks: list[float],
    fontsize: float,
    font: str,
) -> "object":
    """One 1x1 axes, one line per condition, points connected."""
    from matplotlib import pyplot as plt

    set_font(font)

    line_width = LINE_WIDTH_MM * POINTS_PER_MM
    marker_edge = MARKER_EDGE_MM * POINTS_PER_MM
    # The stroke straddles the marker path, so face + edge = MARKER_SIZE_MM.
    marker_size = (MARKER_SIZE_MM - MARKER_EDGE_MM) * POINTS_PER_MM

    fig, ax = plt.subplots(figsize=(figsize_cm[0] / 2.54, figsize_cm[1] / 2.54),
                           facecolor=SURFACE, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=0, hspace=0)
    ax.set_facecolor(SURFACE)

    for index, name in enumerate(condition_names):
        series = summary[summary["condition"] == index].sort_values("timepoint")
        color = PALETTE[index % len(PALETTE)]
        style = LINESTYLES[(index // len(PALETTE)) % len(LINESTYLES)]
        errors = series["sd"].to_numpy(dtype=float) if show_errorbars else None

        ax.errorbar(
            times,
            series["mean"].to_numpy(dtype=float),
            yerr=errors,
            label=name,
            color=color,
            linestyle=style,
            linewidth=line_width,
            marker="o",
            markersize=marker_size,
            markeredgecolor=SURFACE,
            markeredgewidth=marker_edge,
            capsize=marker_size / 2,
            elinewidth=line_width,
            zorder=3 + index,
        )

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(line_width)

    ax.set_xticks(times, labels=time_labels)
    ax.set_yticks(yticks, labels=[f"{tick:g}" for tick in yticks])
    ax.tick_params(colors=INK_MUTED, labelsize=fontsize, length=2, width=line_width, pad=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_color(INK_SECONDARY)

    span = float(times[-1] - times[0]) or 1.0
    ax.set_xlim(times[0] - 0.06 * span, times[-1] + 0.06 * span)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(xlabel, color=INK_SECONDARY, fontsize=fontsize, labelpad=2)
    ax.set_ylabel(ylabel, color=INK_SECONDARY, fontsize=fontsize, labelpad=2)
    if title:
        ax.set_title(title, color=INK, fontsize=fontsize, loc="left", pad=3)

    if show_legend:
        legend = ax.legend(
            frameon=False,
            loc=legend_loc,
            fontsize=fontsize,
            handlelength=1.4,
            handletextpad=0.4,
            labelspacing=0.2,
            borderaxespad=0.2,
        )
        for text in legend.get_texts():
            text.set_color(INK_SECONDARY)

    _size_plot_area(fig, ax, figsize_cm)
    return fig


def set_font(font: str) -> None:
    """Render the figure in ``font``, saying so if it is not installed.

    Matplotlib falls back to its own DejaVu Sans for a missing font and only
    mentions it in a warning, which is easy to miss when a journal has asked for
    a specific face.
    """
    from matplotlib import font_manager, rcParams

    installed = {f.name for f in font_manager.fontManager.ttflist}
    if font not in installed:
        close = sorted(name for name in installed if font.split()[0].lower() in name.lower())
        print(
            f"note: font {font!r} is not installed; matplotlib will substitute one."
            + (f" Installed and similar: {', '.join(close)}." if close else ""),
            file=sys.stderr,
        )
    # Set both, so the family the font belongs to resolves to it either way.
    rcParams["font.family"] = font
    rcParams["mathtext.fontset"] = "custom"
    for group in ("font.serif", "font.sans-serif"):
        rcParams[group] = [font] + list(rcParams[group])


def _size_plot_area(fig, ax, figsize_cm: tuple[float, float], pad: float = 0.02) -> None:
    """Make the axes itself measure figsize_cm, with the labels around it.

    ``--figsize`` is the plotting area, not the file: tick labels, axis titles
    and the legend live outside it. Sizing the canvas instead would leave a
    5x3 cm figure with under half of that for the data once 9 pt labels are
    placed -- and panels are laid out by their plot area anyway.

    The margins each label needs are measured from a first draw, then the canvas
    is resized to the plot area plus those margins and the axes pinned inside it.
    """
    target_w, target_h = figsize_cm[0] / 2.54, figsize_cm[1] / 2.54
    for _ in range(2):  # tick labels can change width when the axes is resized
        fig.canvas.draw()
        inches = fig.dpi_scale_trans.inverted()
        axes_box = ax.get_window_extent().transformed(inches)
        full_box = fig.get_tightbbox()
        left = max(0.0, axes_box.x0 - full_box.x0) + pad
        right = max(0.0, full_box.x1 - axes_box.x1) + pad
        bottom = max(0.0, axes_box.y0 - full_box.y0) + pad
        top = max(0.0, full_box.y1 - axes_box.y1) + pad

        width, height = target_w + left + right, target_h + bottom + top
        fig.set_layout_engine("none")
        fig.set_size_inches(width, height)
        ax.set_position([left / width, bottom / height, target_w / width, target_h / height])


# --------------------------------------------------------------------------- #
# command line
# --------------------------------------------------------------------------- #

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize and plot a cycloheximide chase from ImageJ measurement tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "raw", type=Path, help="ImageJ measurement table for the protein of interest"
    )
    parser.add_argument(
        "control", type=Path, help="ImageJ measurement table for the PGK1 loading control"
    )
    parser.add_argument("conditions", type=int, help="number of conditions")
    parser.add_argument("-t", "--timepoints", type=int, default=4, help="time points per condition")
    parser.add_argument("-r", "--replicates", type=int, default=1, help="replicates per condition")
    parser.add_argument(
        "-c", "--column", default="Mean", help="measurement column to quantify (e.g. Mean, IntDen)"
    )
    parser.add_argument(
        "-b",
        "--background",
        type=Path,
        metavar="TSV",
        help="two-row background table: row 1 for the raw table, row 2 for the PGK1 control",
    )
    parser.add_argument(
        "--order",
        default="ctr",
        choices=["ctr", "crt", "tcr", "trc", "rct", "rtc"],
        help="row order in the tables: c=condition, t=time point, r=replicate, "
        "leftmost varies slowest",
    )
    parser.add_argument("--names", nargs="+", metavar="NAME", help="condition names, in table order")
    parser.add_argument(
        "--times",
        nargs="+",
        type=float,
        metavar="T",
        help="time of each time point (default: 30 min steps, i.e. 0 30 60 90)",
    )
    parser.add_argument("--time-unit", default="min", help="unit for --times")
    parser.add_argument("--title", default="", help="plot title (default: none)")
    parser.add_argument(
        "--ylabel", default="% substrate remaining", help="y-axis label"
    )
    parser.add_argument(
        "--figsize", nargs=2, type=float, default=[3.0, 3.0], metavar=("W", "H"),
        help="size of the plot area in cm (labels and legend sit outside it)",
    )
    parser.add_argument(
        "--legend-loc", default="none",
        help="add a legend at this position, e.g. 'lower left', 'upper right'",
    )
    parser.add_argument(
        "--yticks", nargs="+", type=float, default=[0.0, 50.0, 100.0], metavar="PCT",
        help="y-axis ticks, in percent",
    )
    parser.add_argument("--fontsize", type=float, default=9.0, help="text size in pt")
    parser.add_argument(
        "--font", default="Liberation Serif", help="font family for all text in the figure"
    )
    parser.add_argument("-o", "--out", type=Path, default=Path("chx_chase"), help="output prefix")
    parser.add_argument("--format", default="png", choices=["png", "pdf", "svg"], help="figure format")
    parser.add_argument("--dpi", type=int, default=300, help="figure resolution")
    parser.add_argument("--show", action="store_true", help="open the figure in a window as well")
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    for name, value in (
        ("conditions", args.conditions),
        ("timepoints", args.timepoints),
        ("replicates", args.replicates),
    ):
        if value < 1:
            raise InputError(f"--{name} must be at least 1, got {value}")
    if args.timepoints < 2:
        raise InputError("need at least 2 time points to plot a chase")
    if args.names is not None and len(args.names) != args.conditions:
        raise InputError(f"--names needs {args.conditions} entries, got {len(args.names)}")
    if args.times is not None and len(args.times) != args.timepoints:
        raise InputError(f"--times needs {args.timepoints} entries, got {len(args.times)}")


def run(args: argparse.Namespace) -> int:
    validate_args(args)

    expected = args.conditions * args.timepoints * args.replicates
    out_prefix = args.out
    if out_prefix.parent != Path(""):
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

    tables = {}
    for role, path in (("raw", args.raw), ("PGK1", args.control)):
        table = read_imagej_table(path)
        if len(table) != expected:
            raise InputError(
                f"{path.name} has {len(table)} rows but "
                f"{args.conditions} conditions x {args.timepoints} time points x "
                f"{args.replicates} replicates = {expected} were expected"
            )
        tables[role] = table

    signal = pick_measurement(tables["raw"], args.column, args.raw)
    control = pick_measurement(tables["PGK1"], args.column, args.control)

    background = (0.0, 0.0)
    if args.background is not None:
        signal, control, background = subtract_background(
            signal, control, args.background, args.column
        )
        print(
            f"background subtracted ({args.column}): raw -{background[0]:g}, "
            f"PGK1 -{background[1]:g}"
        )

    layout = build_layout(args.conditions, args.timepoints, args.replicates, args.order)
    tidy = normalize(layout, signal, control)
    summary = summarize(tidy)

    condition_names = args.names or [f"Condition {i + 1}" for i in range(args.conditions)]
    if args.times is not None:
        times = np.asarray(args.times, dtype=float)
    else:
        times = 30.0 * np.arange(args.timepoints, dtype=float)
    time_labels = [f"{t:g}" for t in times]
    xlabel = f"Time ({args.time_unit})" if args.time_unit else "Time"

    if args.conditions > len(PALETTE):
        print(
            f"note: {args.conditions} conditions exceed the {len(PALETTE)}-color palette; "
            "colors repeat with a changed line style.",
            file=sys.stderr,
        )

    # Readable copies of the numbers behind the plot. Signal and control are
    # background-subtracted; percent_remaining is what gets plotted.
    time_column = f"time_{args.time_unit}" if args.time_unit else "time"
    labeled = tidy.assign(
        condition_name=[condition_names[i] for i in tidy["condition"]],
        **{time_column: [times[i] for i in tidy["timepoint"]]},
        replicate_no=tidy["replicate"] + 1,
        signal_background=background[0],
        control_background=background[1],
    ).rename(columns={"percent": "percent_remaining"})[
        [
            "condition_name", time_column, "replicate_no",
            "signal", "signal_background", "control", "control_background",
            "ratio", "percent_remaining",
        ]
    ]
    labeled_summary = summary.assign(
        condition_name=[condition_names[i] for i in summary["condition"]],
        **{time_column: [times[i] for i in summary["timepoint"]]},
    ).rename(columns={"mean": "mean_percent", "sd": "sd_percent"})[
        ["condition_name", time_column, "n", "mean_percent", "sd_percent"]
    ]

    per_replicate_path = out_prefix.with_name(out_prefix.name + "_per_replicate.csv")
    summary_path = out_prefix.with_name(out_prefix.name + "_summary.csv")
    figure_path = out_prefix.with_name(f"{out_prefix.name}.{args.format}")
    labeled.to_csv(per_replicate_path, index=False, float_format="%.6g")
    labeled_summary.to_csv(summary_path, index=False, float_format="%.6g")

    import matplotlib

    if not args.show:
        matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    figure = plot_chase(
        summary=summary,
        condition_names=condition_names,
        times=times,
        time_labels=time_labels,
        xlabel=xlabel,
        ylabel=args.ylabel,
        title=args.title,
        show_errorbars=args.replicates > 1,
        # A lone unnamed series needs no legend box; the title carries it.
        show_legend=(
            args.legend_loc.lower() != "none"
            and (args.conditions > 1 or args.names is not None)
        ),
        figsize_cm=tuple(args.figsize),
        legend_loc=args.legend_loc,
        yticks=args.yticks,
        fontsize=args.fontsize,
        font=args.font,
    )
    # No bbox_inches="tight" here: it would crop the canvas and the saved figure
    # would no longer be the requested size.
    figure.savefig(figure_path, dpi=args.dpi, facecolor=SURFACE)

    print(f"wrote {figure_path}")
    print(f"wrote {per_replicate_path}")
    print(f"wrote {summary_path}")
    if args.replicates == 1:
        print("note: one replicate, so no standard deviation was calculated.")

    if args.show:
        plt.show()
    plt.close(figure)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        return run(args)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
