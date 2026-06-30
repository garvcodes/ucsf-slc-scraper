"""
Generate publication-quality figures for the SLC meta-analysis.

Renders the 13 manuscript figures (9 main + 4 supplemental) with a custom
typographic and color system: Inter font (loaded from ./fonts), refined slate
palette with accent highlights, title/subtitle/caption hierarchy, and direct
labeling on single-series charts.

Reads:  raw_output_web_scraping.xlsx
Writes: figures/*.png

Dependencies: pandas openpyxl matplotlib seaborn
Fonts:        ./fonts/Inter-*.ttf
"""

from pathlib import Path
import io
import re
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import seaborn as sns

PROJECT_DIR = Path(__file__).parent
INPUT_PATH  = PROJECT_DIR / "raw_output_web_scraping.xlsx"
GENE_LIST   = PROJECT_DIR / "new_slc_genes.csv"   # full set of genes queried (coverage denominator)
FOLDS_PATH  = PROJECT_DIR / "slc_folds_tmscore.csv"  # per-PDB fold calls from superimpose.py (US-align)
OUT_DIR     = PROJECT_DIR / "figures"
FONT_DIR    = PROJECT_DIR / "fonts"

DPI         = 220
YEAR_CUTOFF = None  # e.g. 2025 -> drop entries with release year >= 2025
CAPTION     = "Source: RCSB PDB"

# --- PLOS ONE submission profile ---
# When rendering in PLOS mode the figures are emitted as TIFFs that meet the
# journal's spec: Arial 8-12 pt, 300 dpi, RGB, LZW-compressed, width <= 7.5 in,
# height <= 8.75 in, sequential Fig#.tif naming, and NO embedded title/caption
# (those live in the manuscript text — see figure_captions.md).
OUT_DIR_PLOS = PROJECT_DIR / "figures_plos"
PLOS_DPI     = 300
PLOS_MAX_W   = 7.2    # inches — design clamp (tight-bbox padding pushes the saved size up toward 7.5)
PLOS_MAX_H   = 8.5    # inches
PLOS_MAX_W_PX = 2250  # hard pixel ceilings from the PLOS spec (7.5 x 8.75 in @ 300 dpi)
PLOS_MAX_H_PX = 2625
_PLOS        = False  # toggled by render_all()

# --- Colour system ---
#
# The whole deck is built on the Okabe-Ito palette, the colourblind-safe
# categorical set Nature and other journals recommend (orange/sky/green/blue/
# vermillion/purple + neutral gray). Single-hue *sequential* ramps (blue) are
# used for ordered data, and a single accent is used to highlight the leading
# bar. Restraint here is what reads as "journal" rather than "spreadsheet".

# Neutral text ramp
INK        = "#1A1A1A"   # near-black, for titles
TEXT       = "#3D3D3D"   # body text
MUTED      = "#6B6B6B"   # secondary
FAINT      = "#9A9A9A"   # captions

BG         = "#FFFFFF"
GRID       = "#ECECEC"   # ultra-faint value-axis gridlines
SPINE      = "#BFBFBF"

# Okabe-Ito base hues
OI_ORANGE  = "#E69F00"
OI_SKY     = "#56B4E9"
OI_GREEN   = "#009E73"
OI_YELLOW  = "#F0E442"
OI_BLUE    = "#0072B2"
OI_VERM    = "#D55E00"
OI_PURPLE  = "#CC79A7"
OI_GRAY    = "#999999"

# Bar colours
NEUTRAL    = "#C8D2DC"   # calm blue-gray for de-emphasized bars
PRIMARY    = OI_BLUE     # single-series emphasis
ACCENT     = OI_VERM     # highlight the leading bar

# Method: cryo-EM = orange, X-ray = blue, NMR = green (Okabe-Ito, colourblind-safe).
METHOD_PALETTE = {
    "ELECTRON MICROSCOPY": OI_ORANGE,   # cryo-EM
    "X-RAY DIFFRACTION":   OI_BLUE,     # X-ray
    "SOLUTION NMR":        OI_GREEN,    # NMR
}

ORGANISM_TYPE_PALETTE = {
    "Eukaryotic":  OI_BLUE,
    "Prokaryotic": OI_ORANGE,
    "Virus":       OI_VERM,
    "Synthetic":   OI_GRAY,
}

# Nominal categorical palette (kept <= 8 distinct hues per the journal guidance);
# pooled "Other" buckets get NEUTRAL gray, appended by the caller.
QUAL_PALETTE = [OI_BLUE, OI_ORANGE, OI_GREEN, OI_SKY, OI_VERM, OI_PURPLE, OI_YELLOW, "#6E6E6E"]


def seq_palette(n):
    """Single-hue blue sequential ramp for *ordered* data (resolution bins, or
    size-ranked donut slices). Darkest first, so the largest / sharpest reads
    strongest. A monochrome gradient avoids the rainbow look journals flag."""
    return sns.color_palette("Blues", n_colors=n + 1)[1:][::-1]


def qual_colors(n, with_other=False):
    """First `n` qualitative colours (cycling if needed). When `with_other`
    is set, the final colour is replaced with NEUTRAL gray for a pooled bucket."""
    base = [QUAL_PALETTE[i % len(QUAL_PALETTE)] for i in range(n)]
    if with_other and base:
        base[-1] = NEUTRAL
    return base


# --- Font + style setup ---

def load_fonts():
    """Register local Inter font files with matplotlib."""
    for ttf in FONT_DIR.glob("Inter-*.ttf"):
        fm.fontManager.addfont(str(ttf))


def setup_style(plos=False):
    load_fonts()
    font_family = "Arial" if plos else "Inter"
    dpi = PLOS_DPI if plos else DPI
    plt.rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "figure.facecolor": BG,
        "axes.facecolor": BG,
        "savefig.facecolor": BG,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.3,

        "axes.edgecolor": SPINE,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,

        "font.family": font_family,
        "font.weight": "regular",

        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "axes.titlecolor": INK,
        "axes.titlepad": 14,
        "axes.titlelocation": "left",

        "axes.labelsize": 11,
        "axes.labelcolor": TEXT,
        "axes.labelweight": "regular",
        "axes.labelpad": 10,

        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "xtick.major.pad": 6,
        "ytick.major.pad": 6,

        "legend.frameon": False,
        "legend.fontsize": 10,
        "legend.title_fontsize": 10,
        "legend.labelcolor": TEXT,

        "axes.grid": False,
        "axes.axisbelow": True,
    })


# --- Data loading ---

def classify_kingdom(k):
    """Map the scraper's KINGDOM column (NCBI superkingdom) to the figure labels."""
    if pd.isna(k) or not str(k).strip():
        return None
    if k in ("Bacteria", "Archaea"):
        return "Prokaryotic"
    if k == "Eukaryota":
        return "Eukaryotic"
    if k == "Viruses":
        return "Virus"
    return None


def senior_author(s):
    if pd.isna(s):
        return None
    parts = [p.strip() for p in str(s).split(";") if p.strip()]
    return parts[-1] if parts else None


def norm_oligomer(s):
    """Collapse case-variant oligomeric-state labels (e.g. 'Dimeric' -> 'dimeric')."""
    if pd.isna(s) or not str(s).strip():
        return None
    return str(s).strip().lower()


def slc_family(gene):
    """Map a gene symbol to its family token.

    SLC / SLCO symbols collapse to their family number (SLC6A1 -> 'SLC6',
    SLCO1B1 -> 'SLCO1'); SLC-related symbols that don't use the SLC prefix
    (XPR1, NPC1, SV2A, CTNS, MFSD2A …) keep their own letter+number stem so
    they are not all dumped into one bucket.
    """
    g = str(gene).strip()
    m = re.match(r"(SLCO\d+|SLC\d+)", g)
    if m:
        return m.group(1)
    m = re.match(r"([A-Z]+\d*)", g)
    return m.group(1) if m else g


def load_data(path, year_cutoff=None):
    df = pd.read_excel(path)
    df["release_year"] = pd.to_datetime(df["RELEASE DATE"], errors="coerce").dt.year
    df["RESOLUTION"]   = pd.to_numeric(df["RESOLUTION"], errors="coerce")
    bins = list(range(1, 10))
    labels = [f"{b}-{b+1}" for b in bins[:-1]]
    df["resolution_bin"] = pd.cut(df["RESOLUTION"], bins=bins, labels=labels, right=False)
    # Prefer the materialized PROK_EUK column from the scraper; fall back to deriving
    # it from KINGDOM for older output files that predate that column.
    if "PROK_EUK" in df.columns:
        df["organism_type"] = df["PROK_EUK"].replace("", None)
    else:
        df["organism_type"] = df["KINGDOM"].apply(classify_kingdom)
    df["senior_author"]  = df["STRUCTURE AUTHOR"].apply(senior_author)
    if year_cutoff is not None:
        df = df[df["release_year"] < year_cutoff]
    return df


# --- Plot chrome ---

def panel_tag(ax, letter):
    """Bold panel label (A, B, ...) at the top-left of a sub-panel, model-paper style."""
    ax.text(-0.10, 1.08, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", color=INK, va="top", ha="left")


def panel_title(ax, text):
    ax.set_title(text, fontsize=9, loc="left", weight="bold", color=INK, pad=6)


def panel_legend(ax, **kw):
    opts = dict(fontsize=7.5, frameon=False, loc="upper left", handletextpad=0.4,
                labelspacing=0.25)
    opts.update(kw)
    ax.legend(**opts)


def thin_xticks(ax, every=3, fontsize=7):
    """Keep every Nth x-tick label (declutters dense year axes in small panels)."""
    for i, lab in enumerate(ax.get_xticklabels()):
        lab.set_visible(i % every == 0)
        lab.set_fontsize(fontsize)
        lab.set_rotation(45)
        lab.set_ha("right")


def with_chrome(fig, ax, title, subtitle=None, caption=CAPTION):
    """Apply the title / subtitle / caption layout. Replaces ax.set_title.

    In PLOS mode the title, subtitle and source caption are NOT drawn — PLOS ONE
    requires them in the manuscript text, not the image file — so the figure is
    laid out edge-to-edge instead. The titles still flow into figure_captions.md.
    """
    if _PLOS:
        fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.99])
        return
    # Title at top-left, in figure coords
    fig.text(0.04, 0.97, title, fontsize=16, weight="bold",
             color=INK, ha="left", va="top")
    if subtitle:
        fig.text(0.04, 0.925, subtitle, fontsize=11,
                 color=MUTED, ha="left", va="top")
    if caption:
        fig.text(0.985, 0.02, caption, fontsize=9,
                 color=FAINT, ha="right", va="bottom")
    # Reserve top for title/subtitle, bottom for caption
    top = 0.86 if subtitle else 0.89
    fig.tight_layout(rect=[0.02, 0.05, 0.98, top])


def highlight_top(values, base=NEUTRAL, accent=ACCENT, n=1):
    """Per-bar colors with the top n bars highlighted."""
    ranked = pd.Series(values).rank(method="min", ascending=False)
    return [accent if r <= n else base for r in ranked]


def order_methods(columns):
    return [m for m in METHOD_PALETTE if m in columns] + \
           [m for m in columns if m not in METHOD_PALETTE]


def add_value_labels(ax, fontsize=9, color=TEXT, padding=4, weight="medium"):
    for container in ax.containers:
        ax.bar_label(container, padding=padding, fontsize=fontsize,
                     color=color, weight=weight)


# --- Clean primitives (horizontal bars / stacked area), journal house style ---

def _strip_axes(ax, keep="left"):
    """Minimal chrome: drop all spines/ticks except the one anchoring the labels.
    Horizontal-bar charts keep the left (category) edge; everything else goes so
    the bars and their end-labels carry the figure."""
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(sp == keep)
    if keep == "left":
        ax.tick_params(axis="x", length=0)
        ax.set_xticks([])
        ax.tick_params(axis="y", length=0)


def hbar(ax, labels, values, *, highlight_n=1, italic=False, fontsize=8,
         value_fontsize=8, xmax=None, bar_height=0.72):
    """Single-series horizontal bar chart, largest at the top, leading bar accented,
    each value labeled at the bar end. The readable alternative to rotated x-labels.

    `xmax` forces a shared x-scale (so bars in side-by-side panels stay comparable)."""
    values = np.asarray(values, dtype=float)
    y = np.arange(len(labels))[::-1]
    colors = highlight_top(values, n=highlight_n)
    ax.barh(y, values, color=colors, height=bar_height, edgecolor="none", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=fontsize,
                       fontstyle="italic" if italic else "normal")
    top = xmax if xmax is not None else values.max() * 1.13
    pad = top * 0.012
    for yi, v in zip(y, values):
        ax.text(v + pad, yi, f"{int(round(v))}", va="center", ha="left",
                fontsize=value_fontsize, color=TEXT)
    ax.set_xlim(0, top)
    ax.set_ylim(-0.7, len(labels) - 0.3)
    _strip_axes(ax)


def hbar_stacked(ax, labels, pivot, palette, *, fontsize=8, value_fontsize=8):
    """Stacked horizontal bars (e.g. per-transporter counts split by method), largest
    total at the top, total labeled at the bar end."""
    y = np.arange(len(labels))[::-1]
    left = np.zeros(len(labels))
    for col in pivot.columns:
        vals = pivot[col].to_numpy(dtype=float)
        ax.barh(y, vals, left=left, color=palette.get(col, NEUTRAL), height=0.72,
                edgecolor="white", linewidth=0.6, label=col, zorder=3)
        left += vals
    totals = pivot.sum(axis=1).to_numpy(dtype=float)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=fontsize)
    pad = totals.max() * 0.015
    for yi, t in zip(y, totals):
        ax.text(t + pad, yi, f"{int(round(t))}", va="center", ha="left",
                fontsize=value_fontsize, color=TEXT)
    ax.set_xlim(0, totals.max() * 1.13)
    ax.set_ylim(-0.7, len(labels) - 0.3)
    _strip_axes(ax)


def stacked_area(ax, years, pivot, colors):
    """Smooth stacked-area timeline — calmer than stacked bars with white gridlines.
    `pivot` is years x categories; columns are drawn bottom-to-top in order."""
    ys = [pivot[c].to_numpy(dtype=float) for c in pivot.columns]
    ax.stackplot(years, ys, colors=colors, labels=list(pivot.columns),
                 edgecolor="white", linewidth=0.4)
    ax.set_xlim(min(years), max(years))
    ax.set_ylim(bottom=0)
    ax.margins(x=0)
    # Sparse, unrotated year ticks at 5-year marks read far cleaner than every bar.
    lo = (int(min(years)) + 4) // 5 * 5
    ticks = list(range(lo, int(max(years)) + 1, 5))
    ax.set_xticks(ticks)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=0)


def ranked_seq_colors(n, with_other=False):
    """Blue sequential ramp for `n` size-ranked slices (darkest = largest); a pooled
    'Other' slice, if present, is the trailing NEUTRAL gray."""
    k = n - (1 if with_other else 0)
    cols = list(seq_palette(k))
    if with_other:
        cols.append(NEUTRAL)
    return cols


def _on_color(c):
    """Black or white, whichever stays legible on background colour `c`."""
    r, g, b = mcolors.to_rgb(c)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return INK if luminance > 0.6 else "white"


def donut(ax, counts, colors, label_threshold=0.03):
    """Draw a donut chart of `counts` (a Series) onto `ax`, returning the wedges.

    Slices at or above `label_threshold` of the total are labeled in place with
    their percentage; everything is also spelled out in the legend the caller
    attaches. The hole is left open for a center total.
    """
    total = int(counts.sum())
    wedges, _ = ax.pie(
        counts.values, colors=colors, startangle=90, counterclock=False,
        wedgeprops=dict(width=0.42, edgecolor=BG, linewidth=2.2),
    )
    for wedge, value, c in zip(wedges, counts.values, colors):
        frac = value / total if total else 0
        if frac < label_threshold:
            continue
        ang = (wedge.theta2 + wedge.theta1) / 2
        x = 0.79 * np.cos(np.deg2rad(ang))
        y = 0.79 * np.sin(np.deg2rad(ang))
        ax.text(x, y, f"{frac * 100:.0f}%", ha="center", va="center",
                fontsize=9.5, weight="semibold", color=_on_color(c))
    ax.text(0, 0, f"{total}", ha="center", va="center",
            fontsize=22, weight="bold", color=INK)
    ax.text(0, -0.16, "structures", ha="center", va="center",
            fontsize=9.5, color=MUTED)
    ax.set(aspect="equal")
    return wedges


def add_donut_inset(host_ax, counts, colors, caption=None,
                    loc=(0.0, 0.56, 0.30, 0.42)):
    """Embed a compact donut of `counts` in a corner of `host_ax`.

    Shares the host chart's category colors, so the host's own legend doubles as
    the inset's key — the inset just shows the aggregate composition at a glance.
    Kept deliberately spare: a wide ring with the major slice shares labeled and a
    single centred total. The old version also stacked a text caption inside the
    hole, which collided with the total and the slice labels at this size; that is
    gone (the panel title already says what the inset summarizes). `caption` is
    accepted for backward compatibility but no longer drawn.
    """
    ax = host_ax.inset_axes(loc)
    total = int(counts.sum())
    # Deliberately spare: a colour-coded proportion ring (its key is the host panel's
    # own legend) with just the centred total. Per-slice percentages live in the
    # dedicated Fig4/Fig5 donuts and the manuscript text — crammed into a ring this
    # small they clipped against the inset boundary and read as noise.
    ax.pie(counts.values, colors=colors, startangle=90, counterclock=False,
           wedgeprops=dict(width=0.34, edgecolor=BG, linewidth=1.3))
    ax.text(0, 0, f"{total}", ha="center", va="center",
            fontsize=11, weight="bold", color=INK, clip_on=False)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set(aspect="equal")
    return ax


def donut_legend(ax, counts, colors):
    """Right-side legend: 'Label — n (xx%)' for each slice, largest first."""
    total = int(counts.sum())
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=9,
                   markerfacecolor=c, markeredgecolor="none",
                   label=f"{lab}   {v}  ({v / total * 100:.0f}%)")
        for lab, v, c in zip(counts.index, counts.values, colors)
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, fontsize=10, handletextpad=0.6, labelspacing=0.7)


def save_fig(fig, name):
    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / f"{name}.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"saved {out}")


def save_fig_plos(fig, label):
    """Save `fig` as a PLOS ONE-compliant TIFF named `label` (e.g. 'Fig1').

    Conforms to the journal's raster spec: <= 7.5 x 8.75 in, 300 dpi, RGB (alpha
    flattened onto white), LZW-compressed TIFF, hard pixel ceilings enforced.
    Title/caption are suppressed upstream (they belong in the manuscript text).
    """
    OUT_DIR_PLOS.mkdir(exist_ok=True)
    w, h = fig.get_size_inches()
    scale = min(PLOS_MAX_W / w, PLOS_MAX_H / h)   # fit the design box, preserve aspect
    fig.set_size_inches(w * scale, h * scale)

    buf = io.BytesIO()
    fig.savefig(buf, format="tiff", dpi=PLOS_DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    buf.seek(0)

    im = Image.open(buf)
    # Flatten any alpha onto white -> RGB (PLOS rejects RGBA/CMYK).
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.split()[-1])
        im = bg
    else:
        im = im.convert("RGB")
    # Enforce the hard pixel ceilings (tight bbox can nudge past 7.5 in).
    if im.width > PLOS_MAX_W_PX or im.height > PLOS_MAX_H_PX:
        f = min(PLOS_MAX_W_PX / im.width, PLOS_MAX_H_PX / im.height)
        im = im.resize((round(im.width * f), round(im.height * f)), Image.LANCZOS)

    out = OUT_DIR_PLOS / f"{label}.tif"
    im.save(out, format="TIFF", compression="tiff_lzw", dpi=(PLOS_DPI, PLOS_DPI))
    print(f"saved {out}  ({im.width} x {im.height}px, {im.width / PLOS_DPI:.2f} x {im.height / PLOS_DPI:.2f} in)")


# --- Main figures ---

def fig_structures_by_year(df):
    """Total unique PDB structures deposited per release year.

    This is the raw deposition-activity view (e.g. the cryo-EM-driven surge in
    2023-2025). Distinct from cumulative SLC *coverage* below: many of these recent
    structures are repeat depositions of already-characterized transporters.
    """
    sub = df.drop_duplicates("PDB ID").dropna(subset=["release_year"])
    counts = sub["release_year"].astype(int).value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(counts.index.astype(str), counts.values,
           color=highlight_top(counts.values), edgecolor="white",
           linewidth=0.6, width=0.72)
    ax.set_xlabel("Release Year")
    ax.set_ylabel("Structures Deposited")
    add_value_labels(ax)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    with_chrome(fig, ax,
        "SLC Structures Deposited per Year",
        f"Unique PDB depositions by release year  ·  {int(counts.sum())} total structures")
    return fig


def fig_cumulative_slc_coverage(df):
    """Cumulative count of distinct SLCs with at least one deposited structure, by year.

    Plotted cumulatively (not as per-year first-appearances) so the curve reflects
    growing structural coverage and does not read as "zero new structures" in years
    where every newly-solved SLC had already been solved before — the misreading the
    per-year-debut version invited.
    """
    sub = df.dropna(subset=["release_year"])
    first_year = sub.groupby("GENE")["release_year"].min().astype(int)
    debuts = first_year.value_counts().sort_index()
    years = range(int(debuts.index.min()), int(debuts.index.max()) + 1)
    cumulative = debuts.reindex(years, fill_value=0).cumsum()

    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = [str(y) for y in cumulative.index]
    ax.fill_between(x, cumulative.values, color=PRIMARY, alpha=0.12)
    ax.plot(x, cumulative.values, color=PRIMARY, linewidth=2.4,
            marker="o", markersize=4, markerfacecolor=PRIMARY, markeredgecolor="white")
    # Label the final point with the running total.
    ax.annotate(f"{int(cumulative.iloc[-1])}", (x[-1], cumulative.iloc[-1]),
                textcoords="offset points", xytext=(0, 8), ha="center",
                fontsize=10, weight="bold", color=INK)
    ax.set_xlabel("Release Year")
    ax.set_ylabel("Distinct SLCs Solved (cumulative)")
    ax.set_ylim(bottom=0)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    # Denominator = full gene list queried, not just genes that happen to have a structure.
    total_genes = df["GENE"].nunique()
    if GENE_LIST.exists():
        total_genes = pd.read_csv(GENE_LIST).iloc[:, 0].dropna().str.strip().nunique()
    with_chrome(fig, ax,
        "Cumulative Structural Coverage of SLC Transporters",
        f"Distinct SLCs with ≥ 1 deposited structure  ·  {int(cumulative.iloc[-1])} of {total_genes} genes queried")
    return fig


def fig_new_slcs_by_year(df):
    """Number of SLCs solved for the first time in each release year (debuts only).

    Per-year first-appearance counts: each gene is counted once, in the year its
    earliest structure was released. Distinct from fig_structures_by_year (which
    counts every unique PDB, including repeat depositions) and from the cumulative
    coverage curve (the running total of these same debuts). Recent years read low
    here by design — most genes debuted earlier — so pair it with the deposition
    and cumulative views rather than reading it alone.
    """
    sub = df.dropna(subset=["release_year"])
    first_year = sub.groupby("GENE")["release_year"].min().astype(int)
    debuts = first_year.value_counts().sort_index()
    years = range(int(debuts.index.min()), int(debuts.index.max()) + 1)
    debuts = debuts.reindex(years, fill_value=0)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(debuts.index.astype(str), debuts.values,
           color=highlight_top(debuts.values), edgecolor="white",
           linewidth=0.6, width=0.72)
    ax.set_xlabel("Release Year")
    ax.set_ylabel("Newly-Solved SLCs")
    add_value_labels(ax)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    with_chrome(fig, ax,
        "New SLCs Solved per Year",
        f"First-time structural debuts by release year  ·  {int(debuts.sum())} distinct SLCs")
    return fig


def fig_slc_debut_timeline(df):
    """Vertical timeline naming every SLC in the year it was first solved.

    Same debut events as fig_new_slcs_by_year, but spelled out: each gene appears
    once, beside the year of its earliest structure, colored by the organism type
    of that first structure. Years with no new SLC are still drawn so the temporal
    spacing (and the recent acceleration) reads honestly.
    """
    sub = df.dropna(subset=["release_year"]).copy()
    sub["release_year"] = sub["release_year"].astype(int)
    # The gene's earliest structure: debut year + organism type at debut.
    debut = (sub.loc[sub.groupby("GENE")["release_year"].idxmin()]
                [["GENE", "release_year", "organism_type"]])
    by_year = {y: list(g.sort_values("GENE")[["GENE", "organism_type"]]
                       .itertuples(index=False, name=None))
               for y, g in debut.groupby("release_year")}
    ymin, ymax = int(debut["release_year"].min()), int(debut["release_year"].max())

    # Names sit in a fixed aligned grid (clean columns) rather than wrapping at
    # variable widths. Each year gets exactly the vertical room its names need, so
    # busy recent years stop colliding and sparse early years stay compact.
    COLS = 8
    spine_x, name_x0, name_x1 = 0.075, 0.165, 0.99
    colw = (name_x1 - name_x0) / COLS
    ROW_H, GAP, EMPTY = 1.0, 0.62, 0.5

    layout, yc = [], 0.0
    for y in range(ymin, ymax + 1):
        genes = by_year.get(y, [])
        if genes:
            rows = int(np.ceil(len(genes) / COLS))
            layout.append((y, yc, rows, genes))
            yc += (rows - 1) * ROW_H + 1.0 + GAP
        else:
            layout.append((y, yc, 0, []))
            yc += EMPTY
    total = yc

    fig, ax = plt.subplots(figsize=(11, min(13.0, total * 0.225 + 1.4)))
    ax.set_xlim(0, 1)
    ax.set_ylim(total - 0.4, -1.4)   # earliest at top
    ax.axis("off")

    y_first = layout[0][1]
    y_last = layout[-1][1]
    ax.plot([spine_x, spine_x], [y_first, y_last], color=SPINE, linewidth=1.5, zorder=1)

    band = False
    for (y, y0, rows, genes) in layout:
        has = bool(genes)
        if has:
            if band:   # subtle alternating band groups each populated year
                ax.axhspan(y0 - 0.5, y0 + (rows - 1) * ROW_H + 0.5,
                           xmin=0.02, xmax=0.996, color="#F5F7F9", zorder=0)
            band = not band
        ax.scatter([spine_x], [y0], s=46 if has else 13,
                   color=PRIMARY if has else SPINE, zorder=3,
                   edgecolor=BG, linewidth=1.1 if has else 0)
        # Label every populated year; for the compressed runs of empty years only
        # label 5-year marks so the gap years (e.g. 1995-2003) don't stack up.
        if has or y % 5 == 0:
            ax.text(spine_x - 0.016, y0, str(y), ha="right", va="center",
                    fontsize=9.5, color=INK if has else FAINT,
                    weight="bold" if has else "regular")
        # Count of SLCs first solved that year, in the gutter next to the names
        # (this folds in the old "new SLCs per year" bar chart).
        if has:
            ax.text(spine_x + 0.030, y0, f"{len(genes)}", ha="left", va="center",
                    fontsize=8.5, color=ACCENT, weight="bold")
        for i, (name, otype) in enumerate(genes):
            r, c = divmod(i, COLS)
            ax.text(name_x0 + c * colw, y0 + r * ROW_H, name, ha="left", va="center",
                    fontsize=8, color=ORGANISM_TYPE_PALETTE.get(otype, MUTED),
                    weight="medium")

    present = [t for t in ORGANISM_TYPE_PALETTE if t in set(debut["organism_type"].dropna())]
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=7,
                          markerfacecolor=ORGANISM_TYPE_PALETTE[t], markeredgecolor="none",
                          label=t) for t in present]
    if handles:
        ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.99, 0.005),
                  frameon=False, fontsize=9, handletextpad=0.4, labelspacing=0.5,
                  ncol=len(handles))

    with_chrome(fig, ax,
        "Timeline of First-Solved SLC Transporters",
        f"Each SLC at its structural debut, colored by organism type; bold number = SLCs "
        f"first solved that year  ·  {len(debut)} distinct SLCs")
    return fig


def fig_top_slcs(df, top_n=10):
    counts = df.drop_duplicates("PDB ID")["GENE"].value_counts().head(top_n)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(counts.index, counts.values,
           color=highlight_top(counts.values), edgecolor="white",
           linewidth=0.6, width=0.7)
    ax.set_xlabel("SLC Transporter")
    ax.set_ylabel("Structures")
    add_value_labels(ax)
    plt.setp(ax.get_xticklabels(), rotation=0)
    with_chrome(fig, ax,
        f"Top {top_n} SLC Transporters by Number of Structures",
        f"Unique PDB depositions  ·  {df['PDB ID'].nunique()} total structures")
    return fig


def fig_top_authors(df, top_n=10):
    counts = (df.dropna(subset=["senior_author"])
                .drop_duplicates("PDB ID")["senior_author"]
                .value_counts().head(top_n))
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(counts.index, counts.values,
           color=highlight_top(counts.values), edgecolor="white",
           linewidth=0.6, width=0.7)
    ax.set_xlabel("Senior Author")
    ax.set_ylabel("Structures")
    add_value_labels(ax)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")
    with_chrome(fig, ax,
        f"Top {top_n} Senior Authors by Unique SLC Structures",
        "Last author on each PDB deposition")
    return fig


def fig_top_organisms(df, top_n=10):
    counts = df.drop_duplicates("PDB ID")["SOURCE ORGANISM"].value_counts().head(top_n)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(counts.index, counts.values,
           color=highlight_top(counts.values), edgecolor="white",
           linewidth=0.6, width=0.7)
    ax.set_xlabel("Source Organism")
    ax.set_ylabel("Structures")
    add_value_labels(ax)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", fontstyle="italic")
    with_chrome(fig, ax,
        f"Top {top_n} Source Organisms",
        "Homo sapiens dominates by an order of magnitude")
    return fig


def fig_organism_type_by_year(df):
    sub = df.drop_duplicates("PDB ID").dropna(subset=["release_year", "organism_type"]).copy()
    sub["release_year"] = sub["release_year"].astype(int)
    pivot = (sub.groupby(["release_year", "organism_type"]).size()
                .unstack(fill_value=0)
                .reindex(columns=list(ORGANISM_TYPE_PALETTE.keys()), fill_value=0))
    pivot = pivot.loc[:, (pivot != 0).any(axis=0)]   # drop categories with no structures
    fig, ax = plt.subplots(figsize=(11, 5.5))
    pivot.plot(kind="bar", stacked=True, ax=ax,
               color=[ORGANISM_TYPE_PALETTE[c] for c in pivot.columns],
               edgecolor="white", linewidth=0.5, width=0.78)
    ax.set_xlabel("Release Year")
    ax.set_ylabel("Structures")
    ax.legend(title="", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    # Inset donut: overall prokaryotic/eukaryotic composition across all years.
    totals = pivot.sum(axis=0)
    add_donut_inset(ax, totals, [ORGANISM_TYPE_PALETTE[c] for c in totals.index])
    with_chrome(fig, ax,
        "Prokaryotic vs Eukaryotic SLC Structures by Year",
        f"Eukaryotic structures dominate  ·  {int(pivot.values.sum())} structures shown")
    return fig


def fig_method_by_year(df):
    sub = df.drop_duplicates("PDB ID").dropna(subset=["release_year"]).copy()
    sub["release_year"] = sub["release_year"].astype(int)
    pivot = sub.groupby(["release_year", "METHOD"]).size().unstack(fill_value=0)
    pivot = pivot[order_methods(pivot.columns)]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    pivot.plot(kind="bar", stacked=True, ax=ax,
               color=[METHOD_PALETTE.get(c, "#777") for c in pivot.columns],
               edgecolor="white", linewidth=0.5, width=0.78)
    ax.set_xlabel("Release Year")
    ax.set_ylabel("Structures")
    ax.legend(title="", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    # Inset donut: overall method composition across all years.
    totals = pivot.sum(axis=0)
    add_donut_inset(ax, totals, [METHOD_PALETTE.get(c, NEUTRAL) for c in totals.index])
    with_chrome(fig, ax,
        "Structure Determination Methodology by Year",
        "Cryo-EM took over after 2019; X-ray dominated the preceding era")
    return fig


def fig_top_slcs_by_method(df, top_n=10):
    unique = df.drop_duplicates("PDB ID")
    top_genes = unique["GENE"].value_counts().head(top_n).index
    sub = unique[unique["GENE"].isin(top_genes)]
    pivot = sub.groupby(["GENE", "METHOD"]).size().unstack(fill_value=0).loc[top_genes]
    pivot = pivot[order_methods(pivot.columns)]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    pivot.plot(kind="bar", stacked=True, ax=ax,
               color=[METHOD_PALETTE.get(c, "#777") for c in pivot.columns],
               edgecolor="white", linewidth=0.5, width=0.72)
    ax.set_xlabel("SLC Transporter")
    ax.set_ylabel("Structures")
    ax.legend(title="", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=0)
    with_chrome(fig, ax,
        f"Top {top_n} SLC Transporters by Method",
        "Method composition for the most-deposited transporters")
    return fig


def fig_resolution_by_method(df):
    sub = df.drop_duplicates("PDB ID").dropna(subset=["resolution_bin"])
    pivot = sub.groupby(["resolution_bin", "METHOD"], observed=True).size().unstack(fill_value=0)
    pivot = pivot[order_methods(pivot.columns)]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    pivot.plot(kind="bar", stacked=True, ax=ax,
               color=[METHOD_PALETTE.get(c, "#777") for c in pivot.columns],
               edgecolor="white", linewidth=0.5, width=0.7)
    ax.set_xlabel("Resolution (Å)")
    ax.set_ylabel("Structures")
    ax.legend(title="", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=0)
    with_chrome(fig, ax,
        "Structure Resolution by Method",
        "X-ray claims the high-resolution tail; EM packs the 3–5 Å band")
    return fig


def fig_resolution_by_year(df):
    sub = df.drop_duplicates("PDB ID").dropna(subset=["release_year"]).copy()
    sub["release_year"] = sub["release_year"].astype(int)
    sub["resolution_bin"] = sub["resolution_bin"].cat.add_categories(["NA"]).fillna("NA")
    bin_order = [f"{b}-{b+1}" for b in range(1, 9)] + ["NA"]
    pivot = (sub.groupby(["release_year", "resolution_bin"], observed=True).size()
                .unstack(fill_value=0).reindex(columns=bin_order, fill_value=0))
    pivot = pivot.loc[:, (pivot != 0).any(axis=0)]
    colors = seq_palette(len(pivot.columns))
    fig, ax = plt.subplots(figsize=(11, 5.5))
    pivot.plot(kind="bar", stacked=True, ax=ax,
               color=colors, edgecolor="white", linewidth=0.5, width=0.78)
    ax.set_xlabel("Release Year")
    ax.set_ylabel("Structures")
    ax.legend(title="Resolution (Å)", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    with_chrome(fig, ax,
        "Structure Resolution Range by Year",
        "Sub-3 Å structures grew with the cryo-EM wave")
    return fig


# --- Supplemental figures (full-list cuts) ---

def fig_all_slcs(df):
    """Full per-gene ranking as two side-by-side horizontal-bar columns (ranks 1-n/2
    and the rest), so all ~85 transporter names stay upright and readable instead of
    being crushed into rotated x-tick labels."""
    counts = df.drop_duplicates("PDB ID")["GENE"].value_counts()
    genes, vals = list(counts.index), counts.values
    n = len(genes)
    half = (n + 1) // 2
    xmax = vals.max() * 1.16   # shared scale so the two columns stay comparable

    rows = half
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, max(7.0, rows * 0.205)))
    hbar(axL, genes[:half], vals[:half], highlight_n=1,
         fontsize=7, value_fontsize=7, xmax=xmax, bar_height=0.66)
    hbar(axR, genes[half:], vals[half:], highlight_n=0,
         fontsize=7, value_fontsize=7, xmax=xmax, bar_height=0.66)
    # keep both columns the same row pitch even though the 2nd has one fewer gene
    axR.set_ylim(axL.get_ylim())
    axL.set_xlabel("Structures", fontsize=9)
    axR.set_xlabel("Structures", fontsize=9)
    with_chrome(fig, axL,
        "All SLC Transporters by Number of Structures",
        f"{n} unique genes  ·  ranked by deposited structures")
    return fig


# --- Pie / donut figures (categorical composition) ---

def fig_oligomeric_pie(df):
    """Share of structures by oligomeric (assembly) state."""
    states = df.drop_duplicates("PDB ID")["OLIGOMERIC STATE"].apply(norm_oligomer)
    counts = states.value_counts()
    # Collapse the rare long tail into 'other' so the donut stays legible.
    keep = counts[counts >= 0.02 * counts.sum()]
    other = counts[counts < 0.02 * counts.sum()].sum()
    counts = keep.copy()
    if other:
        counts["other"] = other
    colors = ranked_seq_colors(len(counts), with_other=bool(other))
    fig, ax = plt.subplots(figsize=(9, 5.6))
    donut(ax, counts, colors)
    donut_legend(ax, counts, colors)
    with_chrome(fig, ax,
        "Structures by Oligomeric State",
        f"Monomeric and dimeric assemblies predominate  ·  {int(counts.sum())} structures")
    return fig


def fig_family_pie(df, top_n=10):
    """Share of structures by SLC family, top families named and the rest pooled."""
    fam = df.drop_duplicates("PDB ID")["GENE"].apply(slc_family)
    counts = fam.value_counts()
    top = counts.head(top_n)
    other = counts.iloc[top_n:].sum()
    if other:
        top = pd.concat([top, pd.Series({f"Other ({len(counts) - top_n} families)": other})])
    colors = ranked_seq_colors(len(top), with_other=bool(other))
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    donut(ax, top, colors)
    donut_legend(ax, top, colors)
    with_chrome(fig, ax,
        f"Structures by SLC Family (top {top_n})",
        f"{counts.size} families represented  ·  {int(counts.sum())} structures")
    return fig


# --- Structural fold classification (US-align / TM-score) ---

# Human-readable names for the fold reference codes used in superimpose.py.
FOLD_LABELS = {
    "MFS":  "MFS",
    "LeuT": "LeuT / APC",
    "Glt":  "Glutamate transporter",
    "CNT":  "CNT",
    "MitC": "Mitochondrial carrier",
    "UraA": "UraA / SLC4-26",
    "NhaA": "NhaA",
    "NCX":  "NCX",
    "AmtB": "AmtB",
    "YiiP": "YiiP / CDF",
    "DMT":  "DMT",
}


def load_folds():
    """Per-PDB fold calls from superimpose.py, or None if not yet generated."""
    if not FOLDS_PATH.exists():
        return None
    f = pd.read_csv(FOLDS_PATH)
    f["PDB ID"] = f["PDB ID"].astype(str).str.upper()
    return f


def _confident_folds(df):
    """Merge fold calls onto unique PDBs and return the confident subset + the
    ordered (fold, label, colour) triples shared by the fold figures."""
    folds = load_folds()
    u = df.drop_duplicates("PDB ID").copy()
    u["PDB ID"] = u["PDB ID"].astype(str).str.upper()
    m = u.merge(folds, on="PDB ID", how="inner")
    conf = m[m["FOLD_CONFIDENT"] == True].copy()              # noqa: E712 (pandas mask)
    counts = conf["FOLD"].value_counts()
    order = list(counts.index)
    colors = qual_colors(len(order))
    return conf, counts, order, colors


def fig_fold_distribution(df):
    """Structural-fold distribution as a horizontal bar chart and a matching donut,
    sharing one per-fold colour scheme: the bars give absolute counts and name each
    fold, the donut gives the proportional split and overall total. Folds are assigned
    by US-align structural alignment against canonical references (TM-score >= 0.5)."""
    conf, counts, order, colors = _confident_folds(df)
    labels = [FOLD_LABELS.get(f, f) for f in order]

    fig, (axBar, axPie) = plt.subplots(1, 2, figsize=(7.4, 3.9),
                                       gridspec_kw={"width_ratios": [1.05, 0.95]})

    # Left: horizontal bars, one colour per fold (the donut shares these colours,
    # so the bar labels double as the donut's key — no separate legend needed).
    y = np.arange(len(order))[::-1]
    vmax = counts.values.max()
    axBar.barh(y, counts.values, color=colors, height=0.74, edgecolor="none", zorder=3)
    axBar.set_yticks(y)
    axBar.set_yticklabels(labels, fontsize=8.5)
    for yi, v in zip(y, counts.values):
        axBar.text(v + vmax * 0.015, yi, str(int(v)), va="center", ha="left",
                   fontsize=8, color=TEXT)
    axBar.set_xlim(0, vmax * 1.15)
    axBar.set_ylim(-0.7, len(order) - 0.3)
    _strip_axes(axBar)
    axBar.set_xlabel("Structures", fontsize=8.5)
    panel_title(axBar, "Structures per fold")

    # Right: donut of the same counts/colours (in-slice percentages + centre total).
    counts_lab = pd.Series(counts.values, index=labels)
    donut(axPie, counts_lab, colors)
    panel_title(axPie, "Proportional split")

    fig.tight_layout(w_pad=2.0)
    return fig


def fig_fold_tmscore(df):
    """Per-fold TM-score distribution (match quality) against the 0.5 same-fold line —
    the evidence that confidently assigned structures are genuine fold members."""
    conf, counts, order, colors = _confident_folds(df)
    labels = [FOLD_LABELS.get(f, f) for f in order]

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    y = np.arange(len(order))[::-1]
    rng = np.random.RandomState(0)
    for yi, f, c in zip(y, order, colors):
        vals = conf.loc[conf["FOLD"] == f, "TM_SCORE_TO_FOLD_REF"].to_numpy(dtype=float)
        jit = yi + (rng.rand(len(vals)) - 0.5) * 0.55
        ax.scatter(vals, jit, s=11, color=c, alpha=0.5, edgecolor="none", zorder=3)
        ax.scatter([np.median(vals)], [yi], s=55, color=c, zorder=4,
                   edgecolor="white", linewidth=1.1)
    ax.axvline(0.5, color=MUTED, linestyle="--", linewidth=1.0, zorder=2)
    ax.text(0.5, len(order) - 0.3, "same-fold threshold", fontsize=7.5,
            color=MUTED, ha="center", va="top")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(0.4, 1.02)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel("TM-score to fold reference  (large dot = per-fold median)", fontsize=9)
    ax.tick_params(length=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    with_chrome(fig, ax,
        "Structural-fold match quality by TM-score",
        "Each structure's TM-score to its assigned fold reference")
    return fig


# --- Composite multi-panel figures (model-paper style) ---

def fig_structures_and_timeline(df):
    """(A) structures per transporter by method (horizontal); (B) deposition + coverage."""
    u = df.drop_duplicates("PDB ID")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.9))

    # Panel A: top transporters as horizontal stacked bars (method), largest on top.
    top_genes = u["GENE"].value_counts().head(12).index
    sub = u[u["GENE"].isin(top_genes)]
    pivot = (sub.groupby(["GENE", "METHOD"]).size().unstack(fill_value=0)
                .reindex(top_genes))
    pivot = pivot[order_methods(pivot.columns)]
    hbar_stacked(axA, list(top_genes), pivot, METHOD_PALETTE,
                 fontsize=8, value_fontsize=7.5)
    axA.set_xlabel("Structures", fontsize=8.5)
    panel_title(axA, "Structures per transporter")
    handles = [plt.Rectangle((0, 0), 1, 1, color=METHOD_PALETTE[c]) for c in pivot.columns]
    labels = ["Cryo-EM" if c == "ELECTRON MICROSCOPY" else
              "X-ray" if c == "X-RAY DIFFRACTION" else "NMR" for c in pivot.columns]
    axA.legend(handles, labels, fontsize=7.5, frameon=False, loc="lower right",
               handlelength=1.0, handletextpad=0.5, labelspacing=0.3,
               bbox_to_anchor=(1.0, 0.02))
    panel_tag(axA, "A")

    # Panel B: depositions per year (area) + cumulative distinct-SLC coverage (line).
    s = u.dropna(subset=["release_year"]).copy()
    s["release_year"] = s["release_year"].astype(int)
    years = list(range(int(s["release_year"].min()), int(s["release_year"].max()) + 1))
    per_year = s["release_year"].value_counts().reindex(years, fill_value=0)
    debuts = s.groupby("GENE")["release_year"].min().value_counts().reindex(years, fill_value=0)
    cumulative = debuts.cumsum()
    axB.fill_between(years, per_year.values, color=NEUTRAL, alpha=0.85, zorder=2,
                     linewidth=0)
    axB.set_ylim(bottom=0)
    axB.set_xlim(min(years), max(years))
    axB.set_ylabel("Structures / year", fontsize=8.5, color=MUTED)
    axB.set_xlabel("Release year", fontsize=8.5)
    axB.tick_params(axis="y", labelsize=8, colors=MUTED, length=0)
    axB.tick_params(axis="x", labelsize=8, length=0)
    lo = (years[0] + 4) // 5 * 5
    axB.set_xticks(list(range(lo, years[-1] + 1, 5)))
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)
    ax2 = axB.twinx()
    ax2.plot(years, cumulative.values, color=PRIMARY, linewidth=2.2, zorder=4)
    ax2.set_ylabel("Cumulative SLCs solved", color=PRIMARY, fontsize=8.5)
    ax2.tick_params(axis="y", colors=PRIMARY, labelsize=8, length=0)
    ax2.set_ylim(bottom=0)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_color(PRIMARY)
    ax2.annotate(f"{int(cumulative.iloc[-1])}", (years[-1], cumulative.iloc[-1]),
                 textcoords="offset points", xytext=(-2, 4), ha="right",
                 fontsize=8, weight="bold", color=PRIMARY)
    panel_title(axB, "Depositions & cumulative coverage")
    panel_tag(axB, "B")

    fig.tight_layout(w_pad=3.0)
    return fig


def fig_timeline_distributions(df):
    """2x2 timeline distributions as clean stacked-area charts (A-C) plus a
    resolution-by-method bar panel (D). No inset pies — one chart type per panel."""
    u = df.drop_duplicates("PDB ID").copy()
    u["release_year"] = u["release_year"].astype("Int64")
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.4))
    (axA, axB), (axC, axD) = axes

    def area_by_year(col, order, palette, ax, title, tag, legend_labels=None):
        s = u.dropna(subset=["release_year", col]).copy()
        s["release_year"] = s["release_year"].astype(int)
        years = list(range(int(s["release_year"].min()), int(s["release_year"].max()) + 1))
        cols = [c for c in order if c in s[col].unique()]
        pivot = (s.groupby(["release_year", col]).size().unstack(fill_value=0)
                    .reindex(index=years, fill_value=0).reindex(columns=cols, fill_value=0))
        stacked_area(ax, years, pivot, [palette[c] for c in pivot.columns])
        ax.set_xlabel("Release year", fontsize=8.5)
        ax.set_ylabel("Structures", fontsize=8.5)
        ax.tick_params(labelsize=8)
        panel_title(ax, title)
        names = legend_labels or list(pivot.columns)
        handles = [plt.Rectangle((0, 0), 1, 1, color=palette[c]) for c in pivot.columns]
        ax.legend(handles, names, fontsize=7.5, frameon=False, loc="upper left",
                  handlelength=1.0, handletextpad=0.5, labelspacing=0.3)
        panel_tag(ax, tag)

    # A: eukaryotic vs prokaryotic source
    area_by_year("organism_type", list(ORGANISM_TYPE_PALETTE), ORGANISM_TYPE_PALETTE,
                 axA, "Source organism type", "A")
    # B: determination method
    area_by_year("METHOD", order_methods(list(METHOD_PALETTE)), METHOD_PALETTE,
                 axB, "Determination method", "B",
                 legend_labels=["Cryo-EM", "X-ray", "NMR"])

    # C: resolution range by year (stacked bars — clearer than an area of near-
    # identical blue bands for reading per-year magnitudes)
    s = u.dropna(subset=["release_year"]).copy()
    s["release_year"] = s["release_year"].astype(int)
    yearsC = list(range(int(s["release_year"].min()), int(s["release_year"].max()) + 1))
    bin_order = [f"{b}-{b+1}" for b in range(1, 9)]
    s2 = s.dropna(subset=["resolution_bin"])
    pivotC = (s2.groupby(["release_year", "resolution_bin"], observed=True).size()
                 .unstack(fill_value=0).reindex(index=yearsC, fill_value=0)
                 .reindex(columns=bin_order, fill_value=0))
    pivotC = pivotC.loc[:, (pivotC != 0).any(axis=0)]
    colorsC = seq_palette(len(pivotC.columns))
    x = np.arange(len(yearsC))
    bottom = np.zeros(len(yearsC))
    for col, color in zip(pivotC.columns, colorsC):
        vals = pivotC[col].to_numpy(dtype=float)
        axC.bar(x, vals, bottom=bottom, color=color, width=0.86,
                edgecolor="white", linewidth=0.3, zorder=3)
        bottom += vals
    axC.set_xlim(-0.6, len(yearsC) - 0.4)
    axC.set_ylim(bottom=0)
    ticks = [i for i, y in enumerate(yearsC) if y % 5 == 0]   # sparse 5-year labels
    axC.set_xticks(ticks)
    axC.set_xticklabels([str(yearsC[i]) for i in ticks])
    for sp in ("top", "right"):
        axC.spines[sp].set_visible(False)
    axC.set_xlabel("Release year", fontsize=8.5)
    axC.set_ylabel("Structures", fontsize=8.5)
    axC.tick_params(labelsize=8, length=0)
    panel_title(axC, "Resolution range")
    handlesC = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colorsC]
    axC.legend(handlesC, [f"{c} Å" for c in pivotC.columns], fontsize=6.8,
               frameon=False, loc="upper left", ncol=2, handlelength=1.0,
               handletextpad=0.4, labelspacing=0.25, columnspacing=0.9)
    panel_tag(axC, "C")

    # D: resolution by method (clean vertical stacked bars)
    sub = u.dropna(subset=["resolution_bin"])
    pivotD = (sub.groupby(["resolution_bin", "METHOD"], observed=True).size()
                 .unstack(fill_value=0))
    pivotD = pivotD[order_methods(pivotD.columns)]
    x = np.arange(len(pivotD.index))
    bottom = np.zeros(len(pivotD.index))
    for c in pivotD.columns:
        vals = pivotD[c].to_numpy(dtype=float)
        axD.bar(x, vals, bottom=bottom, color=METHOD_PALETTE.get(c, NEUTRAL),
                width=0.74, edgecolor="white", linewidth=0.5, zorder=3)
        bottom += vals
    axD.set_xticks(x)
    axD.set_xticklabels([str(b) for b in pivotD.index], rotation=0, fontsize=7.5)
    axD.set_xlabel("Resolution (Å)", fontsize=8.5)
    axD.set_ylabel("Structures", fontsize=8.5)
    axD.tick_params(axis="y", labelsize=8, length=0)
    axD.tick_params(axis="x", length=0)
    for sp in ("top", "right"):
        axD.spines[sp].set_visible(False)
    panel_title(axD, "Resolution by method")
    handlesD = [plt.Rectangle((0, 0), 1, 1, color=METHOD_PALETTE[c]) for c in pivotD.columns]
    axD.legend(handlesD, ["Cryo-EM" if c == "ELECTRON MICROSCOPY" else "X-ray"
                          for c in pivotD.columns], fontsize=7.5, frameon=False,
               loc="upper right", handlelength=1.0, handletextpad=0.5, labelspacing=0.3)
    panel_tag(axD, "D")

    fig.tight_layout(w_pad=3.0, h_pad=3.2)
    return fig


def fig_sources_and_leaders(df):
    """(A) top source organisms; (B) group leaders — both as horizontal bars so the
    long organism/author names stay upright and readable."""
    u = df.drop_duplicates("PDB ID")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    # Panel A: top source organisms (italic species names down the left)
    org = u["SOURCE ORGANISM"].value_counts().head(10)
    org_labels = [s if len(str(s)) <= 26 else str(s)[:24] + "…" for s in org.index]
    hbar(axA, org_labels, org.values, italic=True, fontsize=7.5, value_fontsize=7.5)
    axA.set_xlabel("Structures", fontsize=8.5)
    panel_title(axA, "Source organisms")
    panel_tag(axA, "A")

    # Panel B: group leaders (senior authors)
    auth = u.dropna(subset=["senior_author"])["senior_author"].value_counts().head(10)
    hbar(axB, list(auth.index), auth.values, fontsize=7.5, value_fontsize=7.5)
    axB.set_xlabel("Structures", fontsize=8.5)
    panel_title(axB, "Group leaders")
    panel_tag(axB, "B")

    fig.tight_layout(w_pad=2.5)
    return fig


def fig_workflow(df):
    """Methods schematic: the query → PDB → mine → organize → visualize → interpret
    pipeline (an SLC analog of the ABC-transporter workflow figure).

    Laid out top-to-bottom with a single vertical cursor so every connector is
    [short arrow] then [stage label] in its own band — the arrow never runs through
    the label (the bug in the first version)."""
    u = df.drop_duplicates("PDB ID")
    n_struct, n_genes = u["PDB ID"].nunique(), df["GENE"].nunique()

    fig, ax = plt.subplots(figsize=(6.6, 9.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    cx = 0.5
    BOXFILL = "#F4F6F8"
    cur = [0.972]   # mutable cursor = top edge of the next element to draw

    def connector(label):
        """A short downward arrow followed, in a separate band below it, by the
        blue stage label. Advances the cursor past both."""
        cur[0] -= 0.018                      # gap below the block just drawn
        a0, a1 = cur[0], cur[0] - 0.034      # arrow span
        ax.annotate("", xy=(cx, a1), xytext=(cx, a0),
                    arrowprops=dict(arrowstyle="-|>", color="#9AA5B1", lw=1.7,
                                    shrinkA=0, shrinkB=0))
        ly = a1 - 0.030                      # label sits clearly below the arrowhead
        ax.text(cx, ly, label, ha="center", va="center", color=OI_BLUE,
                fontsize=10.5, weight="bold")
        cur[0] = ly - 0.026                  # clear the label before the next block

    def pill(text, *, w, h=0.058, fill=OI_BLUE, fg="white", fs=12):
        cy = cur[0] - h / 2
        ax.add_patch(mpatches.FancyBboxPatch(
            (cx - w / 2, cy - h / 2), w, h,
            boxstyle="round,pad=0.004,rounding_size=0.022",
            facecolor=fill, edgecolor="none", zorder=3, clip_on=False))
        ax.text(cx, cy, text, ha="center", va="center", color=fg, fontsize=fs,
                weight="bold", zorder=4)
        cur[0] -= h
        return cy

    # 1. Query --------------------------------------------------------------
    pill("SLC transporter gene set", w=0.50, fs=12.5)
    connector("Search")

    # 2. RCSB PDB -----------------------------------------------------------
    h = 0.092
    cy = cur[0] - h / 2
    ax.add_patch(mpatches.FancyBboxPatch(
        (cx - 0.21, cy - h / 2), 0.42, h,
        boxstyle="round,pad=0.004,rounding_size=0.022",
        facecolor=BOXFILL, edgecolor=SPINE, lw=1.1, zorder=3))
    ax.text(cx, cy + 0.014, "RCSB PDB", ha="center", va="center", color=INK,
            fontsize=13.5, weight="bold", zorder=4)
    ax.text(cx, cy - 0.018, "Search + Data APIs", ha="center", va="center",
            color=MUTED, fontsize=9.5, zorder=4)
    cur[0] -= h
    connector("Data mining")

    # 3. Variables mined (boxed for consistency) ----------------------------
    variables = ["Family", "Source organism", "Resolution", "Determination method",
                 "Group leader", "Oligomeric state", "Publication year"]
    h = 0.078
    cy = cur[0] - h / 2
    ax.add_patch(mpatches.FancyBboxPatch(
        (cx - 0.40, cy - h / 2), 0.80, h,
        boxstyle="round,pad=0.004,rounding_size=0.022",
        facecolor="white", edgecolor=SPINE, lw=1.0, zorder=3))
    ax.text(cx, cy + 0.017, "   ·   ".join(variables[:4]), ha="center",
            va="center", color=INK, fontsize=9.5, zorder=4)
    ax.text(cx, cy - 0.017, "   ·   ".join(variables[4:]), ha="center",
            va="center", color=INK, fontsize=9.5, zorder=4)
    cur[0] -= h
    connector("Data organization")

    # 4. Table schema -------------------------------------------------------
    cols = ["SLC gene", "PDB ID", "Year", "Organism", "Prok/Euk",
            "Method", "Resolution", "Oligomer"]
    tw, th = 0.94, 0.052
    cy = cur[0] - th / 2
    x0, cwid = cx - tw / 2, tw / len(cols)
    for i, c in enumerate(cols):
        ax.add_patch(mpatches.Rectangle(
            (x0 + i * cwid, cy - th / 2), cwid, th,
            facecolor=BOXFILL if i % 2 == 0 else "#E7ECF1",
            edgecolor=SPINE, lw=0.7, zorder=3))
        ax.text(x0 + (i + 0.5) * cwid, cy, c, ha="center", va="center",
                color=TEXT, fontsize=7.0, weight="medium", zorder=4)
    cur[0] -= th
    ax.text(cx, cur[0] - 0.016, f"one row per (gene, structure)  ·  {n_struct} "
            f"structures, {n_genes} genes", ha="center", va="center",
            color=FAINT, fontsize=8)
    cur[0] -= 0.034
    connector("Structural superposition")

    # 4b. Structural fold classification (US-align / TM-score) --------------
    h = 0.082
    cy = cur[0] - h / 2
    ax.add_patch(mpatches.FancyBboxPatch(
        (cx - 0.40, cy - h / 2), 0.80, h,
        boxstyle="round,pad=0.004,rounding_size=0.022",
        facecolor=BOXFILL, edgecolor=SPINE, lw=1.0, zorder=3))
    ax.text(cx, cy + 0.020, "Structural fold classification", ha="center",
            va="center", color=INK, fontsize=11.5, weight="bold", zorder=4)
    ax.text(cx, cy - 0.006, "US-align superposition vs canonical fold references",
            ha="center", va="center", color=MUTED, fontsize=8.8, zorder=4)
    ax.text(cx, cy - 0.026, "TM-score ≥ 0.5  →  assigned fold", ha="center",
            va="center", color=OI_BLUE, fontsize=8.8, weight="bold", zorder=4)
    cur[0] -= h
    connector("Data visualization")

    # 5. Visualization glyphs ----------------------------------------------
    h = 0.066
    _workflow_icons(ax, cur[0] - h / 2)
    cur[0] -= h
    connector("Data interpretation")

    # 6. Outcome -----------------------------------------------------------
    out_h = 0.058
    cy_out = pill("Trends in SLC structural biology", w=0.58, h=out_h,
                  fill=OI_ORANGE, fs=12)

    # Fit the view to the actual content so the top and bottom pills are never
    # clipped by the [0,1] axes box (the earlier cut-off bug).
    ax.set_ylim(cy_out - out_h / 2 - 0.015, 1.015)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.995, bottom=0.005)
    return fig


def _workflow_icons(ax, y, h=0.054):
    """Four tiny chart glyphs (bars, timeline, donut, trend) on light L-shaped axes,
    centred vertically on `y` — the 'data visualization' step of the workflow."""
    import matplotlib.lines as mlines
    centers = [0.20, 0.40, 0.60, 0.80]
    w = 0.115
    ybot, ytop = y - h / 2, y + h / 2

    def axes_glyph(c):
        left = c - w / 2
        ax.add_line(mlines.Line2D([left, left], [ybot, ytop], color=SPINE,
                                  lw=1.0, zorder=3))
        ax.add_line(mlines.Line2D([left, c + w / 2], [ybot, ybot], color=SPINE,
                                  lw=1.0, zorder=3))

    def line(c, ys, color):
        xs = np.linspace(c - w / 2 + 0.012, c + w / 2 - 0.008, len(ys))
        yy = [ybot + 0.006 + v * (h - 0.012) for v in ys]
        ax.add_line(mlines.Line2D(xs, yy, color=color, lw=1.8, zorder=4,
                                  solid_capstyle="round", solid_joinstyle="round"))

    for c in centers:
        axes_glyph(c)
    # bars
    c = centers[0]
    bx = np.linspace(c - w / 2 + 0.014, c + w / 2 - 0.014, 4)
    for xi, hh in zip(bx, [0.30, 0.55, 0.78, 1.0]):
        ax.add_patch(mpatches.Rectangle((xi - 0.006, ybot + 0.004), 0.012,
                                        hh * (h - 0.010), facecolor=OI_BLUE,
                                        edgecolor="none", zorder=4))
    # timeline / rising trace
    line(centers[1], [0.10, 0.22, 0.16, 0.5, 0.72, 1.0], OI_ORANGE)
    # donut
    c = centers[2]
    ax.add_patch(mpatches.Wedge((c, y), 0.020, 0, 360, width=0.009,
                                facecolor=NEUTRAL, edgecolor="none", zorder=4))
    ax.add_patch(mpatches.Wedge((c, y), 0.020, 90, 250, width=0.009,
                                facecolor=OI_BLUE, edgecolor="none", zorder=5))
    # trend line
    line(centers[3], [0.18, 0.36, 0.30, 0.62, 0.84, 1.0], OI_GREEN)


def figure_list():
    """Ordered (internal_name, fn, plos_title, interpretation) registry.

    `plos_title` is the manuscript caption title (<= 15 words, PLOS limit); it is
    written to figure_captions.md rather than drawn into the TIFF. `interpretation`
    is a plain-language reading of what the figure shows (the takeaway), also
    written to figure_captions.md. Entries whose internal name starts with 'S' are
    emitted as supplementary (S#_Fig).
    """
    return [
        ("fig0_workflow",               fig_workflow,
         "Workflow for studying the expansion of SLC structural biology",
         "Schematic of the analysis pipeline: a curated SLC gene set is queried "
         "against the RCSB PDB through its Search and Data APIs; per-structure "
         "variables (family, source organism, resolution, method, group leader, "
         "oligomeric state, year) are mined and organized into one row per "
         "(gene, structure); the table is then visualized and interpreted."),
        ("fig1_structures_timeline",    fig_structures_and_timeline,
         "Number of SLC structures per transporter and over time",
         "A handful of transporters dominate the structural record: SLC6A4, SLC1A1, "
         "XPR1, and SLC3A2 are the most-solved, and the high counts for SLC6A4 and "
         "SLC1A1 reflect attributed bacterial/archaeal homologs (panel A). Deposition "
         "activity was negligible before the mid-2000s and has risen steeply since "
         "~2019; the cumulative-coverage line shows that distinct-SLC coverage has "
         "broadened in step with that surge rather than just re-solving the same few "
         "transporters (panel B)."),
        ("fig2_debut_timeline",         fig_slc_debut_timeline,
         "Timeline of SLC transporters at their first solved structure",
         "Each transporter is named in the year its first structure was released and "
         "coloured by organism type, with the count of first-time debuts shown for each "
         "year; the field broadens from a trickle of early prokaryotic-homolog structures "
         "to a rapid recent run of first-time eukaryotic SLC structures."),
        ("fig3_timeline_distributions", fig_timeline_distributions,
         "Timeline distributions of source, method, and resolution",
         "Eukaryotic sources have dominated every recent year and now overwhelmingly "
         "drive new depositions (panel A). The field switched techniques around 2019: "
         "X-ray crystallography led through the 2010s, but cryo-EM has been the "
         "predominant method since (panel B), and the spread of sub-3 A structures "
         "tracks that switch (panel C). Yet the very highest resolutions remain X-ray's "
         "domain — the 1–2 A band is essentially all X-ray, while cryo-EM piles "
         "into the 3–4 A range, the single most common resolution regime (panel D)."),
        ("fig4_sources_leaders",        fig_sources_and_leaders,
         "Source organisms and group leaders in SLC structural biology",
         "Homo sapiens is the source for two-thirds of all structures, an order of "
         "magnitude ahead of any other organism, with the remainder spread across model "
         "mammals and thermophilic prokaryotes used as homologs (panel A; inset: ~79% "
         "eukaryotic vs ~21% prokaryotic). The field is also concentrated by laboratory: "
         "a small number of senior authors — led by Gouaux, Boudker, and Lee — "
         "account for a large share of the output (panel B)."),
        ("fig5_family",                 fig_family_pie,
         "Distribution of SLC structures across SLC families",
         "Structural effort is heavily skewed toward a few families: SLC6 (~15%) and "
         "SLC1 (~13%) alone account for roughly a quarter of all structures, and the ten "
         "largest families cover about two-thirds, leaving a long tail of ~27 sparsely "
         "characterized families pooled as 'Other'."),
        ("fig6_oligomeric",             fig_oligomeric_pie,
         "Distribution of SLC structures by oligomeric state",
         "Most SLC structures are small assemblies: monomers (38%) and dimers (34%) "
         "together make up roughly three-quarters of the dataset, with trimers (19%) and "
         "tetramers (7%) making up most of the rest and higher-order states rare."),
        ("fig7_folds",                  fig_fold_distribution,
         "Structural fold classification of SLC structures",
         "Classifying each structure by US-align structural alignment against a panel of "
         "canonical fold references recovers the expected dominance of a few folds: the "
         "LeuT/APC and MFS folds together account for the largest share of confidently "
         "assigned structures (TM-score >= 0.5), shown as both per-fold counts (bars) and "
         "the proportional split (donut)."),
        ("S1_fold_tmscore",             fig_fold_tmscore,
         "Structural-fold match quality by TM-score",
         "The TM-score of every structure to its assigned fold reference; most sit well "
         "above the 0.5 same-fold threshold (per-fold medians highlighted), confirming "
         "that confident assignments are genuine fold members rather than marginal matches."),
        ("S2_all_slcs",                 fig_all_slcs,
         "Structure counts for all SLC transporters in the dataset",
         "The full per-transporter ranking makes the skew explicit: a few transporters "
         "carry dozens of structures while the majority have only a handful, a "
         "long-tailed distribution typical of structure-focused effort."),
    ]


def render_all(df, plos):
    """Render every figure in either house-style PNG or PLOS TIFF mode."""
    global _PLOS
    _PLOS = plos
    setup_style(plos)
    main_i = sup_i = 0
    mapping = []
    for name, fn, title, interp in figure_list():
        fig = fn(df)
        if plos:
            if name.startswith("S"):
                sup_i += 1
                label = f"S{sup_i}_Fig"
            else:
                main_i += 1
                label = f"Fig{main_i}"
            save_fig_plos(fig, label)
            mapping.append((label, title, interp))
        else:
            save_fig(fig, name)
    return mapping


def write_captions(mapping):
    """Emit figure_captions.md: the title PLOS expects in the manuscript text plus a
    plain-language interpretation (the takeaway) for each figure."""
    lines = ["# Figure captions and interpretations (PLOS ONE)\n",
             "_Titles must appear in the manuscript text, not in the image files. Each "
             "title is <= 15 words; the italic paragraph beneath is a plain-language "
             "reading of the figure (the takeaway), suitable for the legend or Results._\n"]
    for label, title, interp in mapping:
        pretty = label.replace("_", " ")
        lines.append(f"**{pretty}. {title}.**\n")
        lines.append(f"_{interp}_\n")
    (PROJECT_DIR / "figure_captions.md").write_text("\n".join(lines))
    print(f"saved {PROJECT_DIR / 'figure_captions.md'}")


def main():
    df = load_data(INPUT_PATH, year_cutoff=YEAR_CUTOFF)
    print(f"Loaded {len(df)} rows ({df['GENE'].nunique()} unique genes, "
          f"{df['PDB ID'].nunique()} unique PDB IDs)")
    print("\n--- house-style PNGs (figures/) ---")
    render_all(df, plos=False)
    print("\n--- PLOS ONE TIFFs (figures_plos/) ---")
    mapping = render_all(df, plos=True)
    write_captions(mapping)


if __name__ == "__main__":
    main()
