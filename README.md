# SLC transporter structural meta-analysis

Code and data for a meta-analysis of solute carrier (SLC) transporter structures in
the RCSB Protein Data Bank. The pipeline scrapes per-gene structural metadata from
RCSB's APIs, classifies each structure's fold by structural alignment, and renders
the manuscript figures.

## Contents

| File | What it is |
|---|---|
| `rcsb_scraper_api.py` | Stage 1 — scrapes RCSB Search + Data APIs into `raw_output_web_scraping.{csv,xlsx}` |
| `superimpose.py` | Stage 2 — assigns each structure to a fold via US-align (TM-score) |
| `make_figures.py` | Stage 3 — renders all manuscript figures (house-style PNG + PLOS TIFF) |
| `SLC_Data_Cleaning_Graphs.Rmd` | legacy R figure path (kept for reference) |
| `new_slc_genes.csv` | input gene list (89 SLC genes queried) |
| `raw_output_web_scraping.{csv,xlsx}` | scraped dataset (one row per gene–structure) |
| `slc_folds_tmscore.csv` | per-structure fold call + TM-score + RMSD |
| `slc_folds_tmscore_allrefs.csv` | full TM-score matrix (every structure × every reference) |
| `raw_output_with_folds.csv` | scraped dataset with fold calls merged on |
| `slc_classification_verified.csv` | HGNC-verified gene/family classification |
| `figures/` | house-style figure PNGs |
| `figures_plos/` | 300-dpi PLOS submission TIFFs |
| `fonts/` | Inter `.ttf`s required by `make_figures.py` |

## Requirements

- **Python 3** with: `pandas`, `numpy`, `requests`, `openpyxl`, `matplotlib`, `seaborn`, `pillow`

  ```bash
  pip install pandas numpy requests openpyxl matplotlib seaborn pillow
  ```

- **US-align** (only for `superimpose.py`) — a single-file C++ program, built once.
  The binary is platform-specific, so it is **not** committed; build it yourself:

  ```bash
  git clone --depth 1 https://github.com/pylelab/USalign.git
  cd USalign && c++ -O3 -o USalign USalign.cpp
  ```

  Then either put `USalign` on your `PATH`, drop it next to `superimpose.py`, or point
  the `USALIGN` environment variable at it.

## Running the pipeline

The committed intermediate files let you run any stage on its own. For a full
rebuild, run them in order:

```bash
# Stage 1 — scrape RCSB (writes raw_output_web_scraping.{csv,xlsx})
python rcsb_scraper_api.py

# Stage 2 — classify each structure's fold (writes slc_folds_tmscore*.csv,
#           raw_output_with_folds.csv). Requires the US-align binary above and
#           downloads mmCIF files into ./cif_cache (git-ignored).
python superimpose.py

# Stage 3 — render every figure into figures/ and figures_plos/
python make_figures.py
```

Notes:
- Stage 2 must run before Stage 3 if you want the fold figure (Fig 7), which reads
  `slc_folds_tmscore.csv`. The repo ships this file, so Stage 3 works out of the box.
- `make_figures.py` registers each figure in `figure_list()`; main figures are named
  `Fig1…Fig7` and supplements `S1_Fig…S3_Fig` (the PLOS-mode TIFF names).

## Fold classification (`superimpose.py`)

Each structure's longest protein chain is structurally aligned (US-align) against a
panel of canonical fold references; the fold of the highest-scoring reference is
assigned, with TM-score normalised to the reference length. Calls at **TM-score ≥ 0.5**
(the standard same-fold threshold) are treated as confident; structures whose fold is
not in the reference panel, or that are not membrane transporters, are reported as
`Unassigned` rather than forced onto a non-matching reference.
