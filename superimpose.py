"""
superimpose.py
--------------
Assigns each SLC structure to a canonical structural fold and scores how well it
matches, using a *real* structural alignment (US-align / TM-align) rather than the
naive Cα index-pairing of the original version.

For every unique PDB in the scraped dataset this script:
  1. downloads the mmCIF (so large cryo-EM entries, which have no legacy .pdb, are
     not silently dropped),
  2. extracts the longest protein chain (the transporter; skips Fabs/nanobodies/tags),
  3. structurally aligns it against each canonical fold-reference structure with
     US-align, which establishes residue correspondence first and reports a
     length-normalised TM-score,
  4. assigns the fold = reference with the highest TM-score (normalised by the
     reference length), and records that score, the aligned-core RMSD, and whether
     the match is confident (TM-score >= 0.5, the standard "same fold" threshold).

Why this replaces the old approach: superposition RMSD is only meaningful between
*corresponding* residues. The previous code paired Cα atoms by file order
(`ref_cas[:n]` vs `mob_cas[:n]`), which measures sequence ordering, not 3-D shape —
the same structure reversed scored 26 A and two unrelated folds returned a
plausible-looking ~25 A with no error. TM-align/US-align solve the correspondence
problem and give a fold-diagnostic, length-independent score.

Output columns (one row per unique PDB):
  PDB ID, FOLD, FOLD_REF_PDB, TM_SCORE_TO_FOLD_REF, RMSD_TO_FOLD_REF,
  ALIGNED_LENGTH, CHAIN_LENGTH, FOLD_CONFIDENT

Requirements:
  - Python 3 with pandas (+ openpyxl for the .xlsx).  No biopython needed.
  - The `USalign` binary. Build once (single C++ file):
        git clone --depth 1 https://github.com/pylelab/USalign.git
        cd USalign && c++ -O3 -o USalign USalign.cpp
    then put it on PATH, set the USALIGN env var, or drop it next to this script.

Usage:
    python superimpose.py
"""

import os
import re
import sys
import time
import shutil
import subprocess
import urllib.request
import urllib.error
from collections import defaultdict

import pandas as pd

# ── Paths / config ────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV   = os.path.join(PROJECT_DIR, "raw_output_web_scraping.csv")
OUT_CSV     = os.path.join(PROJECT_DIR, "slc_folds_tmscore.csv")
OUT_WIDE    = os.path.join(PROJECT_DIR, "slc_folds_tmscore_allrefs.csv")  # per-reference scores
MERGED_CSV  = os.path.join(PROJECT_DIR, "raw_output_with_folds.csv")
CIF_CACHE   = os.path.join(PROJECT_DIR, "cif_cache")
CA_CACHE    = os.path.join(CIF_CACHE, "ca")

CONFIDENT_TM = 0.50   # standard TM-score threshold for "same fold"
REQUEST_PAUSE = 0.05  # be polite to the RCSB file server

# Locate the US-align binary: $USALIGN, then ./USalign, then PATH.
USALIGN = (os.environ.get("USALIGN")
           or (os.path.join(PROJECT_DIR, "USalign")
               if os.path.exists(os.path.join(PROJECT_DIR, "USalign")) else None)
           or shutil.which("USalign") or shutil.which("USalign.exe"))

# ── Canonical fold reference PDBs ─────────────────────────────────────────────
# Each maps a fold name -> a representative experimental structure for that fold.
# These are the well-established, literature-standard archetypes; verify before
# final submission. Folds with no curated reference yet are simply not scored
# (a structure of such a fold will get a low best-TM and FOLD_CONFIDENT == False).
FOLD_REFERENCES = {
    "MFS":  "4PYP",   # GLUT1 / SLC2A1 (major facilitator superfamily)
    "LeuT": "2A65",   # LeuT, Aquifex aeolicus (APC / LeuT fold)
    "MitC": "1OKC",   # ADP/ATP carrier (mitochondrial carrier fold)
    "UraA": "3QE7",   # UraA uracil transporter (SLC4/23/26 nucleobase-cation fold)
    "Glt":  "1XFH",   # GltPh, Pyrococcus (glutamate-transporter / elevator fold)
    "NhaA": "1ZCD",   # NhaA Na+/H+ antiporter
    "YiiP": "3H90",   # YiiP zinc transporter (CDF fold)
    "NCX":  "3V5U",   # NCX_Mj Na+/Ca2+ exchanger
    "AmtB": "1U7G",   # AmtB ammonium transporter
    "CNT":  "3TIJ",   # vcCNT concentrative nucleoside transporter
    "DMT":  "4GE6",   # nucleotide-sugar transporter (drug/metabolite transporter)
}


# ── mmCIF handling ────────────────────────────────────────────────────────────

def fetch_cif(pdb_id):
    """Download and cache the mmCIF for `pdb_id`. Returns path or None on failure."""
    pdb_id = pdb_id.upper()
    path = os.path.join(CIF_CACHE, f"{pdb_id}.cif")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        urllib.request.urlretrieve(url, path)
        time.sleep(REQUEST_PAUSE)
        return path
    except (urllib.error.HTTPError, urllib.error.URLError, OSError) as e:
        print(f"    [fetch] {pdb_id}: {e}")
        if os.path.exists(path):
            os.remove(path)
        return None


def _atom_site_rows(cif_path):
    """Yield the column index map and the ATOM/HETATM token rows of the _atom_site loop."""
    order, cols, rows, in_loop = [], {}, [], False
    with open(cif_path) as f:
        for ln in f:
            if ln.startswith("_atom_site."):
                order.append(ln.strip().split(".", 1)[1])
                cols = {c: i for i, c in enumerate(order)}
                in_loop = True
                continue
            if in_loop:
                if ln.startswith(("ATOM", "HETATM")):
                    toks = ln.split()
                    if len(toks) >= len(order):
                        rows.append(toks)
                elif rows and (ln.startswith(("#", "loop_", "_", "data_"))):
                    break
    return cols, rows


def _write_ca_pdb(recs, g, out):
    """Write a CA-only PDB from a chain's _atom_site rows; `g` reads a column."""
    with open(out, "w") as o:
        for i, r in enumerate(recs, 1):
            try:
                x, y, z = (float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                           float(g(r, "Cartn_z")))
            except ValueError:
                continue
            resn = g(r, "label_comp_id")[:3]
            try:
                seq = int(g(r, "auth_seq_id"))
            except ValueError:
                seq = i
            o.write(f"ATOM  {i:5d}  CA  {resn:>3} A{seq:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        o.write("END\n")


def chain_ca_files(pdb_id, min_len=120):
    """Extract every protein chain with >= `min_len` Cα as its own CA-only PDB.

    Returns a list of (chain_id, ca_pdb_path, n_residues) sorted longest-first.
    The longest chain is usually the transporter, but in heteromers it may not be
    (e.g. the SLC3A2/4F2hc heavy chain is longer than its SLC7 transporter partner),
    so the caller tries the longest chain first and falls back to the others to find
    the chain that actually carries a transporter fold. If every chain is shorter
    than `min_len` (e.g. an isolated resolved domain), the single longest is returned
    so the entry is still scored.
    """
    cif = fetch_cif(pdb_id)
    if cif is None:
        return []
    try:
        cols, rows = _atom_site_rows(cif)
    except Exception as e:
        print(f"    [parse] {pdb_id}: {e}")
        return []
    need = ("group_PDB", "label_atom_id", "auth_asym_id", "label_comp_id",
            "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z")
    if not all(k in cols for k in need):
        print(f"    [parse] {pdb_id}: missing _atom_site columns")
        return []

    def g(r, k):
        return r[cols[k]]

    chains = defaultdict(list)
    for r in rows:
        if g(r, "group_PDB") != "ATOM":          # polymer only, skip HETATM
            continue
        if g(r, "label_atom_id").strip('"') != "CA":
            continue
        chains[g(r, "auth_asym_id")].append(r)
    if not chains:
        return []

    ordered = sorted(chains, key=lambda c: len(chains[c]), reverse=True)
    out = []
    for ch in ordered:
        if len(chains[ch]) < min_len:
            break                                 # ordered desc, so the rest are smaller
        p = os.path.join(CA_CACHE, f"{pdb_id.upper()}_{ch}_ca.pdb")
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            _write_ca_pdb(chains[ch], g, p)
        out.append((ch, p, len(chains[ch])))
    if not out:                                   # all chains below min_len -> keep longest
        ch = ordered[0]
        p = os.path.join(CA_CACHE, f"{pdb_id.upper()}_{ch}_ca.pdb")
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            _write_ca_pdb(chains[ch], g, p)
        out.append((ch, p, len(chains[ch])))
    return out


def longest_chain_ca(pdb_id):
    """Longest protein chain only, as (ca_pdb_path, n_residues) — used for the
    single-chain fold reference structures."""
    cs = chain_ca_files(pdb_id, min_len=1)
    return (cs[0][1], cs[0][2]) if cs else (None, 0)


# ── US-align ──────────────────────────────────────────────────────────────────

_TM_RE   = re.compile(r"TM-score=\s*([\d.]+).*normalized by length of Structure_2")
_RMSD_RE = re.compile(r"RMSD=\s*([\d.]+)")
_ALN_RE  = re.compile(r"Aligned length=\s*(\d+)")


def usalign(mobile_ca, ref_ca):
    """Run US-align, normalised by the reference length. Returns (tm, rmsd, aligned_len)."""
    try:
        res = subprocess.run([USALIGN, mobile_ca, ref_ca, "-mol", "prot"],
                             capture_output=True, text=True, timeout=120)
    except (subprocess.TimeoutExpired, OSError) as e:
        print(f"    [usalign] {e}")
        return None, None, None
    out = res.stdout
    tm   = _TM_RE.search(out)
    rmsd = _RMSD_RE.search(out)
    aln  = _ALN_RE.search(out)
    return (float(tm.group(1)) if tm else None,
            float(rmsd.group(1)) if rmsd else None,
            int(aln.group(1)) if aln else None)


# ── Main ──────────────────────────────────────────────────────────────────────

def build_reference_cache():
    """Pre-extract the CA file for every fold reference. Returns {fold: (ref_pdb, ca_path)}."""
    print("Loading fold references...")
    ref_ca = {}
    for fold, ref_pdb in FOLD_REFERENCES.items():
        ca, n = longest_chain_ca(ref_pdb)
        if ca:
            ref_ca[fold] = (ref_pdb, ca)
            print(f"  {fold:5s} {ref_pdb}: {n} Ca")
        else:
            print(f"  {fold:5s} {ref_pdb}: FAILED to load — fold skipped")
    return ref_ca


def main():
    if not USALIGN:
        sys.exit("ERROR: US-align binary not found. Set $USALIGN, put `USalign` next "
                 "to this script, or add it to PATH. See the module docstring to build it.")
    os.makedirs(CIF_CACHE, exist_ok=True)
    os.makedirs(CA_CACHE, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)
    df = df[df["PDB ID"].notna() & (df["PDB ID"] != "NO STRUCTURE")]
    pdb_ids = sorted(df["PDB ID"].astype(str).str.upper().unique())
    print(f"{len(pdb_ids)} unique PDB structures to fold-classify "
          f"against {len(FOLD_REFERENCES)} references.\n")

    ref_ca = build_reference_cache()
    if not ref_ca:
        sys.exit("ERROR: no fold references could be loaded.")

    print(f"\nClassifying {len(pdb_ids)} structures...\n")
    rows, wide_rows = [], []
    for k, pdb in enumerate(pdb_ids, 1):
        chains = chain_ca_files(pdb)
        if not chains:
            print(f"  [{k}/{len(pdb_ids)}] {pdb}: could not load — skipped")
            continue

        # Try chains longest-first and keep the best (chain, reference) match,
        # stopping as soon as a chain confidently matches a fold. This lets a
        # heteromer fall through a non-transporter longest chain (e.g. the
        # SLC3A2/4F2hc heavy chain) to the chain that actually carries the fold
        # (its SLC7 light-chain partner).
        best = None   # (tm, rmsd, aln, fold, chain_id, chain_len, scores)
        for chain_id, ca_path, clen in chains:
            scores, cb = {}, None        # cb = best (tm, rmsd, aln, fold) for this chain
            for fold, (ref_pdb, ref_path) in ref_ca.items():
                tm, rmsd, aln = usalign(ca_path, ref_path)
                if tm is None:
                    continue
                scores[fold] = round(tm, 4)
                if cb is None or tm > cb[0]:
                    cb = (tm, rmsd, aln, fold)
            if cb is None:
                continue
            if best is None or cb[0] > best[0]:
                best = (*cb, chain_id, clen, scores)
            if best[0] >= CONFIDENT_TM:   # confident match found; no need for more chains
                break

        if best is None:
            print(f"  [{k}/{len(pdb_ids)}] {pdb}: no alignment produced — skipped")
            continue

        best_tm, best_rmsd, best_aln, best_fold, chain_used, chain_len, scores = best
        confident = best_tm >= CONFIDENT_TM
        rows.append({
            "PDB ID": pdb,
            "FOLD": best_fold if confident else "Unassigned",
            "FOLD_REF_PDB": FOLD_REFERENCES[best_fold],
            "TM_SCORE_TO_FOLD_REF": round(best_tm, 4),
            "RMSD_TO_FOLD_REF": best_rmsd,
            "ALIGNED_LENGTH": best_aln,
            "CHAIN_USED": chain_used,
            "CHAIN_LENGTH": chain_len,
            "FOLD_CONFIDENT": confident,
        })
        wide_rows.append({"PDB ID": pdb, **scores})

        flag = "" if confident else "  (low confidence)"
        print(f"  [{k}/{len(pdb_ids)}] {pdb}: {best_fold} "
              f"TM={best_tm:.3f} RMSD={best_rmsd} chain={chain_used}{flag}")

        if k % 25 == 0:   # checkpoint so a long run is resumable / inspectable
            pd.DataFrame(rows).to_csv(OUT_CSV, index=False)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    try:
        out.to_excel(OUT_CSV.replace(".csv", ".xlsx"), index=False)
    except Exception:
        pass
    pd.DataFrame(wide_rows).to_csv(OUT_WIDE, index=False)

    # Merge fold calls back onto every (gene, PDB) row of the scraped dataset.
    df["__u"] = df["PDB ID"].astype(str).str.upper()
    merged = df.merge(out, left_on="__u", right_on="PDB ID", how="left",
                      suffixes=("", "_fold")).drop(columns=["__u", "PDB ID_fold"])
    merged.to_csv(MERGED_CSV, index=False)

    n_conf = int(out["FOLD_CONFIDENT"].sum())
    print(f"\nDone. {len(out)} structures classified, {n_conf} confident "
          f"(TM >= {CONFIDENT_TM}).")
    print("By fold (confident only):")
    print(out[out["FOLD_CONFIDENT"]]["FOLD"].value_counts().to_string())
    print(f"\nWrote:\n  {OUT_CSV}\n  {OUT_WIDE}\n  {MERGED_CSV}")


if __name__ == "__main__":
    main()
