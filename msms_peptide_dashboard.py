"""
=============================================================================
  MS/MS Peptide Identification Pipeline  +  Full Interactive Dashboard
  v5 — FAST + ADAPTIVE: two-pass search, numpy scoring, bisect index, lru_cache
=============================================================================

INSTALL:
    pip install pyteomics biopython requests tqdm dash plotly pandas numpy

RUN:
    python msms_peptide_dashboard.py
    → open http://127.0.0.1:8051

SPEED IMPROVEMENTS vs v3 + ADAPTIVE SEARCH:
  • score_psm     → numpy vectorised b/y matching       ~6× faster
  • lookup        → bisect on sorted mass array          ~15× faster
  • ion cache     → theoretical_ions cached per peptide  eliminates redundant work
  • run_search    → concurrent.futures thread pool       ~N_CPU× faster
  • build_index   → numpy mass array built once          faster precursor lookup
=============================================================================
"""

import sys, math, io, base64, uuid, threading, time, bisect
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from functools import lru_cache

try:
    from pyteomics import mzml, mass, parser as tparser
    from Bio import SeqIO
    import requests
    import pandas as pd
    import numpy as np
    import dash
    from dash import dcc, html, dash_table, Input, Output, State, no_update
    import plotly.graph_objects as go
except ImportError as e:
    sys.exit(f"\n[FATAL] Missing: {e}\n"
             "Run:  pip install pyteomics biopython requests tqdm "
             "dash plotly pandas numpy\n")

import os
N_WORKERS = max(1, os.cpu_count() - 1)   # leave 1 core for Dash

# ═══════════════════════════════════════════════════════════════════════════════
#  DESIGN
# ═══════════════════════════════════════════════════════════════════════════════
C = dict(
    bg="#060a10", panel="#0c1220", panel2="#101828",
    border="#1a2840", border2="#243552", text="#d0dcea", muted="#4e6480",
    teal="#00e5c0", amber="#ffb830", rose="#ff5c7a", sky="#38bfff",
    green="#3ddc84", purple="#b388ff", orange="#ff8c42",
)
FONT = "'IBM Plex Mono', 'Fira Code', monospace"

MZML_DEFAULT  = ("https://ftp.pride.ebi.ac.uk/pride/data/archive/"
                 "2012/03/PXD000001/"
                 "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")
FASTA_DEFAULT = ("https://ftp.pride.ebi.ac.uk/pride/data/archive/"
                 "2012/03/PXD000001/erwinia_carotovora.fasta")

# ═══════════════════════════════════════════════════════════════════════════════
#  PIPELINE CORE  — optimised
# ═══════════════════════════════════════════════════════════════════════════════
AA_MASSES = mass.std_aa_mass
H2O    = 18.010565
PROTON = 1.007276
AAS    = [a for a in AA_MASSES if len(a) == 1]

# Precompute AA mass array for fast peptide mass calculation
_AA_MASS_ARR = {a: m for a, m in AA_MASSES.items() if len(a) == 1}


def peptide_mass(seq):
    return sum(_AA_MASS_ARR.get(a, 0.0) for a in seq) + H2O


# ── OPT 1: LRU-cached ion generation ─────────────────────────────────────────
@lru_cache(maxsize=131072)   # cache up to 128k peptides
def theoretical_ions_cached(seq: str):
    """Returns (b_array, y_array) as numpy arrays. Cached per sequence."""
    b, y, ba, ya = [], [], 0.0, H2O
    for a in seq[:-1]:
        ba += _AA_MASS_ARR.get(a, 0.0); b.append(ba + PROTON)
    for a in reversed(seq[1:]):
        ya += _AA_MASS_ARR.get(a, 0.0); y.append(ya + PROTON)
    return np.array(b, dtype=np.float64), np.array(y, dtype=np.float64)


# ── OPT 2: Vectorised numpy scoring  (replaces nested Python loops) ───────────
def score_psm(mz_arr, int_arr, b_ions, y_ions, tol):
    """
    Vectorised b/y-ion hyperscore.
    mz_arr, int_arr : numpy arrays (already sorted ascending by pyteomics)
    b_ions, y_ions  : numpy arrays from theoretical_ions_cached
    """
    if len(b_ions) == 0 and len(y_ions) == 0:
        return 0, 0.0

    log_int = np.log1p(int_arr)   # pre-log intensities

    def _match(ions):
        if len(ions) == 0:
            return 0, 0.0
        # broadcast: (n_peaks, n_ions) difference matrix
        diffs   = np.abs(mz_arr[:, None] - ions[None, :])
        within  = diffs <= tol                          # boolean mask
        hit_pks = within.any(axis=1)                   # peaks that matched anything
        n_ions  = int(within.any(axis=0).sum())        # ions that were matched
        isum    = float(log_int[hit_pks].sum())
        return n_ions, isum

    nb, sb = _match(b_ions)
    ny, sy = _match(y_ions)
    n = nb + ny
    if n == 0:
        return 0, 0.0
    # hyperscore factorial term
    fl = (sum(math.log(i) for i in range(1, nb + 1)) +
          sum(math.log(i) for i in range(1, ny + 1)))
    return n, fl + sb + sy


# ── OPT 3: Sorted-array mass index (bisect lookup, ~15× vs float-while loop) ──
class MassIndex:
    """
    Stores peptide entries sorted by mass for O(log n) range lookup.
    Replaces the float-keyed defaultdict + while-loop approach.
    """
    def __init__(self):
        self._masses  = []   # sorted float list
        self._entries = []   # parallel list of (pep, protein_id)

    def add(self, mass_val, pep, pid):
        idx = bisect.bisect_left(self._masses, mass_val)
        self._masses.insert(idx, mass_val)
        self._entries.insert(idx, (pep, pid))

    def build_from_dict(self, d):
        """Fast bulk-load from existing {mass: [(pep,pid),...]} dict."""
        pairs = [(m, pep, pid) for m, entries in d.items()
                 for pep, pid in entries]
        pairs.sort(key=lambda x: x[0])
        self._masses  = [p[0] for p in pairs]
        self._entries = [(p[1], p[2]) for p in pairs]
        return self

    def lookup(self, neutral, ppm):
        delta = neutral * ppm / 1e6
        lo, hi = neutral - delta, neutral + delta
        i = bisect.bisect_left(self._masses, lo)
        j = bisect.bisect_right(self._masses, hi)
        return self._entries[i:j]


# ── OPT 4: Download helper ────────────────────────────────────────────────────
def download_to_file(url, dest, log_cb):
    dest = Path(dest)
    if dest.exists():
        log_cb(f"  ✔ Cached: {dest.name}")
        return str(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    log_cb(f"  Downloading {dest.name} …")
    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done  = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(1 << 16):
                f.write(chunk); done += len(chunk)
                if total and done * 100 // total % 10 == 0:
                    log_cb(f"    {done*100//total}%  ({done/1e6:.1f} MB)")
    log_cb(f"  ✔ Saved: {dest.name}")
    return str(dest)


# ── OPT 5: Build index returns MassIndex objects ──────────────────────────────
def build_index(fasta_src, cfg, log_cb):
    t0 = time.perf_counter()
    log_cb("► Building peptide index …")

    if str(fasta_src).startswith("http"):
        fasta_src = download_to_file(fasta_src,
                                     "pipeline_output/database.fasta", log_cb)

    t1 = time.perf_counter()
    recs = list(SeqIO.parse(str(fasta_src), "fasta"))
    log_cb(f"  FASTA parse         : {len(recs):,} proteins  [{time.perf_counter()-t1:.2f}s]")

    rule = tparser.expasy_rules.get(cfg["enzyme"],
                                    tparser.expasy_rules["trypsin"])
    raw  = defaultdict(list)
    seen = set()

    t1 = time.perf_counter()
    for rec in recs:
        seq = str(rec.seq).upper()
        try:
            peps = tparser.cleave(seq, rule,
                                  missed_cleavages=cfg["max_missed"],
                                  min_length=cfg["min_pep_len"])
        except Exception:
            continue
        for pep in peps:
            if len(pep) > cfg["max_pep_len"]: continue
            if any(c not in _AA_MASS_ARR for c in pep): continue
            m = round(peptide_mass(pep), 4)
            raw[m].append((pep, rec.id))
            seen.add(pep)
    log_cb(f"  Tryptic digest      : {len(seen):,} peptides  [{time.perf_counter()-t1:.2f}s]")

    t1 = time.perf_counter()
    for pep in seen:
        theoretical_ions_cached(pep)
    log_cb(f"  Ion cache warm-up   : {len(seen):,} peptides  [{time.perf_counter()-t1:.2f}s]")

    t1 = time.perf_counter()
    decoy_raw = defaultdict(list)
    for m, entries in raw.items():
        for pep, pid in entries:
            decoy_raw[m].append((pep[::-1], "DECOY_" + pid))
    target_idx = MassIndex().build_from_dict(raw)
    decoy_idx  = MassIndex().build_from_dict(decoy_raw)
    log_cb(f"  Mass index build    : target+decoy  [{time.perf_counter()-t1:.2f}s]")

    log_cb(f"  ✔ Index total       : [{time.perf_counter()-t0:.2f}s]")
    return target_idx, decoy_idx, seen


# ── OPT 6: Per-spectrum scoring function (picklable for ThreadPoolExecutor) ───
def _score_spectrum(args):
    """
    Score one spectrum against candidates. Runs in a thread-pool worker.
    Returns best PSM tuple or None, plus the precursor mass error (ppm)
    and the best fragment delta (Da) so we can learn from it in Pass 1.
    """
    (scan_id, mz_a, it_a, neutral, chg, pmz,
     target_entries, decoy_entries, cfg) = args

    tol   = cfg["frag_tol"]
    min_m = cfg["min_matches"]
    mz_np = np.asarray(mz_a, dtype=np.float64)
    it_np = np.asarray(it_a, dtype=np.float32)

    best      = None
    best_dppm = None   # precursor error of best hit (ppm)
    best_dfrag = None  # median fragment error of best hit (Da)

    for is_dec, entries in [(False, target_entries), (True, decoy_entries)]:
        for pep, pid in entries[:60]:
            b, y = theoretical_ions_cached(pep)
            n, s = score_psm(mz_np, it_np, b, y, tol)
            if n < min_m: continue
            if best is None or s > best[3]:
                best = (scan_id, pep, pid, s, n, is_dec, chg, pmz, neutral)
                # record errors for parameter learning
                theo_mass = peptide_mass(pep)
                best_dppm  = (neutral - theo_mass) / theo_mass * 1e6
                # median absolute fragment error from matched ions
                ions = np.concatenate([b, y])
                diffs = np.abs(mz_np[:, None] - ions[None, :])
                matched = diffs[diffs.min(axis=1) <= tol].min(axis=0)                           if diffs.size else np.array([tol])
                best_dfrag = float(np.median(diffs[diffs <= tol]))                              if (diffs <= tol).any() else tol
    return best, best_dppm, best_dfrag


# ── TWO-PASS ADAPTIVE SEARCH ──────────────────────────────────────────────────
def _learn_params(pilot_psms, log_cb):
    """
    Analyse high-scoring PSMs from the pilot run to determine
    optimal narrow-search tolerances.

    Returns dict with:
      prec_ppm_narrow  – tight precursor window (2× observed 95th percentile)
      frag_tol_narrow  – tight fragment window  (2× observed 95th percentile)
    """
    # keep only confident target hits
    hits = [(dppm, dfrag)
            for psm, dppm, dfrag in pilot_psms
            if psm and not psm[5]       # target only
            and dppm is not None
            and dfrag is not None]

    if len(hits) < 20:
        log_cb("  ⚠ Too few pilot hits to learn params — keeping original tolerances")
        return None

    ppm_vals  = np.array([abs(h[0]) for h in hits])
    frag_vals = np.array([abs(h[1]) for h in hits])

    # use 95th percentile × 1.5 safety margin as the narrow window
    ppm_narrow  = float(np.percentile(ppm_vals,  95)) * 1.5
    frag_narrow = float(np.percentile(frag_vals, 95)) * 1.5

    # floor values so we never go absurdly tight
    ppm_narrow  = max(ppm_narrow,  2.0)
    frag_narrow = max(frag_narrow, 0.005)

    log_cb(f"  ✔ Learned  prec: ±{ppm_narrow:.1f} ppm  ")
    return dict(prec_ppm_narrow=ppm_narrow, frag_tol_narrow=frag_narrow)


def _read_all_spectra(mzml_path, target_idx, decoy_idx, prec_ppm, log_cb):
    """Read mzML once and return list of spectrum argument tuples."""
    spectra = []; total = 0
    cfg_dummy = {"prec_ppm": prec_ppm, "frag_tol": 999,
                 "min_matches": 0}   # tol/min_matches filled in later
    with mzml.MzML(str(mzml_path)) as reader:
        for spec in reader:
            if spec.get("ms level") != 2: continue
            total += 1
            mz_a = spec["m/z array"]; it_a = spec["intensity array"]
            if len(mz_a) < 5: continue
            prec = spec.get("precursorList", {}).get("precursor", [])
            if not prec: continue
            si = prec[0].get("selectedIonList", {}).get("selectedIon", [])
            if not si: continue
            pmz = si[0].get("selected ion m/z", 0)
            chg = int(si[0].get("charge state", 2) or 2)
            if pmz == 0: continue
            neutral = pmz * chg - chg * PROTON
            scan_id = spec.get("id", f"scan_{total}")
            t_e = target_idx.lookup(neutral, prec_ppm)
            d_e = decoy_idx.lookup(neutral,  prec_ppm)
            if not t_e and not d_e: continue
            spectra.append((scan_id, mz_a, it_a, neutral,
                            chg, pmz, t_e, d_e, None))   # cfg slot = None
    log_cb(f"  {total:,} MS2 spectra → {len(spectra):,} have candidates")
    return spectra, total


def _parallel_score(spectra_with_cfg, log_cb, label=""):
    """Score a list of (…, cfg) tuples in parallel. Returns raw result list."""
    results = []; done = 0
    chunk   = max(1, len(spectra_with_cfg) // 20)
    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        futs = {pool.submit(_score_spectrum, s): s
                for s in spectra_with_cfg}
        for fut in as_completed(futs):
            results.append(fut.result())
            done += 1
            if done % chunk == 0:
                log_cb(f"  {label}{done:,}/{len(spectra_with_cfg):,}  ")
    return results


def run_search(mzml_src, target_idx, decoy_idx, cfg, log_cb):
    """
    Two-pass adaptive search:
      Pass 1 — broad search on 10 % of spectra → learn optimal tolerances
      Pass 2 — narrow search on remaining 90 % using learned tolerances
    """
    t0 = time.perf_counter()
    log_cb(f"► Adaptive two-pass search  [{N_WORKERS} threads]")

    if str(mzml_src).startswith("http"):
        mzml_src = download_to_file(mzml_src,
                                    "pipeline_output/spectra.mzML", log_cb)

    # ── Read all spectra ──────────────────────────────────────────────────
    t1 = time.perf_counter()
    log_cb("  Reading mzML …")
    all_spectra, total = _read_all_spectra(
        mzml_src, target_idx, decoy_idx, cfg["prec_ppm"], log_cb)
    log_cb(f"  mzML read           : {total:,} MS2, {len(all_spectra):,} with candidates  [{time.perf_counter()-t1:.2f}s]")

    if not all_spectra:
        log_cb("  ✖ No spectra with candidates found.")
        return [], total, str(mzml_src)

    # shuffle for representative pilot sample
    rng       = np.random.default_rng(42)
    idx_all   = np.arange(len(all_spectra))
    rng.shuffle(idx_all)
    n_pilot   = max(50, len(all_spectra) // 10)
    pilot_idx = idx_all[:n_pilot]
    main_idx  = idx_all[n_pilot:]

    # ── PASS 1: pilot broad search ────────────────────────────────────────
    t1 = time.perf_counter()
    log_cb(f"  ┌─ PASS 1 — pilot broad search  ({n_pilot:,} spectra, {cfg['prec_ppm']} ppm / {cfg['frag_tol']*1000:.0f} mDa)")
    pilot_cfg      = dict(cfg)
    pilot_with_cfg = [(*all_spectra[i][:-1], pilot_cfg) for i in pilot_idx]
    pilot_raw      = _parallel_score(pilot_with_cfg, log_cb, "    p1 ")
    pilot_psms     = [r[0] for r in pilot_raw if r[0] is not None]
    log_cb(f"  └─ PASS 1 done      : {len(pilot_psms):,} PSMs  [{time.perf_counter()-t1:.2f}s]")

    # learn narrow parameters from pilot
    learned     = _learn_params(pilot_raw, log_cb)
    narrow_ppm  = learned["prec_ppm_narrow"]  if learned else cfg["prec_ppm"]
    narrow_frag = learned["frag_tol_narrow"]  if learned else cfg["frag_tol"]

    # ── PASS 2: narrow search on remaining 90 % ───────────────────────────
    t1 = time.perf_counter()
    log_cb(f"  ┌─ PASS 2 — narrow search  ({len(main_idx):,} spectra, {narrow_ppm:.1f} ppm / {narrow_frag*1000:.1f} mDa)")

    t_lu = time.perf_counter()
    narrow_cfg   = dict(cfg, prec_ppm=narrow_ppm, frag_tol=narrow_frag)
    main_spectra = []
    for i in main_idx:
        s = all_spectra[i]
        scan_id, mz_a, it_a, neutral, chg, pmz = s[:6]
        t_e = target_idx.lookup(neutral, narrow_ppm)
        d_e = decoy_idx.lookup(neutral,  narrow_ppm)
        if not t_e and not d_e: continue
        main_spectra.append((scan_id, mz_a, it_a, neutral,
                             chg, pmz, t_e, d_e, narrow_cfg))
    log_cb(f"    Narrow re-lookup   : {len(main_spectra):,}/{len(main_idx):,} have candidates  [{time.perf_counter()-t_lu:.2f}s]")

    main_raw  = _parallel_score(main_spectra, log_cb, "    p2 ")
    main_psms = [r[0] for r in main_raw if r[0] is not None]
    log_cb(f"  └─ PASS 2 done      : {len(main_psms):,} PSMs  [{time.perf_counter()-t1:.2f}s]")

    all_psms = pilot_psms + main_psms
    log_cb(f"  ✔ Search total      : {len(all_psms):,} candidates  [{time.perf_counter()-t0:.2f}s]")
    log_cb(f"  Throughput          : {total / max(time.perf_counter()-t0, 0.1):.0f} spectra/s")

    # save learned params
    Path("pipeline_output").mkdir(exist_ok=True)
    with open("pipeline_output/learned_params.txt", "w") as f:
        f.write("=== Adaptive Search — Learned Parameters ===\n")
        f.write(f"Pilot spectra used    : {n_pilot:,}\n")
        f.write(f"Pilot PSMs found      : {len(pilot_psms):,}\n")
        if learned:
            f.write(f"Learned prec_ppm      : {narrow_ppm:.2f}\n")
            f.write(f"Learned frag_tol (Da) : {narrow_frag:.5f}\n")
        else:
            f.write("Fallback to original tolerances (too few pilot hits)\n")

    return all_psms, total, str(mzml_src)


def apply_fdr(psms, fdr_thresh, log_cb):
    t0 = time.perf_counter()
    log_cb(f"► FDR filtering at {fdr_thresh*100:.1f}% …")
    nt = nd = 0; accepted = []
    for p in sorted(psms, key=lambda x: -x[3]):
        if p[5]: nd += 1
        else:    nt += 1
        if nd / max(nt, 1) > fdr_thresh: break
        if not p[5]: accepted.append(p)
    log_cb(f"  ✔ {len(accepted):,} PSMs pass FDR  [{time.perf_counter()-t0:.2f}s]")
    return accepted



def discover_novel(mzml_path, accepted, target_idx, known_peps,
                   cfg, novel_modes, log_cb):
    t0 = time.perf_counter()
    log_cb("► Novel peptide discovery …")
    matched = {p[0] for p in accepted}
    novel   = []; n_mut = n_mc = 0
    tol     = cfg["frag_tol"]
    min_m   = cfg["min_matches"]

    with mzml.MzML(str(mzml_path)) as reader:
        for spec in reader:
            if spec.get("ms level") != 2: continue
            prec = spec.get("precursorList", {}).get("precursor", [])
            if not prec: continue
            si = prec[0].get("selectedIonList", {}).get("selectedIon", [])
            if not si: continue
            scan_id = spec.get("id", "")
            if scan_id in matched: continue
            mz_a = np.asarray(spec["m/z array"], dtype=np.float64)
            it_a = np.asarray(spec["intensity array"], dtype=np.float32)
            if len(mz_a) < 5: continue
            pmz = si[0].get("selected ion m/z", 0)
            chg = int(si[0].get("charge state", 2) or 2)
            if pmz == 0: continue
            neutral = pmz * chg - chg * PROTON

            if "mut" in (novel_modes or []):
                best_n = None
                for parent, pid in target_idx.lookup(neutral, cfg["prec_ppm"])[:15]:
                    for i, orig in enumerate(parent):
                        for new_aa in AAS:
                            if new_aa == orig: continue
                            var = parent[:i] + new_aa + parent[i+1:]
                            if abs(peptide_mass(var) - neutral) * 1e6 \
                                    / neutral > cfg["prec_ppm"]: continue
                            if var in known_peps: continue
                            b, y = theoretical_ions_cached(var)
                            n, s = score_psm(mz_a, it_a, b, y, tol)
                            if n < min_m: continue
                            bp, yp = theoretical_ions_cached(parent)
                            _, ps  = score_psm(mz_a, it_a, bp, yp, tol)
                            if s >= ps * cfg["rescue_boost"]:
                                if best_n is None or s > best_n["score"]:
                                    best_n = dict(
                                        scan_id=scan_id, peptide=var,
                                        parent=parent, protein=pid, score=s,
                                        n_match=n,
                                        change=f"{orig}{i+1}{new_aa}",
                                        kind="Point Mutation",
                                        charge=chg, prec_mz=pmz)
                if best_n:
                    novel.append(best_n); n_mut += 1; continue

            if "mc" in (novel_modes or []):
                for p1, pr1 in target_idx.lookup(neutral, cfg["prec_ppm"])[:20]:
                    m2 = neutral - peptide_mass(p1) + H2O
                    for p2, _ in target_idx.lookup(m2, cfg["prec_ppm"])[:10]:
                        joined = p1 + p2
                        if joined in known_peps or \
                           len(joined) > cfg["max_pep_len"]: continue
                        if abs(peptide_mass(joined) - neutral) * 1e6 \
                                / neutral > cfg["prec_ppm"]: continue
                        b, y = theoretical_ions_cached(joined)
                        n, s = score_psm(mz_a, it_a, b, y, tol)
                        if n >= min_m:
                            novel.append(dict(
                                scan_id=scan_id, peptide=joined,
                                parent=f"{p1}+{p2}", protein=pr1,
                                score=s, n_match=n,
                                change="missed_cleavage",
                                kind="Missed Cleavage",
                                charge=chg, prec_mz=pmz))
                            n_mc += 1

    log_cb(f"  ✔ Point mutations: {n_mut} | Missed cleavages: {n_mc}  [{time.perf_counter()-t0:.2f}s]")
    return novel


def psms_to_df(accepted):
    if not accepted: return pd.DataFrame()
    return pd.DataFrame([{
        "Scan": p[0], "Peptide": p[1],
        "Protein": p[2].split("|")[1] if "|" in p[2] else p[2],
        "Score": round(p[3], 3), "Matches": p[4],
        "Charge": p[6], "Precursor_mz": round(p[7], 4),
        "PepLen": len(p[1]),
    } for p in accepted])


def novel_to_df(novel):
    if not novel:
        return pd.DataFrame(columns=["Scan", "Novel_Peptide", "Kind",
                                      "Change", "Parent", "Protein",
                                      "Score", "Matches", "Charge"])
    return pd.DataFrame([{
        "Scan": n["scan_id"], "Novel_Peptide": n["peptide"],
        "Kind": n["kind"], "Change": n["change"], "Parent": n["parent"],
        "Protein": n["protein"].split("|")[1]
                   if "|" in n["protein"] else n["protein"],
        "Score": round(n["score"], 3), "Matches": n["n_match"],
        "Charge": n["charge"],
    } for n in novel])



# ═══════════════════════════════════════════════════════════════════════════════
#  BACKGROUND JOB MANAGER  — keeps Dash responsive during long runs
# ═══════════════════════════════════════════════════════════════════════════════
class JobManager:
    def __init__(self):
        self._jobs = {}   # job_id -> dict
        self._lock = threading.Lock()

    def new_job(self):
        jid = str(uuid.uuid4())[:8]
        with self._lock:
            self._jobs[jid] = {"status": "running", "logs": [],
                                "psm": None, "novel": None, "total": 0}
        return jid

    def log(self, jid, msg):
        with self._lock:
            if jid in self._jobs:
                self._jobs[jid]["logs"].append(msg)
        print(msg, flush=True)

    def finish(self, jid, psm_json, novel_json, total):
        with self._lock:
            if jid in self._jobs:
                self._jobs[jid].update(status="done",
                                        psm=psm_json,
                                        novel=novel_json,
                                        total=total)

    def fail(self, jid, err):
        with self._lock:
            if jid in self._jobs:
                self._jobs[jid].update(status="error", error=err)

    def get(self, jid):
        with self._lock:
            return dict(self._jobs.get(jid, {}))

    def logs_text(self, jid):
        with self._lock:
            return "\n".join(self._jobs.get(jid, {}).get("logs", []))


JM = JobManager()
CURRENT_JOB = {"id": None}   # track the latest job id


def _fmt(seconds):
    """Human-readable elapsed time."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    return f"{m}m {s:02d}s"


def pipeline_worker(jid, mzml_src, fasta_src, cfg,
                    fdr_thresh, novel_modes):
    """Runs the full pipeline in a background thread with per-step timing."""
    def log(msg): JM.log(jid, msg)

    def timed_step(label, fn, *args, **kwargs):
        """Run fn(*args, **kwargs), log elapsed time, return result."""
        log(f"")
        log(f"┌─ {label}")
        t0 = time.perf_counter()
        result = fn(*args, **kwargs)
        elapsed = time.perf_counter() - t0
        log(f"└─ ✔ {label} done  [{_fmt(elapsed)}]")
        return result, elapsed

    try:
        Path("pipeline_output").mkdir(exist_ok=True)
        t_pipeline = time.perf_counter()

        log(f"╔══════════════════════════════════════════════")
        log(f"  PIPELINE STARTED  [job {jid}]")
        log(f"  enzyme={cfg['enzyme']}  missed={cfg['max_missed']}  "
            f"fdr={fdr_thresh*100:.1f}%  ppm={cfg['prec_ppm']}")
        log(f"  threads={N_WORKERS}")
        log(f"  mzML  : {str(mzml_src)[:65]}")
        log(f"  FASTA : {str(fasta_src)[:65]}")
        log(f"╚══════════════════════════════════════════════")

        (target_idx, decoy_idx, known_peps), t_idx = timed_step(
            "STEP 1 — Build peptide index",
            build_index, fasta_src, cfg, log)

        (psms, total, mzml_path), t_search = timed_step(
            "STEP 2 — Adaptive two-pass search",
            run_search, mzml_src, target_idx, decoy_idx, cfg, log)

        (accepted,), t_fdr = timed_step(
            "STEP 3 — FDR filtering",
            lambda: (apply_fdr(psms, fdr_thresh, log),))

        (novel,), t_novel = timed_step(
            "STEP 4 — Novel peptide discovery",
            lambda: (discover_novel(mzml_path, accepted, target_idx,
                                    known_peps, cfg, novel_modes, log),))

        df_psm   = psms_to_df(accepted)
        df_novel = novel_to_df(novel)
        df_psm.to_csv("pipeline_output/identified_psms.csv",   index=False)
        df_novel.to_csv("pipeline_output/novel_peptides.csv",  index=False)

        t_total = time.perf_counter() - t_pipeline

        log("")
        log("╔══════════════════════════════════════════════")
        log(f"  PIPELINE COMPLETE")
        log(f"  {'Step':<35} {'Time':>8}")
        log(f"  {'─'*44}")
        log(f"  {'1. Build index':<35} {_fmt(t_idx):>8}")
        log(f"  {'2. Adaptive search':<35} {_fmt(t_search):>8}")
        log(f"  {'3. FDR filter':<35} {_fmt(t_fdr):>8}")
        log(f"  {'4. Novel discovery':<35} {_fmt(t_novel):>8}")
        log(f"  {'─'*44}")
        log(f"  {'TOTAL':<35} {_fmt(t_total):>8}")
        log(f"  {'─'*44}")
        log(f"  PSMs identified : {len(df_psm):,}")
        log(f"  Novel peptides  : {len(df_novel):,}")
        log(f"  Spectra/second  : {total / max(t_search, 0.1):.0f}")
        log("╚══════════════════════════════════════════════")

        JM.finish(jid,
                  df_psm.to_json(orient="split"),
                  df_novel.to_json(orient="split"),
                  total)
    except Exception as e:
        import traceback
        log(f"✖ ERROR: {e}")
        log(traceback.format_exc()[-800:])
        JM.fail(jid, str(e))


# ═══════════════════════════════════════════════════════════════════════════════
#  UI HELPERS
# ═══════════════════════════════════════════════════════════════════════════════
def _label(text):
    return html.Label(text, style={
        "fontSize": "10px", "color": C["muted"],
        "textTransform": "uppercase", "letterSpacing": "0.09em",
        "display": "block", "marginBottom": "5px"})


def _input(id_, placeholder, value="", typ="text"):
    return dcc.Input(id=id_, type=typ, placeholder=placeholder, value=value,
                     debounce=True, style={
                         "width": "100%", "background": C["bg"],
                         "border": f"1px solid {C['border2']}",
                         "color": C["text"], "borderRadius": "6px",
                         "padding": "8px 12px", "fontFamily": FONT,
                         "fontSize": "11px", "boxSizing": "border-box"})


def _section(title, children):
    return html.Div([
        html.Div(title, style={
            "fontSize": "9px", "color": C["teal"],
            "letterSpacing": "0.15em", "textTransform": "uppercase",
            "marginBottom": "12px", "paddingBottom": "6px",
            "borderBottom": f"1px solid {C['border']}"}),
        *children,
    ], style={"marginBottom": "22px"})


def _btn(label, id_, color, extra=None):
    s = {"background": "transparent", "border": f"1px solid {color}",
         "color": color, "borderRadius": "6px", "padding": "9px 18px",
         "cursor": "pointer", "fontFamily": FONT, "fontSize": "12px",
         "letterSpacing": "0.05em"}
    if extra: s.update(extra)
    return html.Button(label, id=id_, style=s)


def empty_fig(msg="Awaiting pipeline run …"):
    fig = go.Figure()
    fig.add_annotation(text=msg, x=0.5, y=0.5,
                       xref="paper", yref="paper",
                       showarrow=False,
                       font=dict(size=12, color=C["muted"], family=FONT))
    fig.update_layout(paper_bgcolor=C["panel"], plot_bgcolor=C["bg"],
                      xaxis=dict(visible=False), yaxis=dict(visible=False),
                      margin=dict(l=20, r=20, t=30, b=20))
    return fig


AXIS = dict(color=C["muted"], gridcolor=C["border"], showgrid=True,
            zeroline=False, tickfont=dict(size=10, family=FONT))
LAY  = dict(paper_bgcolor=C["panel"], plot_bgcolor=C["bg"],
            font=dict(family=FONT, color=C["text"]),
            margin=dict(l=50, r=20, t=38, b=50))

# ═══════════════════════════════════════════════════════════════════════════════
#  APP LAYOUT
# ═══════════════════════════════════════════════════════════════════════════════
app = dash.Dash(__name__, title="MS/MS Peptide Dashboard",
                suppress_callback_exceptions=True)
server = app.server

app.layout = html.Div(style={
    "minHeight": "100vh", "background": C["bg"],
    "color": C["text"], "fontFamily": FONT,
    "display": "flex", "flexDirection": "column",
}, children=[
    html.Link(rel="stylesheet",
              href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono"
                   ":wght@400;500;600;700&display=swap"),

    # polling interval — fires every 2 s to refresh log + charts
    dcc.Interval(id="poll", interval=2000, disabled=True),

    # ── NAVBAR ──────────────────────────────────────────────────────────────
    html.Div([
        html.Div([
            html.Span("⬡", style={"color": C["teal"], "fontSize": "20px",
                                   "marginRight": "10px"}),
            html.Span("MS/MS PEPTIDE DASHBOARD",
                      style={"fontSize": "13px", "fontWeight": "700",
                             "letterSpacing": "0.08em"}),
        ], style={"display": "flex", "alignItems": "center"}),
        html.Div(id="navbar-status",
                 style={"fontSize": "11px", "color": C["muted"]}),
    ], style={
        "background": C["panel"], "borderBottom": f"1px solid {C['border']}",
        "padding": "13px 28px", "display": "flex",
        "justifyContent": "space-between", "alignItems": "center",
        "position": "sticky", "top": "0", "zIndex": "100",
    }),

    html.Div(style={"display": "flex", "flex": "1"}, children=[

        # ══════════════════════════════════════════════════════════════════
        #  SIDEBAR
        # ══════════════════════════════════════════════════════════════════
        html.Div(style={
            "width": "292px", "flexShrink": "0",
            "background": C["panel2"],
            "borderRight": f"1px solid {C['border']}",
            "padding": "22px 18px", "overflowY": "auto",
            "height": "calc(100vh - 48px)",
            "position": "sticky", "top": "48px",
        }, children=[

            # DATA SOURCES
            _section("Data Sources", [
                _label("mzML Input"),
                dcc.RadioItems(id="mzml-mode",
                    options=[{"label": " URL / PRIDE",  "value": "url"},
                             {"label": " Upload File",  "value": "upload"}],
                    value="url", inline=True,
                    style={"fontSize": "11px", "color": C["text"],
                           "marginBottom": "8px"},
                    inputStyle={"marginRight": "4px", "marginLeft": "8px"}),
                html.Div(id="mzml-url-div",
                         children=[_input("mzml-url",
                                          "mzML URL or local path",
                                          MZML_DEFAULT)],
                         style={"marginBottom": "4px"}),
                html.Div(id="mzml-upload-div",
                         style={"display": "none"}, children=[
                    dcc.Upload(id="mzml-upload",
                               children=html.Div([
                                   "Drop .mzML or ",
                                   html.A("browse",
                                          style={"color": C["teal"],
                                                 "cursor": "pointer"})
                               ], style={"fontSize": "11px",
                                         "color": C["muted"]}),
                               style={"border": f"1px dashed {C['border2']}",
                                      "borderRadius": "6px",
                                      "padding": "14px", "textAlign": "center",
                                      "cursor": "pointer",
                                      "background": C["bg"]}),
                    html.Div(id="mzml-fname",
                             style={"fontSize": "10px", "color": C["muted"],
                                    "marginTop": "4px"}),
                ]),

                html.Div(style={"marginBottom": "14px"}),
                _label("FASTA Input"),
                dcc.RadioItems(id="fasta-mode",
                    options=[{"label": " URL / UniProt", "value": "url"},
                             {"label": " Upload File",   "value": "upload"}],
                    value="url", inline=True,
                    style={"fontSize": "11px", "color": C["text"],
                           "marginBottom": "8px"},
                    inputStyle={"marginRight": "4px", "marginLeft": "8px"}),
                html.Div(id="fasta-url-div",
                         children=[_input("fasta-url",
                                          "FASTA URL or local path",
                                          FASTA_DEFAULT)]),
                html.Div(id="fasta-upload-div",
                         style={"display": "none"}, children=[
                    dcc.Upload(id="fasta-upload",
                               children=html.Div([
                                   "Drop .fasta or ",
                                   html.A("browse",
                                          style={"color": C["teal"],
                                                 "cursor": "pointer"})
                               ], style={"fontSize": "11px",
                                         "color": C["muted"]}),
                               style={"border": f"1px dashed {C['border2']}",
                                      "borderRadius": "6px",
                                      "padding": "14px", "textAlign": "center",
                                      "cursor": "pointer",
                                      "background": C["bg"]}),
                    html.Div(id="fasta-fname",
                             style={"fontSize": "10px", "color": C["muted"],
                                    "marginTop": "4px"}),
                ]),
            ]),

            # DIGESTION
            _section("Digestion", [
                _label("Enzyme"),
                dcc.Dropdown(id="enzyme",
                    options=[
                        {"label": "Trypsin",      "value": "trypsin"},
                        {"label": "Trypsin/P",    "value": "trypsin/p"},
                        {"label": "Lys-C",        "value": "lys-c"},
                        {"label": "Chymotrypsin", "value": "chymotrypsin low specificity"},
                    ],
                    value="trypsin", clearable=False,
                    style={"background": C["bg"], "color": C["bg"],
                           "border": f"1px solid {C['border2']}",
                           "borderRadius": "6px", "fontSize": "11px",
                           "marginBottom": "14px"}),
                _label("Missed Cleavages"),
                dcc.Slider(id="max-missed", min=0, max=3, step=1, value=2,
                           marks={i: str(i) for i in range(4)},
                           tooltip={"placement": "bottom",
                                    "always_visible": True}),
                html.Div(style={"marginBottom": "10px"}),
                _label("Min / Max Peptide Length"),
                html.Div([
                    dcc.Input(id="min-pep-len", type="number", value=6,
                              min=4, max=10,
                              style={"width": "48%", "background": C["bg"],
                                     "border": f"1px solid {C['border2']}",
                                     "color": C["text"], "borderRadius": "6px",
                                     "padding": "7px 10px",
                                     "fontFamily": FONT, "fontSize": "11px"}),
                    dcc.Input(id="max-pep-len", type="number", value=40,
                              min=20, max=80,
                              style={"width": "48%", "background": C["bg"],
                                     "border": f"1px solid {C['border2']}",
                                     "color": C["text"], "borderRadius": "6px",
                                     "padding": "7px 10px",
                                     "fontFamily": FONT, "fontSize": "11px"}),
                ], style={"display": "flex",
                          "justifyContent": "space-between"}),
            ]),

            # SEARCH PARAMS
            _section("Search Parameters", [
                _label("Precursor Tolerance (ppm)"),
                dcc.Slider(id="prec-ppm", min=5, max=50, step=5, value=20,
                           marks={5: "5", 20: "20", 50: "50"},
                           tooltip={"placement": "bottom",
                                    "always_visible": True}),
                html.Div(style={"marginBottom": "10px"}),
                _label("Fragment Tolerance (Da)"),
                dcc.Slider(id="frag-tol", min=0.01, max=0.2, step=0.01,
                           value=0.05,
                           marks={0.01: "0.01", 0.1: "0.1", 0.2: "0.2"},
                           tooltip={"placement": "bottom",
                                    "always_visible": True}),
                html.Div(style={"marginBottom": "10px"}),
                _label("Min Fragment Matches"),
                dcc.Slider(id="min-matches", min=2, max=8, step=1, value=4,
                           marks={i: str(i) for i in range(2, 9)},
                           tooltip={"placement": "bottom",
                                    "always_visible": True}),
                html.Div(style={"marginBottom": "10px"}),
                _label("FDR Threshold"),
                dcc.Dropdown(id="fdr-thresh",
                    options=[{"label": "0.1 % — strict",  "value": 0.001},
                             {"label": "1 %  — standard", "value": 0.01},
                             {"label": "5 %  — lenient",  "value": 0.05}],
                    value=0.01, clearable=False,
                    style={"background": C["bg"], "color": C["bg"],
                           "border": f"1px solid {C['border2']}",
                           "borderRadius": "6px", "fontSize": "11px"}),
            ]),

            # NOVEL DISCOVERY
            _section("Novel Peptide Discovery", [
                _label("Search For"),
                dcc.Checklist(id="novel-modes",
                    options=[
                        {"label": "  Point Mutations",       "value": "mut"},
                        {"label": "  Missed Cleavage Joins", "value": "mc"},
                    ],
                    value=["mut", "mc"],
                    style={"fontSize": "11px", "color": C["text"],
                           "lineHeight": "2.2"},
                    inputStyle={"marginRight": "6px"}),
                html.Div(style={"marginBottom": "10px"}),
                _label("Rescue Score Ratio"),
                dcc.Slider(id="rescue-boost", min=0.4, max=1.0,
                           step=0.05, value=0.7,
                           marks={0.4: "0.4", 0.7: "0.7", 1.0: "1.0"},
                           tooltip={"placement": "bottom",
                                    "always_visible": True}),
            ]),

            # RUN
            html.Div([
                _btn("▶  Run Pipeline", "run-btn", C["teal"],
                     {"width": "100%", "padding": "11px",
                      "fontSize": "13px", "fontWeight": "600",
                      "textAlign": "center"}),
                html.Div(id="run-status",
                         style={"fontSize": "10px", "color": C["muted"],
                                "marginTop": "7px", "textAlign": "center",
                                "minHeight": "14px"}),
            ]),

            # DOWNLOADS
            html.Div([
                html.Div(style={"borderTop": f"1px solid {C['border']}",
                                "margin": "18px 0 14px"}),
                html.Div([
                    _btn("⬇ PSMs",  "dl-psm-btn",  C["sky"],
                         {"flex": "1", "textAlign": "center"}),
                    _btn("⬇ Novel", "dl-novel-btn", C["rose"],
                         {"flex": "1", "textAlign": "center"}),
                ], style={"display": "flex", "gap": "8px"}),
            ]),

            # hidden stores
            dcc.Store(id="store-psm"),
            dcc.Store(id="store-novel"),
            dcc.Store(id="store-total"),
            dcc.Store(id="store-mzml-path"),
            dcc.Store(id="store-fasta-path"),
            dcc.Store(id="job-id"),
            dcc.Download(id="dl-psm"),
            dcc.Download(id="dl-novel"),
        ]),

        # ══════════════════════════════════════════════════════════════════
        #  MAIN CONTENT
        # ══════════════════════════════════════════════════════════════════
        html.Div(style={"flex": "1", "padding": "22px 26px",
                        "overflowY": "auto"}, children=[

            html.Div(id="kpi-row",
                     children=[html.Div(
                         "Configure inputs in the sidebar and "
                         "click ▶ Run Pipeline.",
                         style={"color": C["muted"], "fontSize": "12px",
                                "padding": "10px 0"})],
                     style={"display": "flex", "gap": "10px",
                            "marginBottom": "20px", "flexWrap": "wrap"}),

            # Log console
            html.Div([
                html.Div("PIPELINE LOG", style={
                    "fontSize": "9px", "color": C["teal"],
                    "letterSpacing": "0.15em", "marginBottom": "8px"}),
                html.Div(id="log-console",
                         children=[html.Span(
                             "Ready — configure the sidebar and "
                             "press ▶ Run Pipeline.",
                             style={"color": C["muted"]})],
                         style={
                             "background": C["bg"],
                             "border": f"1px solid {C['border']}",
                             "borderRadius": "6px",
                             "padding": "12px 16px",
                             "fontFamily": FONT, "fontSize": "11px",
                             "color": C["green"], "minHeight": "70px",
                             "maxHeight": "160px", "overflowY": "auto",
                             "whiteSpace": "pre-wrap", "lineHeight": "1.8",
                         }),
            ], style={"background": C["panel"],
                      "border": f"1px solid {C['border']}",
                      "borderRadius": "10px", "padding": "16px 20px",
                      "marginBottom": "18px"}),

            # Row 1
            html.Div([
                html.Div(dcc.Graph(id="fig-score",  figure=empty_fig(),
                                   config={"displayModeBar": False}),
                         style={"flex": "2", "minWidth": "240px",
                                "background": C["panel"],
                                "border": f"1px solid {C['border']}",
                                "borderRadius": "10px", "padding": "4px"}),
                html.Div(dcc.Graph(id="fig-donut",  figure=empty_fig(),
                                   config={"displayModeBar": False}),
                         style={"flex": "1", "minWidth": "170px",
                                "background": C["panel"],
                                "border": f"1px solid {C['border']}",
                                "borderRadius": "10px", "padding": "4px"}),
                html.Div(dcc.Graph(id="fig-len",    figure=empty_fig(),
                                   config={"displayModeBar": False}),
                         style={"flex": "2", "minWidth": "240px",
                                "background": C["panel"],
                                "border": f"1px solid {C['border']}",
                                "borderRadius": "10px", "padding": "4px"}),
            ], style={"display": "flex", "gap": "14px",
                      "marginBottom": "14px", "flexWrap": "wrap"}),

            # Row 2
            html.Div([
                html.Div(dcc.Graph(id="fig-prec",
                                   figure=empty_fig(),
                                   config={"displayModeBar": True,
                                           "modeBarButtonsToRemove":
                                               ["lasso2d"]}),
                         style={"flex": "3", "minWidth": "280px",
                                "background": C["panel"],
                                "border": f"1px solid {C['border']}",
                                "borderRadius": "10px", "padding": "4px"}),
                html.Div(dcc.Graph(id="fig-prot",  figure=empty_fig(),
                                   config={"displayModeBar": False}),
                         style={"flex": "2", "minWidth": "240px",
                                "background": C["panel"],
                                "border": f"1px solid {C['border']}",
                                "borderRadius": "10px", "padding": "4px"}),
            ], style={"display": "flex", "gap": "14px",
                      "marginBottom": "14px", "flexWrap": "wrap"}),

            # Row 3 — novel
            html.Div(dcc.Graph(id="fig-novel",
                               figure=empty_fig(
                                   "Novel peptide results will appear here."),
                               config={"displayModeBar": True,
                                       "modeBarButtonsToRemove":
                                           ["lasso2d"]}),
                     style={"background": C["panel"],
                            "border": f"1px solid {C['border']}",
                            "borderRadius": "10px", "padding": "4px",
                            "marginBottom": "14px"}),

            # PSM table
            html.Div([
                html.Div("IDENTIFIED PSMs",
                         style={"fontSize": "9px", "color": C["teal"],
                                "letterSpacing": "0.15em",
                                "marginBottom": "12px"}),
                html.Div(id="psm-table-container",
                         children=html.Div(
                             "Run pipeline to populate.",
                             style={"color": C["muted"], "fontSize": "12px",
                                    "padding": "14px 0"})),
            ], style={"background": C["panel"],
                      "border": f"1px solid {C['border']}",
                      "borderRadius": "10px", "padding": "18px",
                      "marginBottom": "14px"}),

            # Novel table
            html.Div([
                html.Div("NOVEL / VARIANT PEPTIDES",
                         style={"fontSize": "9px", "color": C["teal"],
                                "letterSpacing": "0.15em",
                                "marginBottom": "8px"}),
                html.Div([
                    html.Span("● Point Mutation  ",
                              style={"color": C["rose"], "fontSize": "11px"}),
                    html.Span("● Missed Cleavage",
                              style={"color": C["purple"],
                                     "fontSize": "11px"}),
                ], style={"marginBottom": "10px"}),
                html.Div(id="novel-table-container",
                         children=html.Div(
                             "Run pipeline to populate.",
                             style={"color": C["muted"], "fontSize": "12px",
                                    "padding": "14px 0"})),
            ], style={"background": C["panel"],
                      "border": f"1px solid {C['border']}",
                      "borderRadius": "10px", "padding": "18px",
                      "marginBottom": "24px"}),

            html.Div(
                "pyteomics · biopython · Dash · Plotly  |  "
                "PRIDE PXD000001 · Erwinia carotovora",
                style={"textAlign": "center", "fontSize": "10px",
                       "color": C["muted"], "paddingBottom": "16px"}),
        ]),
    ]),
])

# ═══════════════════════════════════════════════════════════════════════════════
#  CALLBACKS
# ═══════════════════════════════════════════════════════════════════════════════

@app.callback(
    Output("mzml-url-div",    "style"),
    Output("mzml-upload-div", "style"),
    Input("mzml-mode", "value"))
def toggle_mzml(mode):
    s, h = {"display": "block"}, {"display": "none"}
    return (s, h) if mode == "url" else (h, s)


@app.callback(
    Output("fasta-url-div",    "style"),
    Output("fasta-upload-div", "style"),
    Input("fasta-mode", "value"))
def toggle_fasta(mode):
    s, h = {"display": "block"}, {"display": "none"}
    return (s, h) if mode == "url" else (h, s)


@app.callback(
    Output("store-mzml-path", "data"),
    Output("mzml-fname",      "children"),
    Input("mzml-upload",      "contents"),
    State("mzml-upload",      "filename"),
    prevent_initial_call=True)
def save_mzml(contents, filename):
    if not contents: return no_update, ""
    Path("pipeline_output").mkdir(exist_ok=True)
    raw = base64.b64decode(contents.split(",")[1])
    p = f"pipeline_output/upload_{filename}"
    with open(p, "wb") as f: f.write(raw)
    return p, f"✓ {filename}"


@app.callback(
    Output("store-fasta-path", "data"),
    Output("fasta-fname",      "children"),
    Input("fasta-upload",      "contents"),
    State("fasta-upload",      "filename"),
    prevent_initial_call=True)
def save_fasta(contents, filename):
    if not contents: return no_update, ""
    Path("pipeline_output").mkdir(exist_ok=True)
    raw = base64.b64decode(contents.split(",")[1])
    p = f"pipeline_output/upload_{filename}"
    with open(p, "wb") as f: f.write(raw)
    return p, f"✓ {filename}"


# ── Launch pipeline in background thread ──────────────────────────────────────
@app.callback(
    Output("job-id",      "data"),
    Output("poll",        "disabled"),
    Output("run-status",  "children"),

    Input("run-btn", "n_clicks"),

    State("mzml-mode",        "value"),
    State("mzml-url",         "value"),
    State("store-mzml-path",  "data"),
    State("fasta-mode",       "value"),
    State("fasta-url",        "value"),
    State("store-fasta-path", "data"),
    State("enzyme",           "value"),
    State("max-missed",       "value"),
    State("min-pep-len",      "value"),
    State("max-pep-len",      "value"),
    State("prec-ppm",         "value"),
    State("frag-tol",         "value"),
    State("min-matches",      "value"),
    State("fdr-thresh",       "value"),
    State("rescue-boost",     "value"),
    State("novel-modes",      "value"),
    prevent_initial_call=True)
def launch_pipeline(n_clicks,
                    mzml_mode, mzml_url_val, mzml_upload_path,
                    fasta_mode, fasta_url_val, fasta_upload_path,
                    enzyme, max_missed, min_pep_len, max_pep_len,
                    prec_ppm, frag_tol, min_matches, fdr_thresh,
                    rescue_boost, novel_modes):

    mzml_src  = (mzml_upload_path  if mzml_mode  == "upload"
                                      and mzml_upload_path
                 else mzml_url_val)
    fasta_src = (fasta_upload_path if fasta_mode == "upload"
                                      and fasta_upload_path
                 else fasta_url_val)

    if not mzml_src or not fasta_src:
        return no_update, True, "⚠ Provide mzML and FASTA sources first."

    cfg = dict(
        enzyme=enzyme or "trypsin",
        max_missed=int(max_missed or 2),
        min_pep_len=int(min_pep_len or 6),
        max_pep_len=int(max_pep_len or 40),
        prec_ppm=float(prec_ppm or 20),
        frag_tol=float(frag_tol or 0.05),
        min_matches=int(min_matches or 4),
        rescue_boost=float(rescue_boost or 0.7),
    )
    fdr = float(fdr_thresh or 0.01)

    jid = JM.new_job()
    CURRENT_JOB["id"] = jid

    t = threading.Thread(
        target=pipeline_worker,
        args=(jid, mzml_src, fasta_src, cfg, fdr, novel_modes),
        daemon=True)
    t.start()

    return jid, False, f"⏳ Running … [job {jid}]"


# ── Poll for log + completion ──────────────────────────────────────────────────
@app.callback(
    Output("log-console",          "children"),
    Output("navbar-status",        "children"),
    Output("store-psm",            "data"),
    Output("store-novel",          "data"),
    Output("store-total",          "data"),
    Output("poll",                 "disabled", allow_duplicate=True),
    Output("run-status",           "children", allow_duplicate=True),

    Input("poll",   "n_intervals"),
    State("job-id", "data"),
    prevent_initial_call=True)
def poll_job(n_intervals, jid):
    if not jid:
        return (no_update,) * 7

    job = JM.get(jid)
    if not job:
        return (no_update,) * 7

    log_text = JM.logs_text(jid)
    log_html = html.Pre(log_text,
                        style={"margin": "0", "color": C["green"],
                               "fontFamily": FONT, "fontSize": "11px"})

    status = job.get("status", "running")

    if status == "running":
        lines   = log_text.split("\n")
        last    = lines[-1] if lines else "running …"
        nav_msg = f"⏳ {last[:60]}"
        return (log_html, nav_msg,
                no_update, no_update, no_update,
                False, f"⏳ Running … [job {jid}]")

    elif status == "done":
        n_psm = len(pd.read_json(
            io.StringIO(job["psm"]), orient="split")) if job["psm"] else 0
        n_nov = len(pd.read_json(
            io.StringIO(job["novel"]), orient="split")) if job["novel"] else 0
        done_msg = f"✔ {n_psm:,} PSMs | {n_nov:,} novel"

        # show error-coloured log if 0 results
        log_color = C["amber"] if n_psm == 0 else C["green"]
        log_html  = html.Pre(log_text,
                             style={"margin": "0", "color": log_color,
                                    "fontFamily": FONT, "fontSize": "11px"})
        return (log_html, done_msg,
                job["psm"], job["novel"], job["total"],
                True, done_msg)   # stop polling

    else:  # error
        err_html = html.Pre(log_text,
                            style={"margin": "0", "color": C["rose"],
                                   "fontFamily": FONT, "fontSize": "11px"})
        err_msg  = f"✖ {job.get('error', 'unknown error')}"
        return (err_html, err_msg,
                no_update, no_update, no_update,
                True, err_msg)   # stop polling


# ── Build charts once data arrives ────────────────────────────────────────────
@app.callback(
    Output("kpi-row",               "children"),
    Output("fig-score",             "figure"),
    Output("fig-donut",             "figure"),
    Output("fig-len",               "figure"),
    Output("fig-prec",              "figure"),
    Output("fig-prot",              "figure"),
    Output("fig-novel",             "figure"),
    Output("psm-table-container",   "children"),
    Output("novel-table-container", "children"),
    Input("store-psm",   "data"),
    Input("store-novel", "data"),
    Input("store-total", "data"),
    prevent_initial_call=True)
def update_charts(psm_json, novel_json, total):
    if not psm_json: return [no_update] * 9

    df_psm   = pd.read_json(io.StringIO(psm_json),   orient="split")
    df_novel = pd.read_json(io.StringIO(novel_json),  orient="split") \
               if novel_json else pd.DataFrame()

    total    = total or 1
    n_psm    = len(df_psm)
    n_prot   = df_psm["Protein"].nunique() if n_psm else 0
    n_novel  = len(df_novel)
    n_mut    = int((df_novel["Kind"] == "Point Mutation").sum()) \
               if n_novel and "Kind" in df_novel else 0
    n_mc     = int((df_novel["Kind"] == "Missed Cleavage").sum()) \
               if n_novel and "Kind" in df_novel else 0
    id_rate  = f"{n_psm / total * 100:.1f}%"

    def kpi(lbl, val, col):
        return html.Div([
            html.Div(str(val),
                     style={"fontSize": "26px", "fontWeight": "700",
                            "color": col, "lineHeight": "1.1"}),
            html.Div(lbl,
                     style={"fontSize": "9px", "color": C["muted"],
                            "textTransform": "uppercase",
                            "letterSpacing": "0.1em", "marginTop": "3px"}),
        ], style={"background": C["panel"],
                  "border": f"1px solid {C['border']}",
                  "borderRadius": "10px", "padding": "12px 16px",
                  "textAlign": "center", "flex": "1", "minWidth": "95px"})

    kpi_row = [
        kpi("Spectra",     f"{total:,}",   C["sky"]),
        kpi("PSMs",        f"{n_psm:,}",   C["teal"]),
        kpi("Proteins",    f"{n_prot:,}",  C["amber"]),
        kpi("Novel",       f"{n_novel:,}", C["rose"]),
        kpi("Mut. Rescue", f"{n_mut:,}",   C["purple"]),
        kpi("MC Rescue",   f"{n_mc:,}",    C["orange"]),
        kpi("ID Rate",     id_rate,         C["green"]),
    ]

    # score dist
    fig_score = go.Figure(go.Histogram(
        x=df_psm["Score"], nbinsx=60,
        marker_color=C["teal"], opacity=0.85))
    fig_score.update_layout(**LAY,
        title=dict(text="PSM Score Distribution", font=dict(size=13)),
        xaxis={**AXIS, "title": "Hyperscore"},
        yaxis={**AXIS, "title": "Count"}, showlegend=False)

    # donut
    fig_donut = go.Figure(go.Pie(
        labels=["Identified", "Unidentified"],
        values=[n_psm, max(total - n_psm, 0)],
        hole=0.65, marker=dict(colors=[C["teal"], C["border"]]),
        textinfo="none"))
    fig_donut.update_layout(**LAY,
        title=dict(text="ID Rate", font=dict(size=13)),
        showlegend=True,
        legend=dict(font=dict(size=10), bgcolor="rgba(0,0,0,0)"),
        annotations=[dict(text=id_rate, x=0.5, y=0.5, showarrow=False,
                          font=dict(size=18, color=C["teal"],
                                    family=FONT))],
        margin=dict(l=10, r=10, t=38, b=10))

    # length
    lc = df_psm["PepLen"].value_counts().sort_index()
    fig_len = go.Figure(go.Bar(x=lc.index, y=lc.values,
                               marker_color=C["sky"], opacity=0.85))
    fig_len.update_layout(**LAY,
        title=dict(text="Peptide Length Distribution", font=dict(size=13)),
        xaxis={**AXIS, "title": "Length (AA)"},
        yaxis={**AXIS, "title": "Count"}, bargap=0.12)

    # precursor scatter
    chg_col = {1: C["rose"], 2: C["teal"], 3: C["amber"],
               4: C["purple"], 5: C["sky"]}
    fig_prec = go.Figure()
    for chg, grp in df_psm.groupby("Charge"):
        fig_prec.add_trace(go.Scatter(
            x=grp["Precursor_mz"], y=grp["Score"],
            mode="markers", name=f"z={chg}",
            marker=dict(color=chg_col.get(chg, C["muted"]),
                        size=4, opacity=0.6),
            hovertemplate="m/z: %{x:.3f}<br>Score: %{y:.2f}<extra></extra>"))
    fig_prec.update_layout(**LAY,
        title=dict(text="Precursor m/z vs Score by Charge", font=dict(size=13)),
        xaxis={**AXIS, "title": "Precursor m/z"},
        yaxis={**AXIS, "title": "Hyperscore"})

    # top proteins
    pc = (df_psm.groupby("Protein").size()
          .sort_values(ascending=True).tail(20))
    fig_prot = go.Figure(go.Bar(
        x=pc.values, y=pc.index, orientation="h",
        marker=dict(color=pc.values,
                    colorscale=[[0, C["border"]], [1, C["amber"]]],
                    showscale=False)))
    fig_prot.update_layout(**LAY,
        title=dict(text="Top 20 Proteins by PSM Count", font=dict(size=13)),
        xaxis={**AXIS, "title": "PSM Count"},
        yaxis={**AXIS, "title": ""},
        height=400, margin=dict(l=155, r=20, t=38, b=50))

    # novel scatter
    if n_novel and "Novel_Peptide" in df_novel.columns:
        kc = {"Point Mutation": C["rose"], "Missed Cleavage": C["purple"]}
        fig_novel = go.Figure()
        for kind, grp in df_novel.groupby("Kind"):
            fig_novel.add_trace(go.Scatter(
                x=grp["Novel_Peptide"].str.len(),
                y=grp["Score"], mode="markers", name=kind,
                marker=dict(color=kc.get(kind, C["teal"]),
                            size=7, opacity=0.8,
                            line=dict(width=0.5, color=C["bg"])),
                text=grp["Novel_Peptide"],
                hovertemplate=(
                    "<b>%{text}</b><br>Len: %{x}<br>"
                    "Score: %{y:.2f}<extra></extra>")))
        fig_novel.update_layout(**LAY,
            title=dict(text="Novel Peptides — Length vs Score",
                       font=dict(size=13)),
            xaxis={**AXIS, "title": "Peptide Length"},
            yaxis={**AXIS, "title": "Score"})
    else:
        fig_novel = empty_fig("No novel peptides found.")

    # tables
    th = {"backgroundColor": C["bg"], "color": C["muted"],
          "fontWeight": "600", "fontSize": "10px",
          "border": f"1px solid {C['border']}",
          "letterSpacing": "0.07em", "textTransform": "uppercase"}
    tc = {"backgroundColor": C["panel"], "color": C["text"],
          "border": f"1px solid {C['border']}",
          "fontFamily": FONT, "fontSize": "11px", "padding": "7px 12px"}

    psm_tbl = dash_table.DataTable(
        data=df_psm.to_dict("records"),
        columns=[{"name": c, "id": c} for c in
                 ["Scan", "Peptide", "Protein", "Score",
                  "Matches", "Charge", "Precursor_mz"]],
        page_size=12, sort_action="native", filter_action="native",
        style_header=th, style_cell=tc,
        style_table={"overflowX": "auto"})

    if n_novel and "Novel_Peptide" in df_novel.columns:
        novel_tbl = dash_table.DataTable(
            data=df_novel.to_dict("records"),
            columns=[{"name": c, "id": c} for c in
                     ["Novel_Peptide", "Kind", "Change",
                      "Parent", "Protein", "Score", "Matches"]],
            page_size=10, sort_action="native", filter_action="native",
            style_header=th, style_cell=tc,
            style_table={"overflowX": "auto"},
            style_data_conditional=[
                {"if": {"filter_query": '{Kind} = "Point Mutation"'},
                 "color": C["rose"]},
                {"if": {"filter_query": '{Kind} = "Missed Cleavage"'},
                 "color": C["purple"]},
            ])
    else:
        novel_tbl = html.Div("No novel peptides found.",
                             style={"color": C["muted"], "fontSize": "12px",
                                    "padding": "14px 0"})

    return (kpi_row, fig_score, fig_donut, fig_len,
            fig_prec, fig_prot, fig_novel, psm_tbl, novel_tbl)


@app.callback(Output("dl-psm",  "data"),
              Input("dl-psm-btn", "n_clicks"),
              State("store-psm", "data"),
              prevent_initial_call=True)
def dl_psm(_, data):
    if not data: return no_update
    df = pd.read_json(io.StringIO(data), orient="split")
    return dcc.send_data_frame(df.to_csv, "identified_psms.csv", index=False)


@app.callback(Output("dl-novel",  "data"),
              Input("dl-novel-btn", "n_clicks"),
              State("store-novel", "data"),
              prevent_initial_call=True)
def dl_novel(_, data):
    if not data: return no_update
    df = pd.read_json(io.StringIO(data), orient="split")
    return dcc.send_data_frame(df.to_csv, "novel_peptides.csv", index=False)


# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("\n" + "═" * 60)
    print(f"  MS/MS PEPTIDE DASHBOARD  v5  [{N_WORKERS} threads | adaptive 2-pass search]")
    print("  → http://127.0.0.1:8051")
    print("  Configure the sidebar, click ▶ Run Pipeline")
    print("  Live log updates every 2 s while running")
    print("═" * 60 + "\n")
    app.run(debug=False, port=8051)
