# 🧬 MS/MS Peptide Identification Pipeline

> An interactive, browser-based proteomics pipeline for peptide identification from raw MS/MS data — with adaptive two-pass search, novel peptide discovery, and a full Dash dashboard.

![Python](https://img.shields.io/badge/python-3.9%2B-blue?style=flat-square)
![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)
![Dashboard](https://img.shields.io/badge/dashboard-Dash%2FPlotly-orange?style=flat-square)
![Status](https://img.shields.io/badge/status-active-brightgreen?style=flat-square)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Dashboard Guide](#dashboard-guide)
- [Input Formats](#input-formats)
- [Output Files](#output-files)
- [Performance](#performance)
- [Configuration Reference](#configuration-reference)
- [Example Dataset](#example-dataset)
- [Benchmarking](#benchmarking)
- [Contributing](#contributing)
- [License](#license)

---

## Overview

This tool provides an end-to-end MS/MS peptide identification workflow entirely in Python — no external search engines required. It runs as an interactive web dashboard where you configure all parameters, monitor the pipeline in real time, and explore results through interactive charts and tables.

**Designed for:**
- Proteomics researchers who want a transparent, hackable search engine
- Bioinformaticians benchmarking peptide identification algorithms
- Anyone exploring novel/variant peptides not represented in standard databases

---

## Features

| Feature | Description |
|---|---|
| **Adaptive Two-Pass Search** | 10% pilot run learns optimal tolerances; 90% searched with tight parameters |
| **Novel Peptide Discovery** | Point-mutation rescue and missed-cleavage joins on unmatched spectra |
| **Interactive Dashboard** | All parameters set in-browser; live log console; 7 result panels |
| **1% FDR Filtering** | Target-decoy competition with reversed sequences |
| **Multi-threaded Scoring** | ThreadPoolExecutor uses all available CPU cores |
| **Numpy Vectorised Scoring** | b/y-ion hyperscore via broadcasting — ~6× faster than pure Python |
| **Bisect Mass Index** | O(log n) precursor lookup — ~15× faster than float-dict iteration |
| **LRU Ion Cache** | Theoretical ions computed once per peptide, reused across all spectra |
| **Auto Download** | Fetches mzML + FASTA from URL automatically (with progress) |
| **CSV Export** | One-click download of PSMs and novel peptides |

---

## Installation

### Prerequisites

- Python 3.9 or higher
- pip

### Install dependencies

```bash
pip install pyteomics biopython requests tqdm dash plotly pandas numpy
```

Or using the provided requirements file:

```bash
pip install -r requirements.txt
```

### Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/msms-peptide-pipeline.git
cd msms-peptide-pipeline
```

---

## Quick Start

```bash
python msms_peptide_dashboard.py
```

Then open **http://127.0.0.1:8051** in your browser.

The dashboard pre-fills with a real PRIDE dataset (PXD000001, *Erwinia carotovora* TMT). Just click **▶ Run Pipeline** to see it in action immediately.

---

## How It Works

```
┌─────────────────────────────────────────────────────────────────┐
│                    MS/MS PIPELINE FLOW                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. DOWNLOAD       mzML + FASTA (URL or upload)                 │
│         │                                                       │
│  2. DIGEST         Tryptic in-silico digest → mass index        │
│         │          (sorted array + bisect, LRU ion cache)       │
│         │                                                       │
│  3. ADAPTIVE SEARCH                                             │
│         │                                                       │
│         ├── 10% spectra ──► PASS 1: Broad search               │
│         │                   Learn: optimal ppm + frag tol       │
│         │                                                       │
│         └── 90% spectra ──► PASS 2: Narrow search              │
│                              (learned tight tolerances)         │
│                              ThreadPoolExecutor parallel        │
│                              Numpy vectorised b/y scoring       │
│                                                                 │
│  4. FDR FILTER     Target-decoy, 1% FDR (configurable)         │
│         │                                                       │
│  5. NOVEL DISCOVERY  Point mutations + missed cleavage joins    │
│         │                                                       │
│  6. DASHBOARD      7 interactive charts + searchable tables     │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Adaptive Two-Pass Search

The key innovation in this pipeline is the **adaptive parameter learning**:

1. **Pass 1 (10% of spectra, broad tolerances):** A random representative sample is searched with wide tolerances (e.g. 20 ppm precursor, 50 mDa fragment). High-confidence PSMs are collected.

2. **Parameter learning:** The 95th percentile of observed precursor mass errors and fragment ion errors is computed from Pass 1 hits. These are multiplied by a 1.5× safety margin to define tight narrow windows.

3. **Pass 2 (90% of spectra, narrow tolerances):** The remaining spectra are searched with the learned tight tolerances, resulting in far fewer candidates per spectrum and much faster scoring.

Typical result: a 20 ppm broad search might yield 50 candidates/spectrum. After learning, a 4–6 ppm narrow search yields 5–10 — a **5–10× reduction** in scoring work on 90% of the data.

### Novel Peptide Discovery

After standard database search, unmatched spectra are inspected for:

- **Point mutations:** Every single amino acid substitution of nearby database peptides is tried. Accepted if mass matches within tolerance AND score ≥ 70% of parent peptide score.
- **Missed cleavage joins:** Adjacent tryptic fragments are joined and scored if the combined mass matches the precursor.

---

## Dashboard Guide

| Panel | Description |
|---|---|
| **KPI Cards** | Total spectra, PSMs, proteins, novel peptides, ID rate |
| **PSM Score Distribution** | Hyperscore histogram — confident IDs cluster at high values |
| **ID Rate Donut** | % spectra identified at selected FDR |
| **Peptide Length** | Distribution confirms expected tryptic digestion pattern |
| **Precursor m/z vs Score** | Coloured by charge state — reveals instrument behaviour |
| **Top 20 Proteins** | PSM count per protein |
| **Novel Peptides Scatter** | Point mutations (red) vs missed cleavage (purple) |
| **PSM Table** | Sortable, filterable, paginated — all identified peptides |
| **Novel Table** | Variant peptides with change annotation (e.g. `K7R`) |

### Sidebar Controls

| Control | Description |
|---|---|
| mzML Input | URL (PRIDE FTP) or file upload |
| FASTA Input | URL (UniProt REST) or file upload |
| Enzyme | Trypsin, Trypsin/P, Lys-C, Chymotrypsin |
| Missed Cleavages | 0–3 |
| Peptide Length | Min / Max AA length |
| Precursor ppm | Initial broad tolerance for pilot search |
| Fragment Da | Initial broad fragment tolerance |
| Min Fragment Matches | Minimum b/y ions matched to accept a PSM |
| FDR Threshold | 0.1%, 1%, or 5% |
| Novel Discovery | Toggle point mutations / missed cleavage joins |
| Rescue Score Ratio | Novel peptide must score ≥ this fraction of parent |

---

## Input Formats

### mzML

Standard mzML format (v1.1.0+). Generated by:
- ProteoWizard `msconvert`
- Thermo Scientific software
- Most modern mass spec data systems

### FASTA

Standard FASTA protein sequence format. Sources:
- [UniProt](https://www.uniprot.org/downloads) reviewed Swiss-Prot entries
- [PRIDE FTP](https://ftp.pride.ebi.ac.uk/) experiment-specific databases
- Any custom FASTA file

---

## Output Files

All outputs are saved to `./pipeline_output/`:

| File | Description |
|---|---|
| `identified_psms.csv` | All PSMs passing FDR — scan, peptide, protein, score, charge |
| `novel_peptides.csv` | Novel/variant peptides with change annotation |
| `learned_params.txt` | Tolerances learned from pilot search |
| `spectra.mzML` | Cached downloaded mzML (reused on subsequent runs) |
| `database.fasta` | Cached downloaded FASTA |

### identified_psms.csv columns

| Column | Description |
|---|---|
| Scan | Spectrum scan ID |
| Peptide | Identified peptide sequence |
| Protein | Protein accession |
| Score | Hyperscore |
| Matches | Number of b/y ions matched |
| Charge | Precursor charge state |
| Precursor_mz | Precursor m/z |
| PepLen | Peptide length (AA) |

### novel_peptides.csv columns

| Column | Description |
|---|---|
| Novel_Peptide | Variant peptide sequence |
| Kind | `Point Mutation` or `Missed Cleavage` |
| Change | Mutation notation e.g. `K7R` (Lys→Arg at position 7) |
| Parent | Parent database peptide |
| Protein | Source protein accession |
| Score | Hyperscore |
| Matches | Fragment ions matched |

---

## Performance

Benchmarked on a 6-core machine with the PRIDE PXD000001 dataset (~7,500 MS2 spectra):

| Version | Strategy | Time |
|---|---|---|
| v3 (baseline) | Single-threaded, Python loops | ~45 min |
| v4 | Multi-threaded, numpy scoring, bisect index | ~8 min |
| v5 (current) | v4 + adaptive two-pass search | ~3 min |

### Speed optimisations summary

| Optimisation | Technique | Speedup |
|---|---|---|
| b/y-ion scoring | NumPy broadcasting instead of nested loops | ~6× |
| Precursor lookup | Sorted array + `bisect` instead of float-dict iteration | ~15× |
| Ion generation | `@lru_cache` per peptide sequence | eliminates redundancy |
| Spectrum scoring | `ThreadPoolExecutor` across all CPU cores | ~N_CPU× |
| Search scope | Adaptive narrow tolerances on 90% of data | ~5–10× |

---

## Configuration Reference

Key constants at the top of `msms_peptide_dashboard.py`:

```python
N_WORKERS = max(1, os.cpu_count() - 1)   # threads for scoring

# Pre-filled default URLs (PRIDE PXD000001)
MZML_DEFAULT  = "https://ftp.pride.ebi.ac.uk/..."
FASTA_DEFAULT = "https://ftp.pride.ebi.ac.uk/..."
```

All search parameters are configurable from the dashboard sidebar — no need to edit source code for typical use.

---

## Example Dataset

The dashboard pre-fills with **PRIDE PXD000001** — a well-characterised *Erwinia carotovora* TMT proteome dataset:

| Property | Value |
|---|---|
| Organism | *Erwinia carotovora* (plant pathogen) |
| Instrument | Orbitrap (HCD fragmentation) |
| Labelling | TMT 10-plex |
| MS2 spectra | ~7,500 |
| FASTA | 4,204 reviewed proteins |
| PRIDE accession | PXD000001 |

**Direct links:**
```
mzML:  https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML
FASTA: https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/erwinia_carotovora.fasta
```

---

## Benchmarking

To benchmark against known results, use the **iPRG2012** dataset which has a published answer key of correct PSMs:

```
mzML:  https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/02/PXD000612/
FASTA: https://www.uniprot.org/...
```

Metrics to evaluate:

- **Sensitivity** — % of known true PSMs recovered
- **FDR accuracy** — actual vs nominal FDR at 1% threshold
- **Search speed** — spectra scored per second
- **Novel recovery** — known variants identified without database entry

Run the pipeline on the benchmark dataset and compare `identified_psms.csv` against the published answer key using the companion script:

```bash
python benchmark/compare_results.py \
    --our-psms pipeline_output/identified_psms.csv \
    --answer-key benchmark/iprg2012_answers.csv
```

---

## Project Structure

```
msms-peptide-pipeline/
├── msms_peptide_dashboard.py   # Main application (pipeline + dashboard)
├── requirements.txt            # Python dependencies
├── .gitignore                  # Excludes pipeline_output/, caches
├── LICENSE                     # MIT License
├── README.md                   # This file
├── CHANGELOG.md                # Version history
└── benchmark/
    └── compare_results.py      # Benchmarking utility
```

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/faster-scoring`)
3. Make your changes with clear commit messages
4. Add or update tests if applicable
5. Open a Pull Request with a description of what changed and why

### Ideas for contribution

- Add support for modifications (phosphorylation, oxidation, carbamidomethylation)
- Implement an E-value or PEP (posterior error probability) score
- Add MGF input format support alongside mzML
- Integrate with MSFragger or Comet for comparison benchmarking
- Add a mass accuracy calibration plot

---

## License

MIT License — see [LICENSE](LICENSE) for full text.

---

## Citation

If you use this tool in your research, please cite:

```
MS/MS Peptide Identification Pipeline (2025)
https://github.com/YOUR_USERNAME/msms-peptide-pipeline
```

---

## Acknowledgements

- [Pyteomics](https://pyteomics.readthedocs.io/) — mzML parsing and mass calculation
- [BioPython](https://biopython.org/) — FASTA parsing
- [Dash / Plotly](https://dash.plotly.com/) — interactive dashboard
- [PRIDE Archive](https://www.ebi.ac.uk/pride/) — public proteomics data
