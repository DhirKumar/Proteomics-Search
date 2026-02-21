# Contributing to MS/MS Peptide Identification Pipeline

Thank you for your interest in contributing! This document explains how to
get started, what we're looking for, and how to submit changes.

---

## Getting Started

### 1. Fork and clone

```bash
git clone https://github.com/YOUR_USERNAME/msms-peptide-pipeline.git
cd msms-peptide-pipeline
```

### 2. Create a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate          # Linux / macOS
.venv\Scripts\activate             # Windows
pip install -r requirements.txt
```

### 3. Create a branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/issue-description
```

---

## Code Style

- Follow **PEP 8** — 4-space indentation, max 88 characters per line
- Use **type hints** for all new functions
- Add a **docstring** to every public function
- Keep functions focused — one function, one responsibility

---

## Testing Your Changes

Before submitting, verify:

```bash
# 1. Syntax check
python -m py_compile msms_peptide_dashboard.py

# 2. Run the dashboard and confirm it starts
python msms_peptide_dashboard.py
# open http://127.0.0.1:8051 and click ▶ Run Pipeline

# 3. Benchmark script (optional)
python benchmark/compare_results.py --help
```

---

## What We Welcome

### High priority
- **Modification support** — phosphorylation, oxidation, carbamidomethylation
- **MGF input format** — alongside existing mzML support
- **PEP / E-value scoring** — posterior error probability as alternative to FDR
- **Mass accuracy calibration** — recalibrate precursor masses before search
- **Percolator-style rescoring** — machine learning PSM re-ranking

### Good first issues
- Add unit tests for `score_psm`, `MassIndex.lookup`, `_learn_params`
- Add a `--headless` CLI mode that runs without launching the dashboard
- Support gzipped mzML files (`.mzML.gz`)
- Add charge-state deconvolution for uncharacterised precursors

### Infrastructure
- Docker / conda environment file
- GitHub Actions CI (syntax check + import test on push)
- Sphinx documentation

---

## Pull Request Guidelines

1. **One PR per feature or fix** — keep PRs focused and reviewable
2. **Describe what and why** in the PR description, not just what changed
3. **Update CHANGELOG.md** under an `[Unreleased]` section
4. **Update README.md** if you add new parameters, outputs, or behaviour
5. Reference any related issues with `Fixes #123` or `Relates to #456`

---

## Reporting Bugs

Please open a GitHub issue with:

- Python version (`python --version`)
- OS (Windows / macOS / Linux)
- Full error traceback from the dashboard log console or terminal
- The mzML and FASTA sources you were using (URLs or file description)
- Steps to reproduce

---

## Questions?

Open a GitHub Discussion or issue tagged `question`. We're happy to help.
