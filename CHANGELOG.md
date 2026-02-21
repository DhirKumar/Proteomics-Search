# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [v5.0.0] — 2025

### Added
- **Adaptive two-pass search** — 10% pilot run learns optimal precursor ppm
  and fragment Da tolerances; 90% of spectra searched with tight learned
  parameters for 5–10× reduction in candidate scoring work
- `_learn_params()` — statistical parameter inference from pilot PSMs
  (95th percentile of observed errors × 1.5 safety margin)
- `_read_all_spectra()` — single-pass mzML reader, separates I/O from scoring
- `_parallel_score()` — reusable parallel scoring helper used by both passes
- `learned_params.txt` output file — records pilot statistics and learned
  tolerances for reproducibility

### Changed
- `run_search()` refactored into three composable helper functions
- Pilot sample is randomly shuffled (seed=42) to ensure representativeness
- Progress reporting now shows pass number and percentage separately

---

## [v4.0.0] — 2025

### Added
- **NumPy vectorised b/y-ion scoring** via broadcasting — ~6× faster than
  nested Python loops
- **`MassIndex` class** — sorted array with `bisect` range lookup, replacing
  float-keyed `defaultdict` + while-loop (~15× faster precursor lookup)
- **`@lru_cache`** on `theoretical_ions_cached()` — ions computed once per
  unique peptide sequence, reused across all spectra
- **Ion cache pre-warming** during `build_index` — zero cold-start cost
  when search begins
- **`ThreadPoolExecutor`** parallel spectrum scoring using all CPU cores

### Changed
- `theoretical_ions()` renamed to `theoretical_ions_cached()` and decorated
  with `@lru_cache(maxsize=131072)`
- `lookup()` function replaced by `MassIndex.lookup()` method
- `run_search()` now reads all spectra first (sequential I/O), then scores
  in parallel (CPU-bound)

---

## [v3.0.0] — 2025

### Added
- **Background threading** — pipeline runs in a `daemon` thread so Dash
  stays fully responsive during long searches
- **Live log console** — `dcc.Interval` polls job state every 2 seconds;
  log lines stream into the dashboard in real time
- **`JobManager` class** — thread-safe job state tracking with UUID job IDs
- **Polling auto-stop** — interval disables itself on job completion or error
- Corrected PRIDE URLs (`PXD000001`, `2012/03/` path)
- Separate FASTA from PXD000001 (`erwinia_carotovora.fasta`) instead of
  UniProt REST stream

### Fixed
- Non-responsive Run Pipeline button (was blocking Dash's main thread)
- Wrong PRIDE accession in default URLs (PXD000561 → PXD000001)

---

## [v2.0.0] — 2025

### Added
- Full interactive **Dash dashboard** with 7 result panels
- All pipeline parameters collected in dashboard sidebar (no config editing)
- mzML and FASTA input via URL or drag-and-drop file upload
- KPI stat cards (spectra, PSMs, proteins, novel peptides, ID rate)
- PSM score distribution histogram
- Identification rate donut chart
- Peptide length distribution bar chart
- Precursor m/z vs score scatter (coloured by charge state)
- Top 20 proteins by PSM count horizontal bar chart
- Novel peptides length vs score scatter
- Searchable, sortable PSM table and novel peptides table
- One-click CSV download for PSMs and novel peptides

---

## [v1.0.0] — 2025

### Added
- Initial release — command-line pipeline
- In-silico tryptic digest with configurable missed cleavages
- Float-keyed mass index for precursor lookup
- b/y-ion hyperscore with Python loops
- Target-decoy FDR filtering (reversed sequences)
- Novel peptide discovery: point mutations and missed-cleavage joins
- Download from PRIDE FTP and UniProt REST
- CSV outputs: `identified_psms.csv`, `novel_peptides.csv`
