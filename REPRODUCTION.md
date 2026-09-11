# Reproduction routes (v4.1; scope finalized 11 September 2026)

Use Python 3.12 and install `requirements.txt`. Outputs and raw sources should
be placed outside the release tree. All routes preserve the original 25 July
2026 analysis freeze. A database update is a new analysis version, never a
replacement for a missing historical input. This package is finalized with
the existing verified primary analysis and extensions. The supplementary
AmtDB cross-database audit remains historical and inherited. No AmtDB update
or new three-source run is claimed, and requesting a file from database
authors is not required for this verified scope. See `AS_IS_STATUS.md` for
the final scope decision and publication requirements.

| Route | Starting point | What is regenerated | What remains inherited |
| --- | --- | --- | --- |
| `derived` | Included 489-row project-coded CSV | Main statistical tables/figures 1–5, global sensitivities and v4 extensions | Category definitions, original catalogue metadata and database-coverage audit |
| `aadr` | Exact AADR v66.p1 annotation | Source parsing, marker-specific deduplication, locality normalization, project-coded input, then the same statistical analyses | AmtDB/aYChr cross-database audit and Figure 6; definition/metadata initialization is explicitly logged |
| `raw` (currently blocked) | Exact AADR, AmtDB v1.009 and aYChr files | Optional future three-source pipeline, source catalogues, database audit, main/global analyses, rebuilt input and v4 extensions | Requires the missing original AmtDB CSV; literature/source-scope review and separate calibration simulation are not rerun |

## Verified on 8 September 2026

- Live downloads matched the frozen AADR annotation (13,450,350 bytes) and
  aYChr workbook (272,253 bytes).
- Independent AADR extraction rebuilt the 489-row analytical input byte for
  byte, with 136 sites and 438 mtDNA/229 Y calls.
- A full isolated `derived` run passed every output gate and reproduced all
  32 compared CSV files byte for byte: 20 main tables and 12 extension tables.
  The inherited database-coverage table is excluded from that count.
- All 39 local tests passed with the verified AADR fixture enabled.

The primary and extension calculations are therefore verified from the
included input, with an independent exact-source reconstruction of that same
input. Both regenerated Figure 5 PNGs and the v4 provenance were published
on 10 September 2026.

The complete three-source route has not passed because the exact AmtDB
v1.009 CSV is missing. Independent AADR extraction and the full statistical
rerun do not reproduce the separate cross-database audit. Its saved coverage
table, matching summaries and Figure 6 are retained as historical inherited
results, outside the completed primary-analysis verification. Recorded
validation is in `validation/reproduction_v4_1.json`.

## Download or import frozen sources

    python analysis/fetch_sources.py --resource aadr --resource aychr \
      --output-root ../central-asia-sources --timeout 60 --retries 2 \
      --report ../source_acquisition.json

The downloader uses `data/SOURCES.tsv`. The AADR datafile ID and aYChr commit
are fixed. SHA-256 must match before atomic installation; incomplete bytes
are discarded. A preexisting file with different bytes is left unchanged and
reported as a conflict. Only transient network errors are retried. Selecting
two ready sources does not report that all three sources are ready.

AmtDB v1.009 has no verified historical export URL. Importing it is needed
only for the optional full three-source route; no author request is part of
finishing the current primary-analysis package. If the original becomes
available:

    python analysis/fetch_sources.py --resource amtdb \
      --output-root ../central-asia-sources \
      --local amtdb=/absolute/path/to/original-v1.009.csv \
      --report ../amtdb_acquisition.json

Required SHA-256:
`531e8ee8fae181124f5a9b77b6fe8d677e64e35b815be2a3965020244fe31057`.
The live v1.010 export is not a substitute. No network message is sent to the
database maintainers by this tool. Missing AmtDB is reported as `BLOCKED` and
the all-source command exits nonzero.

## Reproduce the main study from AADR

    python analysis/reproduce.py --mode aadr --fetch \
      --source-root ../central-asia-sources \
      --output ../aadr-reproduction

The AADR extractor verifies the source before parsing. It rebuilds the full
501-person catalogue and the 489-person primary selection, preserving the
136 normalized sites, 438 mtDNA calls and 229 Y calls. The coded input must
match SHA-256
`692a69cf38cc736f96ea5aa6b3b15024a9a49d50c03ae3e1470a8dc475504cc3`.
Source-labelled work files are kept in the external run directory and are not
automatically uploaded anywhere.

The required AADR source SHA-256 is
`98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c`.
The extractor saves the primary source-labelled catalogue and its extraction
manifest. Its working directory must be absent or empty, and its final output
must be a new external file.

## Optional future three-source reproduction

    python analysis/reproduce.py --mode raw --fetch \
      --source-root ../central-asia-sources \
      --amtdb /absolute/path/to/original-v1.009.csv \
      --output ../three-source-reproduction

If sources are already verified in their expected paths, omit `--fetch` and
`--amtdb`. Raw analysis verifies all three hashes again before computation.
The runner does not claim success when acquisition, analysis, output checks,
input reconstruction or reference comparison fails.

## Offline statistical reproduction

    python analysis/reproduce.py --mode derived \
      --output ../derived-reproduction

Only three initialization files are inherited: `results_summary.json`,
`analysis_manifest.json`, and `tables/database_coverage_by_country.csv`. Their
input hashes are recorded in the run report. The first two supply fixed
category definitions and original source metadata, including information not
recoverable from the 489-row input. Statistical summary sections are then
updated by the actual rerun. No old statistical CSV or figure is seeded.

The run report distinguishes inherited content from generated files and
records the complete CSV comparison, floating-point tolerance, environment,
stage logs and exit status. Different floating-point formatting can produce
different file hashes while the numeric comparison still passes; the report
records both byte identity and the maximum numeric difference.

## Tests and CI

    PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v

The integrity workflow discovers every test module, including acquisition,
runner and raw-extraction regressions. Download failures and interruptions
are simulated without network access. One fixture-based test in that same
suite verifies actual AADR extraction; it is skipped unless the source
location is supplied:

    ADNA_AADR_FIXTURE=../central-asia-sources/data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno \
      PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v

That source must pass its frozen hash check. CI does not download upstream
databases or run the full resampling analysis on every pull request. Its
remaining steps validate the saved aggregate/extension outputs and release
checksums; a passing CI run is not a three-source rerun.

## Scope and public release

The runner uses all frozen resampling settings; it has no reduced-replicate
mode that could be mistaken for a release verification. The extension script
itself supports explicit alternate settings for diagnostics, and records them
in its own manifest.

The saved primary inference excludes singleton site-period profiles. Neither
Y encoding supports faster paternal turnover. This update changes acquisition,
reproduction, dependency declarations and failure handling; it does not add
new individuals or new demographic claims.

Project-coded records and site-profile vectors remain potentially linkable
to the public source. The v4.1 candidate preserves the same analytical input
and profile material already present in the public draft PR. Raw catalogues
remain in external working directories. Author approval, archaeological
interpretation review and archival release metadata remain pending. Recovering
the AmtDB snapshot is an optional extension of reproducibility scope; it is
mandatory only before claiming a fresh complete three-source reproduction.
This scope clarification does not merge the draft or create a release.
See `PUBLIC_RELEASE_CHECKLIST.md`.
