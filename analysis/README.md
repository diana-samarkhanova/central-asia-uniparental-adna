# Reproducing the Central Asia uniparental aDNA analysis

This directory contains the frozen scientific workflow for the evidence release
dated 25 July 2026. The GitHub repository is aggregate-only: it includes the
analysis code, figures and aggregate statistical outputs, but no individual
catalogues, exact-ID crosswalks or singleton-revealing site-lineage profiles.
The audited statistical revision and corrected derived outputs are dated
15 August 2026.

Exact reproduction starts from the frozen AADR, AmtDB and aYChr-DB metadata
files. The workflow extracts, deduplicates and harmonizes person-level working
data locally, then produces the aggregate release. It does not reprocess FASTQ,
BAM or other sequencing reads. The estimand is the composition of published
individuals and equal-weighted published localities; it is not a past
population-frequency estimate.

## Environment

Use Python 3.12.13 with the pinned packages in `requirements.txt`:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Alternatively, create the Conda environment from `environment.yml`.

## Verify the committed release

The checksum file uses the standard GNU two-column format. From the repository
root, verify it and run the focused scientific and public-release checks:

```bash
sha256sum -c SHA256SUMS.txt
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
PYTHONDONTWRITEBYTECODE=1 python analysis/test_aggregate_release.py results/aadr-v66p1_2026-07-25
python analysis/verify_release.py
```

`verify_release.py` fails if a checksum is missing or stale; if raw databases,
office/PDF documents, archives or cache files are present; if forbidden
individual/site-lineage tables or transient absolute paths remain; if code
hashes disagree with the manifests; or if the frozen counts and
50,000-replicate paired bootstrap are inconsistent.

Maintainers should generate checksums only after reviewing the intended release
tree:

```bash
python analysis/generate_checksums.py
sha256sum -c SHA256SUMS.txt
python analysis/verify_release.py
```

The generator excludes `SHA256SUMS.txt` itself, `.git` and cache files. It does
not make an unsafe file acceptable: the independent verifier still rejects
raw, document, archive, temporary and cache artifacts.

## Reproduce from the three frozen source databases (blocked pending two locators)

Follow `data/README.md` and `data/SOURCES.tsv` to obtain the exact files. Do not
commit them. The hashes below remain authoritative, but the exact AADR
Dataverse datafile locator and an immutable public AmtDB v1.009 archive URL
have not been independently verified. The aYChr workbook is pinned to commit
`bc770a59ace8cd4c042c6f903d620d93ee751eb0`. Expected local paths and SHA-256
values are:

| Resource | Expected local path | SHA-256 |
|---|---|---|
| AADR v66.p1 | `data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno` | `98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c` |
| AmtDB v1.009 | `data/raw/amtdb_v1_009/amtdb_v1.009_metadata.csv` | `531e8ee8fae181124f5a9b77b6fe8d677e64e35b815be2a3965020244fe31057` |
| aYChr-DB v5 | `data/raw/aychr_db/a-YChr-DB_V5.xlsx` | `e297110a18cba73d4044e8a95c0fae98d7f48633ad6ccfae1cd364e460eb1b3c` |

Verify the inputs before doing any analysis:

```bash
python analysis/verify_inputs.py \
  --aadr data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno \
  --amtdb data/raw/amtdb_v1_009/amtdb_v1.009_metadata.csv \
  --aychr data/raw/aychr_db/a-YChr-DB_V5.xlsx
```

`run_analysis.py` repeats this verification and refuses files whose hashes do
not match the frozen registry. Run the complete extraction and analysis into a
new work directory:

```bash
MPLCONFIGDIR=.mplconfig python analysis/run_analysis.py \
  --aadr data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno \
  --amtdb data/raw/amtdb_v1_009/amtdb_v1.009_metadata.csv \
  --aychr data/raw/aychr_db/a-YChr-DB_V5.xlsx \
  --outdir work/full-reproduction \
  --bootstrap 2000 \
  --paired-bootstrap 50000 \
  --permutations 9999 \
  --callability-resamples 99999 \
  --date-draws 5000 \
  --seed 20260725

MPLCONFIGDIR=.mplconfig python analysis/run_global_sensitivities.py \
  --analysis-output work/full-reproduction \
  --permutations 1999 \
  --seed 20260726

python analysis/test_full_outputs.py work/full-reproduction
```

The frozen settings are 2,000 site-cluster bootstrap replicates, 50,000 paired
replicates per Y encoding, 9,999 primary resamples, 99,999 fixed-margin
callability simulations and 5,000 chronological-bin assignment scenarios. The resampling
implementation assigns one shared draw multiplicity to each `country +
locality` cluster across all periods; the paired analysis shares the same
cluster draw across mtDNA and Y. Named random streams make one procedure
independent of unrelated replicate counts.

Paired tables report the original-sample Δ as the point estimate and retain
the bootstrap median only as a diagnostic.
`bootstrap_two_sided_sign_tail_probability` measures ordinary-bootstrap sign
stability and is not a null-hypothesis *P* value. The ISOGG-prefix analysis is
a nomenclature sensitivity, not a phylogenetic re-call. Files retaining
`date_uncertainty` in their names contain assumption-based bin-assignment
scenarios shared across markers, not calibrated-date posterior draws.

The generated working directory contains person-level intermediates and must
not be committed. `recompute_from_catalogue.py`, `run_global_sensitivities.py`
and `sanitize_release.py` are retained for maintainers who have a locally
generated catalogue; their required person-level inputs are intentionally not
included in GitHub. The post-v66 evidence map and literature audit also include
manual curation steps and are not reconstructed by these commands.

## Input licensing and redistribution

- **AADR v66.p1:** the frozen source registry records CC0 1.0. Cite the AADR
  descriptor, Dataverse DOI and original component studies.
- **AmtDB v1.009:** the registry records CC BY 4.0. Cite the database article
  and preserve attribution.
- **aYChr-DB v5:** the article/supplement is recorded as CC BY 4.0 and the
  workbook is pinned to commit
  `bc770a59ace8cd4c042c6f903d620d93ee751eb0`; the repository did not display a
  separate reuse license during the freeze, so confirm reuse terms before
  version 1.0.

The upstream database files are intentionally excluded from Git. Original code
is MIT-licensed. Documentation, figures and curated output tables are intended
for CC BY 4.0 release after coauthor approval. See `THIRD_PARTY_NOTICES.md`,
`LICENSES/` and `PUBLIC_RELEASE_CHECKLIST.md` for the controlling release notes.

## Coordinate masking and restricted records

Exact archaeological coordinates are not used in the statistical models.
`run_analysis.py` writes local person-level working catalogues with latitude and
longitude rounded to one decimal degree. For a previously generated local
result directory, apply the same policy with:

```bash
MPLCONFIGDIR=.mplconfig python analysis/sanitize_release.py \
  --analysis-output work/full-reproduction \
  --digits 1
```

The sanitizer also regenerates the sampling map and records the masking policy
in `analysis_manifest.json`. Those catalogues still remain outside GitHub.
`verify_release.py` instead enforces the aggregate-only table contract.
Person-level rows from the restricted post-freeze source must remain excluded;
only citation-level summaries may be released.

## Script responsibilities

- `verify_inputs.py`: check frozen upstream hashes before reading data.
- `run_analysis.py`: extract, deduplicate, harmonize, analyze and mask public
  coordinates.
- `recompute_from_catalogue.py`: reproduce statistical results from a locally
  generated person-level analytical catalogue; the input is not committed.
- `run_global_sensitivities.py`: marker-wide sensitivity analyses.
- `sanitize_release.py`: mask an existing output directory and regenerate its
  sampling map.
- `test_aggregate_release.py`: aggregate-only release regression checks.
- `test_full_outputs.py`: full/granular scientific regression checks.
- `test_outputs.py`: backward-compatible alias for the aggregate check; new
  documentation uses the explicit filenames above.
- `generate_checksums.py`: write deterministic GNU-compatible checksums.
- `verify_release.py`: enforce the public-release contract.
