# Central Asia uniparental ancient DNA

Reproducible secondary analysis of published ancient mitochondrial and
Y-chromosome assignments from Kazakhstan, Kyrgyzstan, Tajikistan,
Turkmenistan and Uzbekistan.

Status: public repository with an unmerged pre-submission candidate. Evidence
is frozen on 25 July 2026; statistical/literature/site-normalization corrections
are synchronized to 21 August 2026. Reproduction tooling was updated on
8 September 2026 (v4.1 candidate); its reproducibility scope was finalized
on 11 September 2026 using the existing verified inputs and results. Both
regenerated Figure 5 PNGs and their v4 provenance were published on
10 September 2026. This is not yet a citable public release. See
`AS_IS_STATUS.md` for the completed reproducibility scope and remaining
publication steps.

## Scope

The project harmonizes AADR v66.p1 metadata and describes the composition of
published archaeological individuals and equal-weighted published localities.
The primary analysis and statistical extensions are reproducible from the
included input, independently rebuilt from the exact AADR source. The
marker-specific database audit is retained as a historical supplementary
result; its AmtDB component has not been freshly reproduced. The estimand is
the published evidence base, not the population frequency of lineages in
ancient Central Asia.

The frozen AADR catalogue contains 501 unique archaeological individuals.
The primary interval, 3500 BCE to 1500 CE, contains 489 individuals from
136 normalized country–locality sites, with 438 mitochondrial and 229
Y-chromosome calls. The
country-adjusted all-profile models detect period-associated structure, but
repeated-locality and sensitivity analyses show that it cannot be interpreted
as a region-wide demographic replacement. In 216 paired men from 98 sites,
the broad-L1 contrast Δ(Y−mtDNA) is −0.106381 (95% bootstrap interval −0.267291
to 0.002024). Because that percentile interval crosses zero, the package does
not assign a hypothesis-test *P* value to this contrast. A 12-category AADR
ISOGG-prefix sensitivity gives
Δ=−0.011292 (95% interval −0.154022 to 0.075681); the paired encoding-induced
change from broad L1 is +0.095089 (95% interval 0.023370 to 0.176052). Neither
encoding supports faster Y-chromosome turnover.

The revised primary site-profile estimand requires at least two marker calls
per profile and uses 9,999 null-imposed HC2 cluster-wild resamples. The mtDNA
association is not supported (71 profiles; Holm *P*=0.2224), whereas the Y
association persists (48 analyzed profiles; Holm *P*=0.0140). The >=3 results
are stringent, low-information sensitivities. Figure 5 v4 visually
distinguishes the 9,999-resample revised profile-size rows from 1,999-resample
exploratory checks; exact revised results are in
`analysis_extensions/statistical_extensions_v4/results/stable_profile_cluster_tests.csv`.

## Repository map

- analysis: scientific analysis, sensitivity analysis, input validation and
  integrity checks.
- analysis_extensions/statistical_extensions_v4: revised >=2 primary-profile
  tests (9,999 resamples), >=3 sensitivity, finite-sample TV diagnostics,
  equal-width/time-grid checks and their regression test.
- audits/source_scope_v4: immutable-source resolution, reversible locality-
  string audit and structured Russian-language evidence search.
- data: acquisition instructions, frozen-source registry, and a schema-locked
  AADR-derived analytical input with project-coded records/sites/studies,
  broad marker categories and one-decimal coordinates. Source identifiers,
  locality strings, skeletal fields and terminal haplogroup calls are omitted.
- results/aadr-v66p1_2026-07-25: versioned figures, machine-readable summaries
  and aggregate derived tables. Individual-level catalogues, exact-ID
  crosswalks and singleton-revealing site-lineage profiles are not committed.
- results/post-v66-audit: citation-level audit of studies not fully represented
  in the frozen AADR release.
- literature_audit: targeted-search log and verified bibliography. The search
  is not a completed PRISMA systematic/scoping review.
- tests: focused tests for harmonization, resampling, source acquisition,
  verified extraction, provenance, and preservation of existing outputs.

The unpublished manuscript is deliberately excluded until authorship,
affiliations and coauthor approval are confirmed.

## Quick integrity check

Create a Python 3.12 environment and install the frozen dependencies:

    python -m pip install -r requirements.txt

Then run:

    PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
    PYTHONDONTWRITEBYTECODE=1 python analysis/test_aggregate_release.py results/aadr-v66p1_2026-07-25
    PYTHONDONTWRITEBYTECODE=1 python analysis_extensions/statistical_extensions_v4/test_results.py
    PYTHONDONTWRITEBYTECODE=1 python analysis/verify_release.py

These commands validate the project-coded input and aggregate release,
including the frozen counts, statistical summaries, code hashes and
50,000-replicate paired bootstrap.

The 8 September 2026 local validation passed all 39 tests with the verified
AADR fixture enabled. Independent AADR extraction reproduced the frozen
489-row input byte for byte, and a complete isolated `derived` run reproduced
all 32 compared CSVs byte for byte: 20 main tables and 12 extension tables.
Live AADR and aYChr downloads also matched their frozen hashes. These checks
complete the verified primary-analysis and extension scope of this package.
The missing original AmtDB v1.009 CSV limits fresh reproduction of the
separate cross-database audit; its coverage table and Figure 6 remain
historical, inherited outputs. No new AmtDB version or completed three-source
rerun is claimed. See `validation/reproduction_v4_1.json` for the recorded
validation and `REPRODUCTION.md` for the optional AADR fixture test.

## Reproduce in an isolated output directory

The recommended command downloads the exact AADR annotation, independently
rebuilds the 489-row analytical input, requires a byte-for-byte match to the
frozen input, and then recomputes the main statistics and v4 extensions:

    python analysis/reproduce.py --mode aadr --fetch \
      --source-root ../central-asia-sources \
      --output ../central-asia-aadr-reproduction

Install `requirements.txt` first. It includes the scikit-learn dependency used
by the v4 predictive sensitivity. The runner checks exact package versions;
uses the frozen seeds and replicate counts; records every command, log and
failure; and compares all regenerated statistical CSVs to the reference with
rtol=1e-9 and atol=1e-12. Outputs must go into an empty directory outside this
repository. Existing results are never overwritten by the runner.

For an offline statistical rerun from the included input:

    python analysis/reproduce.py --mode derived \
      --output ../central-asia-derived-reproduction

`aadr` verifies raw AADR extraction; `derived` begins with the included input.
Neither route redoes the AmtDB/aYChr cross-database audit or Figure 6. The
runner explicitly records inherited category definitions, summary metadata
and database-coverage counts. It does not copy previously calculated
statistical CSVs or figures into a new run. Literature/source-scope audits and
the separate wild-bootstrap calibration simulation are outside these reruns.
A request to database maintainers is not required to use or finish this
verified primary-analysis package. See `REPRODUCTION.md` for the optional
three-source route and its unresolved source requirement. A passing local
release verifier is not publication approval.

## Run the statistical extensions separately

The revised >=2 primary-profile analysis and associated TV/time-grid
extensions can also be regenerated in a new external directory:

    PYTHONDONTWRITEBYTECODE=1 python analysis_extensions/statistical_extensions_v4/run_stat_extensions.py \
      --input data/derived/central_asia_analysis_input_v1.csv \
      --outdir ../central-asia-extensions \
      --readme-output ../central-asia-extensions/README.md

    PYTHONDONTWRITEBYTECODE=1 python analysis_extensions/statistical_extensions_v4/test_results.py \
      --input data/derived/central_asia_analysis_input_v1.csv \
      --results ../central-asia-extensions \
      --readme ../central-asia-extensions/README.md

Choose a new extension directory. Use the `derived` runner above for the
complete statistical rerun and comparison with the frozen reference.

The included file is project-coded, not anonymous: combinations of public
AADR-derived attributes can remain linkable to the source resource. It is
therefore limited to fields actually needed by the released analyses.

## Optional three-source reproduction

This broader route is retained for future recovery of the historical source.
It is not a prerequisite for reproducing the primary analysis or extensions.
The saved AmtDB audit and Figure 6 must continue to be identified as
historical inherited results until that source is recovered and checked.
Contacting database authors is an optional recovery method, not a condition
of completing the current package. No newer AmtDB export was adopted.

Follow data/README.md to obtain the exact frozen inputs and verify their
SHA-256 hashes. AADR is pinned to Dataverse version 14.0 (version id 735358),
datafile id 13994518 and its immutable API endpoint. aYChr-DB is pinned to
commit `bc770a59ace8cd4c042c6f903d620d93ee751eb0`. The only remaining
raw-source blocker is AmtDB v1.009: an independent user needs either an
immutable archive URL or a deposit of the original hash-matching CSV. Once it
is available, run:

    python analysis/reproduce.py --mode raw --fetch \
      --source-root ../central-asia-sources \
      --amtdb /absolute/path/to/original-v1.009.csv \
      --output ../central-asia-three-source-reproduction

The full run validates the three frozen input hashes before parsing them.
Figure 4 reports original-sample site-balanced point estimates and bootstrap
percentile intervals; bootstrap medians remain separately labelled
diagnostics. Date outputs are assumption-based bin-assignment scenarios, not
calibrated radiocarbon posterior draws.

The original databases are not committed. Exact input hashes are recorded in
data/SOURCES.tsv and in the analysis manifest.

## Responsible release

The GitHub package contains one deliberately minimized project-coded AADR-
derived analytical input. It contains no source person IDs, locality strings,
skeletal fields, exact coordinates or terminal haplogroup calls. Original
catalogues, exact-ID audit tables and source-labelled site-lineage vectors
remain excluded. The extension includes project-coded site-profile vectors
deterministically derived from the already released analytical input; they add
no source site strings or terminal calls. Other summary tables are limited to
country-period, region-period and model outputs. The sampling map and
analytical input use
coordinates rounded to one decimal degree. Person-level rows derived from a
restricted preprint are not redistributed; only a citation-level audit is
included. See
THIRD_PARTY_NOTICES.md and PUBLIC_RELEASE_CHECKLIST.md before changing
repository visibility.

No new sampling or destructive analysis of human remains was performed.
Interpretation follows the principle that genetic ancestry, archaeological
culture, ethnicity and language are not interchangeable.

## Licensing and citation

Original software in analysis and tests is released under the MIT License.
Original documentation, figures and curated output tables are intended for
CC BY 4.0 release after coauthor approval. Third-party-derived records remain
subject to their source terms. See LICENSES and THIRD_PARTY_NOTICES.md.

Do not mint a DOI or tag version 1.0 until the author list, affiliations,
funding, archaeology review and public-release checklist are complete.
