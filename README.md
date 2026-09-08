# Central Asia uniparental ancient DNA

Reproducible secondary analysis of published ancient mitochondrial and
Y-chromosome assignments from Kazakhstan, Kyrgyzstan, Tajikistan,
Turkmenistan and Uzbekistan.

Status: private pre-submission repository, evidence frozen on 25 July 2026 and
statistical/literature/site-normalization corrections synchronized on 21
August 2026.
This is not yet a citable public release.

## Scope

The project harmonizes AADR v66.p1 metadata, audits marker-specific databases,
and describes the composition of published archaeological individuals and
equal-weighted published localities. The estimand is the published evidence
base, not the population frequency of lineages in ancient Central Asia.

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
association is not reproduced (71 profiles; Holm *P*=0.2224), whereas the Y
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
- tests: focused unit tests for the previously error-prone harmonization and
  resampling rules.

The unpublished manuscript is deliberately excluded until authorship,
affiliations and coauthor approval are confirmed.

## Quick integrity check

Create a Python 3.12 environment and install the frozen dependencies:

    python -m pip install -r requirements.txt

Then run:

    python -m unittest discover -s tests -v
    python analysis/test_aggregate_release.py results/aadr-v66p1_2026-07-25
    python analysis_extensions/statistical_extensions_v4/test_results.py
    python analysis/verify_release.py

These commands validate the project-coded input and aggregate release,
including the frozen counts, statistical summaries, code hashes and
50,000-replicate paired bootstrap.

## Recompute the released analysis from the included input

The public workflow can reproduce the headline analyses without any
person-level source catalogue or private rows. Copy the released result
directory so that the fixed category definitions and non-recomputed audit
table are available, then run:

    cp -R results/aadr-v66p1_2026-07-25 work/reproduced
    MPLCONFIGDIR=.mplconfig python analysis/recompute_from_catalogue.py \
      --catalogue data/derived/central_asia_analysis_input_v1.csv \
      --analysis-output work/reproduced \
      --bootstrap 2000 \
      --paired-bootstrap 50000 \
      --permutations 9999 \
      --callability-resamples 99999 \
      --date-draws 5000 \
      --seed 20260725 \
      --aggregate-only

    MPLCONFIGDIR=.mplconfig python analysis/run_global_sensitivities.py \
      --catalogue data/derived/central_asia_analysis_input_v1.csv \
      --analysis-output work/reproduced \
      --permutations 1999 \
      --seed 20260726

    python analysis/test_aggregate_release.py work/reproduced

The revised >=2 primary-profile analysis and associated TV/time-grid
extensions consume the same included input:

    python analysis_extensions/statistical_extensions_v4/run_stat_extensions.py
    python analysis_extensions/statistical_extensions_v4/test_results.py

The included file is project-coded, not anonymous: combinations of public
AADR-derived attributes can remain linkable to the source resource. It is
therefore limited to fields actually needed by the released analyses.

## Full reproduction

Follow data/README.md to obtain the exact frozen inputs and verify their
SHA-256 hashes. AADR is pinned to Dataverse version 14.0 (version id 735358),
datafile id 13994518 and its immutable API endpoint. aYChr-DB is pinned to
commit `bc770a59ace8cd4c042c6f903d620d93ee751eb0`. The only remaining
raw-source blocker is AmtDB v1.009: an independent user needs either an
immutable archive URL or a deposit of the original hash-matching CSV. Once it
is available, run:

    MPLCONFIGDIR=.mplconfig python analysis/run_analysis.py \
      --aadr data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno \
      --amtdb data/raw/amtdb_v1_009/amtdb_v1.009_metadata.csv \
      --aychr data/raw/aychr_db/a-YChr-DB_V5.xlsx \
      --outdir results/reproduced \
      --bootstrap 2000 \
      --paired-bootstrap 50000 \
      --permutations 9999 \
      --callability-resamples 99999 \
      --date-draws 5000 \
      --seed 20260725

    python analysis/run_global_sensitivities.py \
      --analysis-output results/reproduced \
      --permutations 1999

    python analysis/test_full_outputs.py results/reproduced

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
