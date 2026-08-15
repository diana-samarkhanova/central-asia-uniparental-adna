# Central Asia uniparental ancient DNA

Reproducible secondary analysis of published ancient mitochondrial and
Y-chromosome assignments from Kazakhstan, Kyrgyzstan, Tajikistan,
Turkmenistan and Uzbekistan.

Status: private pre-submission repository, evidence frozen on 25 July 2026 and
statistical/literature audit corrections applied on 15 August 2026.
This is not yet a citable public release.

## Scope

The project harmonizes AADR v66.p1 metadata, audits marker-specific databases,
and describes the composition of published archaeological individuals and
equal-weighted published localities. The estimand is the published evidence
base, not the population frequency of lineages in ancient Central Asia.

The frozen AADR catalogue contains 501 unique archaeological individuals.
The primary interval, 3500 BCE to 1500 CE, contains 489 individuals from
137 localities, with 438 mitochondrial and 229 Y-chromosome calls. The
country-adjusted all-profile models detect period-associated structure, but
repeated-locality and sensitivity analyses show that it cannot be interpreted
as a region-wide demographic replacement. In 216 paired men from 98 sites,
the broad-L1 contrast Δ(Y−mtDNA) is −0.106381 (95% bootstrap interval −0.267291
to 0.002024); 0.05332 is an ordinary-bootstrap sign-stability diagnostic, not
a null-hypothesis *P* value. A 12-category AADR ISOGG-prefix sensitivity gives
Δ=−0.011292 (95% interval −0.154022 to 0.075681); the paired encoding-induced
change from broad L1 is +0.095089 (95% interval 0.023370 to 0.176052). Neither
encoding supports faster Y-chromosome turnover.

## Repository map

- analysis: scientific analysis, sensitivity analysis, input validation and
  integrity checks.
- data: acquisition instructions and frozen-source registry. Raw databases are
  intentionally not redistributed.
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
    python analysis/verify_release.py

These commands validate the aggregate release, including the frozen counts,
statistical summaries, code hashes and 50,000-replicate paired bootstrap.
They do not download the three source databases or rerun the extraction.
Exact regeneration starts from the upstream inputs described below; resulting
individual-level working tables must remain outside Git.

## Full reproduction

Follow data/README.md to obtain the exact frozen inputs and verify their
SHA-256 hashes. As of 15 August 2026 this route remains blocked for an
independent user because the exact AADR Dataverse datafile locator and an
immutable public AmtDB v1.009 archive URL have not been verified. The aYChr-DB
workbook is pinned to commit `bc770a59ace8cd4c042c6f903d620d93ee751eb0`.
Once the two missing locators or redistribution-cleared frozen inputs are
available, run:

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

The GitHub package contains no individual-level ancient-sample tables,
exact-ID audit table or site-period lineage vectors. Aggregate tables are
limited to country-period, region-period and model summaries. The sampling map
was generated from coordinates rounded to one decimal degree; exact source
coordinates are unnecessary for the statistical models and remain available
from the cited source databases. Person-level rows derived from a restricted
preprint are not redistributed; only a citation-level audit is included. See
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
