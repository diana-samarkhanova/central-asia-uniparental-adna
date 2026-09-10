# v4.1 reproduction tooling — 8 September 2026

This candidate continues the audited v4 dataset and statistical definitions.
It does not change the analytical cohort, analysis freeze, lineage encodings,
period boundaries or the interpretation of the manuscript.

## Changes

- Added `fetch_sources.py`: pinned source URLs, strict SHA-256 verification,
  atomic installation, bounded retries, local AmtDB import and machine-readable
  acquisition status. AADR and aYChr downloads have been exercised against the
  actual official files; AmtDB v1.009 remains unresolved.
- Added `extract_aadr_input.py`: rebuilds the primary analytical input directly
  from verified AADR and requires an exact match to the frozen input hash.
- Factored AADR catalogue preparation so extraction and the full raw pipeline
  use the same parsing, marker-specific deduplication, pooling and verified
  locality normalization.
- Fixed record-order differences between raw extraction and saved-catalogue
  recomputation. Category definitions retain their frozen order; records use
  the published catalogue order before assigning random draws.
- Fixed the full-output gate's assumption that every valid run must contain
  a `recomputed_from` field; direct raw runs instead validate all three source
  hashes and their actual source-code provenance.
- Added `reproduce.py`: isolated output directories, exact dependency checks,
  frozen input validation, explicit inherited metadata, stage logs, complete
  numeric CSV comparison and nonzero failure exits.
- Added explicit input/output/README paths to the statistical extensions.
  Isolated runs no longer rewrite a README beside the source script.
- Added missing `scikit-learn==1.8.0` to both environment definitions.
- Fixed aggregate checks to validate actual input provenance by resolved path
  and frozen hash, rather than requiring one hard-coded relative pathname.
- Added regression tests for source download failures, interrupted transfers,
  source/order/provenance defects, output preservation and result comparison.
- Corrected the README's stale private-repository label and retained the newer
  public locality-audit summary already present in GitHub.

## Validation completed

- All 39 tests passed locally with the verified AADR fixture enabled.
- AADR and aYChr downloads matched their pinned official files and frozen
  SHA-256 values.
- Independent AADR extraction reproduced 501 canonical individuals and the
  489-person primary cohort from 136 sites, with 438 mtDNA/229 Y calls.
  The coded input matched SHA-256
  `692a69cf38cc736f96ea5aa6b3b15024a9a49d50c03ae3e1470a8dc475504cc3`.
- A full isolated `derived` rerun passed every output gate and all 32 CSV
  comparisons: 20 main tables and 12 extension tables. Every compared CSV
  was byte-identical. The inherited database-coverage audit is excluded.
- CI discovers all test modules and checks saved output integrity. Actual
  AADR extraction is enabled by its optional source fixture.

See `validation/reproduction_v4_1.json`. This evidence covers independent
AADR reconstruction and the full statistical rerun; complete three-source
reproduction remains blocked by AmtDB v1.009.

## Boundaries

The AADR-only route does not recompute the separate AmtDB/aYChr cross-database
coverage audit. The all-three-source route stops until the original AmtDB CSV
matches its frozen hash. A newer live database is never silently substituted.

This tooling update neither completes author/archaeology approvals nor creates
a merged release, tag or DOI. Detailed raw catalogues are generated only into
external work directories and are not automatically published. The project-
coded analytical input remains byte-identical to v4.

## Figure publication — 10 September 2026

Following explicit owner approval, this update adds the verified regenerated
Figure 5 PNGs and the v4 figure's provenance to the public draft PR. The v4
figure generator records input-table, script and output hashes, raw P values
and rendering settings. All statistical CSVs remain byte-identical to v4.

The earlier code-only commit preserved the old figures while figure publication
awaited approval. This update completes that deferred figure publication; the
remaining AmtDB and final-release requirements above still apply.
