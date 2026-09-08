# Project-coded analytical input

`central_asia_analysis_input_v1.csv` is a minimal AADR-derived input for the
catalogue-level public recomputation route. It contains 489 primary-analysis
rows and preserves only the variables required by the released models and
sensitivity analyses.

The file omits AADR genetic/individual IDs, skeletal identifiers and elements,
locality strings, terminal mtDNA/Y calls, free-text group labels, exact
coordinates and the exact-ID cross-database audit. Records, sites and studies
use release-specific project codes; only precomputed marker-specific kin-
representative flags are retained. Coordinates are rounded to 0.1 degree. No
CKZ/private post-v66 row is present.

This is **project-coded, not anonymous**. Public combinations of country,
period, broad lineage, study and rounded geography may permit linkage back to
the CC0 AADR source. The table must therefore be described as a reproducibility
input rather than as anonymized human-subject data. It is derived from the
frozen AADR v66.p1 public annotation recorded in `../SOURCES.tsv`; cite AADR
and the original component studies.

The author-side generation command is:

```bash
python analysis/build_public_analysis_input.py \
  --input results/aadr-v66p1_2026-07-25/tables/aadr_primary_analysis_catalogue.csv \
  --output data/derived/central_asia_analysis_input_v1.csv
```

The public recomputation workflow consumes the generated file directly and
does not require any source identifier crosswalk.
