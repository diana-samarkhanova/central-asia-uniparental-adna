# Public release checklist

Scope clarification — 11 September 2026: the existing primary AADR analysis
and statistical extensions have been verified, including independent recovery
of the frozen 489-row input, 32 byte-identical CSV comparisons and 39 local
tests with the AADR fixture. Both regenerated Figure 5 PNGs and the v4
provenance are already public. The AmtDB cross-database audit is retained as
historical inherited output. No newer AmtDB version or fresh complete
three-source run is claimed. Author contact is not a condition of completing
this reproducibility package. See `AS_IS_STATUS.md` for the scope record.

- Confirm article author order, affiliations, corresponding author and ORCIDs.
- Confirm CRediT contributions, funding, acknowledgements and competing
  interests.
- Obtain archaeology-specialist review of period and cultural interpretations.
- Confirm package-specific Poseidon reuse terms and immutable commit.
- Keep the verified AADR v66.p1 locator pinned: Dataverse version 14.0
  (version id 735358), datafile id 13994518, immutable API endpoint, size,
  MD5 and SHA-256 recorded in `data/SOURCES.tsv`.
- Keep the AmtDB coverage/crosswalk audit and Figure 6 labelled as historical
  inherited outputs, outside the verified fresh primary-analysis scope.
  Recovery of the exact AmtDB v1.009 snapshot is optional future work. Before
  claiming complete three-source reproduction, supply its hash-matching
  bytes or verified immutable archive and pass the full raw-route checks.
  Do not substitute a newer database version under the historical label.
- Keep Ulytau person-level preprint rows excluded unless reuse permission is
  documented.
- Confirm that the sole person-row file is
  `data/derived/central_asia_analysis_input_v1.csv`, that it passes the exact
  project-code/schema/count checks, and that no source identifier, locality,
  skeletal field, terminal call, exact coordinate, CKZ/private row or exact-ID
  crosswalk is tracked. The extension's site-profile vectors must use only the
  released `SITE####`/`STUDY###` project codes and broad categories.
- Confirm that the public input and sampling map use one-decimal coordinates
  and retain the project-coded-not-anonymous linkage warning.
- Preserve the recorded full primary/extension rerun evidence and compare
  release checksums against the actual candidate. Do not describe passing
  saved-output checks or CI as a new complete three-source run.
- Complete institutional Scopus/Web of Science exports and independent
  screening before claiming PRISMA systematic-review status.
- Add CITATION.cff and Zenodo metadata only after authorship approval.
- Create a versioned release, archive it in Zenodo and add the DOI to Data and
  Code Availability.
- Keep interpretations limited to published individuals and sites.
