# Public release checklist

- Confirm article author order, affiliations, corresponding author and ORCIDs.
- Confirm CRediT contributions, funding, acknowledgements and competing
  interests.
- Obtain archaeology-specialist review of period and cultural interpretations.
- Confirm package-specific Poseidon reuse terms and immutable commit.
- Keep the verified AADR v66.p1 locator pinned: Dataverse version 14.0
  (version id 735358), datafile id 13994518, immutable API endpoint, size,
  MD5 and SHA-256 recorded in `data/SOURCES.tsv`.
- Supply and verify an immutable public AmtDB v1.009 metadata-export URL, or a
  redistribution-cleared frozen snapshot matching the recorded SHA-256.
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
- Repeat the full frozen-input run and compare all release checksums.
- Complete institutional Scopus/Web of Science exports and independent
  screening before claiming PRISMA systematic-review status.
- Add CITATION.cff and Zenodo metadata only after authorship approval.
- Create a versioned release, archive it in Zenodo and add the DOI to Data and
  Code Availability.
- Keep interpretations limited to published individuals and sites.
