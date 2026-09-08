# A2 — immutable source resolution

Audit date: **2026-08-21**

## AADR v66.p1: blocker can be removed

The target file was resolved through the Harvard Dataverse Native API, not by
guessing a download path.

| Field | Verified value |
|---|---|
| Dataset persistent identifier | `doi:10.7910/DVN/FFIDCW` |
| Dataset landing page | <https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi%3A10.7910%2FDVN%2FFFIDCW> |
| Exact version API | <https://dataverse.harvard.edu/api/datasets/:persistentId/versions/14.0?persistentId=doi:10.7910/DVN/FFIDCW> |
| Dataset version | `14.0` |
| Dataset-version ID | `735358` |
| State and release time | `RELEASED`; `2026-06-08T20:36:49Z` |
| Target filename | `v66.p1_2M.aadr.PUB.anno` |
| Dataverse datafile ID | `13994518` |
| Direct datafile-ID access URL | <https://dataverse.harvard.edu/api/access/datafile/13994518> |
| File version / datasetVersionId | `1` / `735358` |
| Restricted | `false` |
| Size | `13,450,350` bytes |
| Dataverse MD5 | `02a75f75de319829e89dd10a0d0f62c5` |
| Independently computed MD5 | `02a75f75de319829e89dd10a0d0f62c5` |
| Independently computed SHA-256 | `98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c` |
| Existing `SOURCES.tsv` SHA-256 | `98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c` |
| Version-level licence | `CC0 1.0`; SPDX `CC0-1.0`; <http://creativecommons.org/publicdomain/zero/1.0> |

The direct URL returned bytes whose MD5 and SHA-256 both match the independent
expectations. The URL should be recorded together with dataset DOI, exact
released version, datafile ID, size, and both checksums. This is sufficient to
remove the present `BLOCKED` text from the AADR row.

### Integration-ready AADR registry values

```text
release = Dataverse v14.0; files v66.p1; released 2026-06-08
landing_url = https://doi.org/10.7910/DVN/FFIDCW
immutable_file_url_or_commit = https://dataverse.harvard.edu/api/access/datafile/13994518
sha256 = 98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c
license = CC0 1.0 (SPDX: CC0-1.0)
provenance_note = dataset version ID 735358; datafile ID 13994518; file version 1; 13,450,350 bytes; MD5 02a75f75de319829e89dd10a0d0f62c5
```

Recommended availability wording:

> The AADR v66.p1 2M annotation file was obtained from Harvard Dataverse
> dataset version 14.0 (dataset DOI 10.7910/DVN/FFIDCW; dataset-version ID
> 735358; datafile ID 13994518). Its SHA-256 is
> 98eec5d897a6feedd274c91b476c4d416e60d12d171409018925522679ba316c.
> The released dataset version declares CC0 1.0.

## AmtDB v1.009: blocker cannot yet be removed

### What was verified

- The AmtDB changelog exposed at
  <https://amtdb.org/help?version=v1.009#v1.009> identifies **v1.009,
  28/02/2024**.
- The live homepage <https://amtdb.org/> now reports **v1.010, 11/08/2026,
  3,759 samples**. A current download therefore is not the registered v1.009
  input.
- The live site footer declares **CC BY 4.0** and links to
  <https://creativecommons.org/licenses/by/4.0/>. That licence permits sharing
  and adaptation subject to attribution, licence link, and indication of
  changes.
- Public web and archive discovery located historical v1.009 pages, but no
  frozen metadata-download object with a stable archive/DOI and no public file
  that could be verified against the registered SHA-256
  `531e8ee8fae181124f5a9b77b6fe8d677e64e35b815be2a3965020244fe31057`.
- The exact raw v1.009 CSV is absent from the inspected workspace and attached
  reproducibility packages. Only a downstream regional derivative is present.

### Consequence

The missing-file blocker is a provenance/recoverability problem, not chiefly a
licensing problem. CC BY 4.0 does not make an absent historical byte stream
reproducible.

Choose one of the following honest routes:

1. **Preferred:** recover the exact original CSV from the authors' machine or
   backup, verify the registered SHA-256, and deposit that unchanged file in a
   versioned Zenodo/repository release. Cite Ehler et al. (2019), name AmtDB
   v1.009 and the retrieval date, link CC BY 4.0, state whether the file is
   unchanged, and credit the database maintainers.
2. Ask the maintainers for an official archived v1.009 export. The primary
   contact listed by AmtDB is Edvard Ehler (`edvard.ehler@img.cas.cz`).
3. If the exact export cannot be recovered, publish only a **frozen regional
   derivative** whose rows, fields, creation logic, source attribution, and
   checksum are documented. State expressly that the full v1.009 export is not
   publicly recoverable and that no independent full-database reconstruction
   is claimed.
4. Alternatively, switch to v1.010 and rerun every dependent extraction,
   deduplication, table, figure, and statistic. Do not relabel a v1.010 export
   as v1.009.

Recommended limitation wording if route 3 is used:

> The legacy mitochondrial layer uses a frozen, attributed regional derivative
> prepared from AmtDB v1.009 (28 February 2024). During the reproducibility
> audit, the live AmtDB service had advanced to v1.010 and no immutable public
> archive of the full v1.009 metadata export could be verified. Accordingly,
> we release and checksum only the regional derivative and do not claim that
> the complete historical database export can be reconstructed from the live
> service.

