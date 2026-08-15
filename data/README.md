# Frozen input acquisition

The full analysis uses three public metadata resources. They are intentionally
excluded from Git because they are upstream-maintained datasets.

1. Open AADR DOI 10.7910/DVN/FFIDCW, select **Dataverse version 14.0 / AADR
   v66.p1**, and obtain `v66.p1_2M.aadr.PUB.anno`. Before release, record the
   Dataverse datafile ID and direct immutable URL; neither was independently
   verified during this audit.
2. Obtain the **frozen AmtDB v1.009 (2024-02-28) metadata CSV**. Do not silently
   substitute the live download: AmtDB moved to v1.010 on 2026-08-11, and no
   immutable public archive URL for the v1.009 export was verified. Use a
   locally retained/coauthor snapshot only if its SHA-256 matches
   `SOURCES.tsv`; otherwise raw-source reproduction is blocked.
3. Download `a-YChr-DB - V5.xlsx` from the verified upstream commit
   `bc770a59ace8cd4c042c6f903d620d93ee751eb0`:
   https://github.com/eelhaik/aYDB/blob/bc770a59ace8cd4c042c6f903d620d93ee751eb0/a-YChr-DB%20-%20V5.xlsx
   and rename it to `a-YChr-DB_V5.xlsx` locally.
4. Place the files at the paths recorded in SOURCES.tsv.
5. Verify them before analysis:

       python analysis/verify_inputs.py \
         --aadr data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno \
         --amtdb data/raw/amtdb_v1_009/amtdb_v1.009_metadata.csv \
         --aychr data/raw/aychr_db/a-YChr-DB_V5.xlsx

The expected hashes identify the files actually used; a matching filename or
landing page alone is insufficient. Stop if any hash differs. The AADR landing
page is the authoritative release record, but this repository does not invent
an unverified Dataverse file ID. Likewise, the interactive AmtDB download
endpoint is not an immutable v1.009 archive. Both locators are explicit release
blockers in `SOURCES.tsv`.
