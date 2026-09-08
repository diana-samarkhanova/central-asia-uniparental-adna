# Released analytical input and frozen-source acquisition

`derived/central_asia_analysis_input_v1.csv` is the schema-locked input for
the public recomputation route. It contains 489 project-coded AADR-derived
rows from 136 normalized sites, 438 broad mtDNA calls and 229 broad Y calls.
It omits source person IDs, locality strings, skeletal fields, exact
coordinates and terminal haplogroup calls; see `derived/README.md` for the
field-level contract and linkage warning.

The instructions below are needed only for a new extraction from the three
upstream databases.

`python analysis/fetch_sources.py --output-root ../central-asia-sources`
downloads and verifies every source with a pinned `download_url` in
`SOURCES.tsv`; AmtDB remains explicitly blocked unless the original CSV is
provided with `--local amtdb=PATH`. Use `--resource aadr --resource aychr` to
acquire the two resolved sources. See `../REPRODUCTION.md` for a complete
isolated run and the independent raw-AADR extraction route.

The full analysis uses three public metadata resources. They are intentionally
excluded from Git because they are upstream-maintained datasets.

1. Open AADR DOI 10.7910/DVN/FFIDCW, select **Dataverse version 14.0 (version
   id 735358) / AADR v66.p1**, and obtain `v66.p1_2M.aadr.PUB.anno` from the
   immutable datafile endpoint
   `https://dataverse.harvard.edu/api/access/datafile/13994518`. The verified
   file is 13,450,350 bytes, has MD5 `02a75f75de319829e89dd10a0d0f62c5`,
   matches the SHA-256 in `SOURCES.tsv`, and is distributed as CC0 1.0.
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
landing page alone is insufficient. Stop if any hash differs. The AADR DOI,
Dataverse version and datafile endpoint are now pinned. The remaining upstream
identity blocker is AmtDB: its interactive download endpoint is not an
immutable v1.009 archive, so the original hash-matching CSV must be deposited
or an immutable archival locator supplied before end-to-end raw-source
reproduction is claimed.
