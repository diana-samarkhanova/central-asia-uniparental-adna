# Released analytical input and frozen-source acquisition

`derived/central_asia_analysis_input_v1.csv` is the schema-locked input for
the public recomputation route. It contains 489 project-coded AADR-derived
rows from 136 normalized sites, 438 broad mtDNA calls and 229 broad Y calls.
It omits source person IDs, locality strings, skeletal fields, exact
coordinates and terminal haplogroup calls; see `derived/README.md` for the
field-level contract and linkage warning.

Scope finalized 11 September 2026: the included input and the full primary
analysis/extension computations are verified. Independent extraction from the
exact AADR annotation recovered the same 489-row input byte for byte. A
complete statistical rerun matched all 32 compared CSVs, and 39 local tests
passed with the AADR fixture; see `validation/reproduction_v4_1.json` and
`AS_IS_STATUS.md`.

No AmtDB update is adopted. The existing AmtDB v1.009 coverage/crosswalk audit
and Figure 6 remain historical inherited results; their original CSV is not
currently available for a fresh source audit. Contacting maintainers or
recovering that file is optional future work, not a condition for reproducing
the primary AADR analysis and extensions.

The instructions below are needed only for independent source extraction.
The `aadr` route requires AADR alone; the optional complete `raw` route needs
all three historical inputs.

Acquire the two resolved sources with:

    python analysis/fetch_sources.py --resource aadr --resource aychr \
      --output-root ../central-asia-sources

The downloader verifies pinned `download_url` entries in `SOURCES.tsv`.
Selecting only AADR and aYChr does not claim that all three sources are ready.
For the optional raw route, AmtDB remains explicitly blocked unless the
original CSV is provided with `--local amtdb=PATH`. See `../REPRODUCTION.md`
for the complete statistical rerun and independent raw-AADR extraction route.

The primary inferential analysis uses AADR; AmtDB and aYChr support the
separate database audit and do not add records to the primary cohort. Raw
copies of these upstream resources are intentionally excluded from Git.

1. Open AADR DOI 10.7910/DVN/FFIDCW, select **Dataverse version 14.0 (version
   id 735358) / AADR v66.p1**, and obtain `v66.p1_2M.aadr.PUB.anno` from the
   immutable datafile endpoint
   `https://dataverse.harvard.edu/api/access/datafile/13994518`. The verified
   file is 13,450,350 bytes, has MD5 `02a75f75de319829e89dd10a0d0f62c5`,
   matches the SHA-256 in `SOURCES.tsv`, and is distributed as CC0 1.0.
2. For the optional three-source route only, obtain the **frozen AmtDB
   v1.009 (2024-02-28) metadata CSV**. Do not silently
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
reproduction is claimed. This limitation does not block the verified primary
analysis and extension package. No new download or successful three-source
run is implied by retaining these acquisition instructions.
