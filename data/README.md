# Frozen input acquisition

The full analysis uses three public metadata resources. They are intentionally
excluded from Git because they are upstream-maintained datasets.

1. Download AADR v66.p1 from DOI 10.7910/DVN/FFIDCW and select the
   v66.p1 2M public annotation file.
2. Download the v1.009 metadata CSV from https://amtdb.org/.
3. Download a-YChr-DB - V5.xlsx from
   https://github.com/eelhaik/aYDB and rename it to
   a-YChr-DB_V5.xlsx locally.
4. Place the files at the paths recorded in SOURCES.tsv.
5. Verify them before analysis:

       python analysis/verify_inputs.py \
         --aadr data/raw/aadr_v66p1/v66.p1_2M.aadr.PUB.anno \
         --amtdb data/raw/amtdb_v1_009/amtdb_v1.009_metadata.csv \
         --aychr data/raw/aychr_db/a-YChr-DB_V5.xlsx

The AADR landing page is the authoritative release record. A stable direct
datafile identifier was not recorded in the frozen audit, so this repository
does not invent one. Record the Dataverse file ID and retrieval timestamp
before version 1.0.
