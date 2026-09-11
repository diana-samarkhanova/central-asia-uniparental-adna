# Source, Russian-language coverage, and locality audit

Audit date: **2026-08-21** (Asia/Almaty)

Historical audit record. For the owner's 11 September 2026 decision to
complete the verified existing package with the AmtDB audit explicitly
inherited, see `../../AS_IS_STATUS.md`. The missing v1.009 export still limits
fresh three-source reconstruction; recovery is optional future work.

Scope: read-only resolution of review points A2, C11, and B10 against the
current frozen catalogue. No manuscript, analysis code, or release-package
file was changed.

## Integration decisions

| Issue | Decision | What may be said in the manuscript/release |
|---|---|---|
| AADR v66.p1 locator | **Resolved** | The exact released Dataverse version and datafile ID were independently resolved; the downloaded bytes reproduce the registered SHA-256. |
| AADR licence | **Resolved** | Dataset version 14.0 declares CC0 1.0 (`CC0-1.0`). |
| AmtDB v1.009 locator | **Still blocked** | No public immutable metadata export for v1.009 was found. The live service is v1.010. Do not claim that v1.009 is externally recoverable. |
| AmtDB licence | **Verified, but not a substitute for the missing file** | The site declares CC BY 4.0. If the exact old file matching the registered hash is recovered, redistribution in a versioned repository is permissible with attribution, licence link, and change statement. |
| Russian-language evidence layer | **Three reference-ready records; one conditional grey record** | Add Botai 2017, Berel 2020, and Tutkaul/Kaylu 2023 to the bibliography with the duplicate/quality qualifications below. Keep Urzhar 2016 outside the quantitative core unless the original full text is recovered. |
| Locality strings | **One safe alias; all numbered candidates remain distinct** | Preserve `locality_raw`; merge only `Biestamak` into `Bestamak`. Do not infer site identity from rounded coordinates or lexical similarity. |

Detailed evidence is in:

- `source_resolution.md`
- `russian_language_search.md`
- `locality_audit.md`
