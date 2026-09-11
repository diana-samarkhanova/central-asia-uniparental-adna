# Verified existing package — 11 September 2026

The owner selected completion on the existing verified analysis without a
request to database maintainers. This finalization retains v4.1 data, code and
figures. It does not change the source release to AmtDB v1.010.

| Component | Evidence and status |
|---|---|
| Primary AADR input | Independently rebuilt from exact AADR v66.p1; all 489 rows and the complete file SHA-256 match the released analytical input. |
| Main statistical results | All 20 compared main CSVs were reproduced byte for byte in the recorded frozen rerun. |
| v4 statistical extensions | All 12 extension CSVs were reproduced byte for byte in the same recorded rerun. |
| Figure 5 | Both regenerated PNGs were published with owner approval; v4 provenance identifies the actual source tables and generator. |
| AmtDB secondary audit | Historical v1.009 results retained from the earlier analysis. The original metadata CSV is missing; the cross-database audit, its coverage table and Figure 6 were not freshly reproduced. |
| Supplementary calibration | Retained from the earlier analysis; the separate simulation was not rerun. |

The recorded run used 489 individuals from 136 sites, with 438 mtDNA calls and
229 Y calls. AmtDB and aYChr records are excluded from primary AADR inference.
For the revised primary profiles with at least two calls, the Holm-adjusted
P values are 0.2224 for mtDNA and 0.0140 for Y. The analyses do not establish
that Y-chromosome turnover is faster than mtDNA turnover.

Use the documented `derived` route to reproduce the statistical results, or
the `aadr` route to additionally reconstruct the analytical input from the
original AADR annotation. Neither route is presented as a fresh three-source
reconstruction. Evidence is recorded in `validation/reproduction_v4_1.json`;
commands and retained metadata are described in `REPRODUCTION.md`.

Recovery of the historical AmtDB CSV is an optional future extension of this
package's reproducibility scope. The raw-source route still requires its exact
registered hash and stops if the file is absent or different. A live database
export must not be relabelled as v1.009.

The reproducibility package is complete within this stated scope. The PR
remains a draft pending authorship metadata, co-author approval, archaeology
interpretation review and the final package licensing review. No merge,
release tag or archival DOI is asserted here.
