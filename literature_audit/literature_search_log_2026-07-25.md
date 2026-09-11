# Ancient uniparental DNA in Central Asia: reproducible search log

**Audit date:** 2026-07-25 (Asia/Almaty)  
**Targeted update check:** 2026-08-15 (Asia/Almaty)  
**Primary geography:** Kazakhstan, Kyrgyzstan, Tajikistan, Turkmenistan, Uzbekistan  
**Scope:** ancient human mtDNA and Y-chromosome evidence, including genome-wide studies with extractable uniparental calls  
**Design:** targeted evidence search and verified seed bibliography; not a completed PRISMA systematic/scoping review

## Status and limitations

This is a partly reproducible **targeted search protocol and verified seed
bibliography**, not a completed PRISMA screening flow. PubMed and Europe PMC
result-set counts were measured on 2026-07-25, but records were not fully
exported, deduplicated or dual-screened. The two searches within each database
overlap and must not be summed. Counts for subscription databases and
Russian-language indexes are recorded as `NA—not executed/exported`, not
estimated. The 2026-08-15 targeted update checked known key records and
metadata corrections; it did not rerun the four API searches, so the captured
counts below remain explicitly tied to 2026-07-25.

No PRISMA flow diagram, exhaustive study count, risk-of-bias assessment or
claim of systematic-review completeness should be derived from this log until
all planned sources are executed and the exported records are screened as
specified below.

## Search log

| Source | Search date | Search | Result count | Use |
|---|---|---:|---:|---|
| PubMed E-utilities | 2026-07-25 | P1 | 119 | Core ancient-DNA search |
| PubMed E-utilities | 2026-07-25 | P2 | 162 | Marker/legacy search |
| Europe PMC REST API | 2026-07-25 | E1 | 115 | Core ancient-DNA search |
| Europe PMC REST API | 2026-07-25 | E2 | 176 | Marker/legacy search |
| Crossref REST API | 2026-07-25 | C1 | 2,118,296 | DOI/title discovery only; ranked query, not a Boolean systematic-search count |
| Crossref REST API | 2026-07-25 | C2 | 8,964 | DOI/title discovery only; ranked query, not a Boolean systematic-search count |
| Google Scholar-style web search | 2026-07-25 | G1–G4 | NA—reliable reproducible count unavailable | Citation chasing/discovery only |
| Scopus | 2026-07-25 | S1/S2 | NA—not executed/exported | Must be run before submission |
| Web of Science Core Collection | 2026-07-25 | S1/S2 | NA—not executed/exported | Must be run before submission |
| eLIBRARY.ru | 2026-07-25 | R1/R2 | NA—not executed/exported | Russian-language grey/legacy literature |
| CyberLeninka | 2026-07-25 | R1/R2 | NA—not executed/exported | Russian-language grey/legacy literature |
| bioRxiv/Research Square | 2026-07-25 | S1/S2 | NA—not executed/exported | Preprints must remain a separate evidence stream |

### Targeted update after the count freeze

The following records were verified or corrected on 2026-08-15. They are not
included in the result counts above:

- Askapuli et al. (2026), Golden Horde elite genomes,
  DOI `10.1073/pnas.2531003123`, together with the terminal-Y interpretation
  exchange by Zhabagin and Sabitov, DOI `10.1073/pnas.2607193123`, and the
  reply by Askapuli, DOI `10.1073/pnas.2609024123`;
- Yang et al. (2026), eastern Tianshan, DOI `10.1093/molbev/msag057`;
- Moots et al. (2026), *Ancient DNA in motion*, arXiv `2608.09399`, a preprint
  submitted after the 2026-07-25 search cutoff and retained only as contextual
  methodological guidance.

Jeong et al. (2019), Kumar et al. (2021) and Wang et al. (2021) were already
present in the seed bibliography and were retained as key regional/contextual
records.

### P1 — PubMed core

```text
("ancient DNA"[Title/Abstract] OR archaeogenom*[Title/Abstract] OR archaeogenetic*[Title/Abstract] OR palaeogenom*[Title/Abstract] OR paleogenom*[Title/Abstract] OR paleogenetic*[Title/Abstract] OR "ancient genome"[Title/Abstract] OR "ancient genomes"[Title/Abstract])
AND
(Kazakhstan[Title/Abstract] OR Kazakh*[Title/Abstract] OR Kyrgyzstan[Title/Abstract] OR Kyrgyz*[Title/Abstract] OR Kirgiz*[Title/Abstract] OR Uzbekistan[Title/Abstract] OR Uzbek*[Title/Abstract] OR Tajikistan[Title/Abstract] OR Tajik*[Title/Abstract] OR Turkmenistan[Title/Abstract] OR Turkmen*[Title/Abstract] OR "Central Asia"[Title/Abstract] OR "Inner Asia"[Title/Abstract] OR Turan[Title/Abstract] OR "Eurasian steppe"[Title/Abstract] OR BMAC[Title/Abstract] OR "Bactria-Margiana"[Title/Abstract] OR Andronovo[Title/Abstract] OR Saka[Title/Abstract] OR Scyth*[Title/Abstract] OR Wusun[Title/Abstract] OR Usun[Title/Abstract] OR Kangju[Title/Abstract] OR Xiongnu[Title/Abstract] OR Turkic[Title/Abstract] OR "Golden Horde"[Title/Abstract] OR Botai[Title/Abstract] OR Berel[Title/Abstract] OR Gonur[Title/Abstract] OR Sarazm[Title/Abstract] OR Koken[Title/Abstract] OR "Boz-Barmak"[Title/Abstract])
AND
("1800/01/01"[Date - Publication] : "2026/07/25"[Date - Publication])
```

### P2 — PubMed marker/legacy

```text
(mtDNA[Title/Abstract] OR mitochondrial[Title/Abstract] OR mitogenom*[Title/Abstract] OR "hypervariable region"[Title/Abstract] OR HVR-I[Title/Abstract] OR HVS-I[Title/Abstract] OR Y-SNP[Title/Abstract] OR Y-STR[Title/Abstract] OR "Y chromosome"[Title/Abstract] OR uniparental[Title/Abstract])
AND
(ancient[Title/Abstract] OR archaeological[Title/Abstract] OR burial[Title/Abstract] OR cemetery[Title/Abstract] OR kurgan[Title/Abstract] OR remains[Title/Abstract])
AND
(Kazakhstan[Title/Abstract] OR Kazakh*[Title/Abstract] OR Kyrgyzstan[Title/Abstract] OR Kyrgyz*[Title/Abstract] OR Kirgiz*[Title/Abstract] OR Uzbekistan[Title/Abstract] OR Uzbek*[Title/Abstract] OR Tajikistan[Title/Abstract] OR Tajik*[Title/Abstract] OR Turkmenistan[Title/Abstract] OR Turkmen*[Title/Abstract] OR "Central Asia"[Title/Abstract] OR "Inner Asia"[Title/Abstract] OR Turan[Title/Abstract] OR "Eurasian steppe"[Title/Abstract] OR BMAC[Title/Abstract] OR "Bactria-Margiana"[Title/Abstract] OR Andronovo[Title/Abstract] OR Saka[Title/Abstract] OR Scyth*[Title/Abstract] OR Wusun[Title/Abstract] OR Usun[Title/Abstract] OR Kangju[Title/Abstract] OR Xiongnu[Title/Abstract] OR Turkic[Title/Abstract] OR "Golden Horde"[Title/Abstract] OR Botai[Title/Abstract] OR Berel[Title/Abstract] OR Gonur[Title/Abstract] OR Sarazm[Title/Abstract] OR Koken[Title/Abstract] OR "Boz-Barmak"[Title/Abstract])
AND
("1800/01/01"[Date - Publication] : "2026/07/25"[Date - Publication])
```

### E1/E2 — Europe PMC

For E1, prefix every P1 concept and geography term with `TITLE_ABS:` and replace the date expression with:

```text
FIRST_PDATE:[1800-01-01 TO 2026-07-25]
```

Exact E1:

```text
(TITLE_ABS:"ancient DNA" OR TITLE_ABS:archaeogenom* OR TITLE_ABS:archaeogenetic* OR TITLE_ABS:palaeogenom* OR TITLE_ABS:paleogenom* OR TITLE_ABS:paleogenetic* OR TITLE_ABS:"ancient genome" OR TITLE_ABS:"ancient genomes") AND (TITLE_ABS:Kazakhstan OR TITLE_ABS:Kazakh* OR TITLE_ABS:Kyrgyzstan OR TITLE_ABS:Kyrgyz* OR TITLE_ABS:Kirgiz* OR TITLE_ABS:Uzbekistan OR TITLE_ABS:Uzbek* OR TITLE_ABS:Tajikistan OR TITLE_ABS:Tajik* OR TITLE_ABS:Turkmenistan OR TITLE_ABS:Turkmen* OR TITLE_ABS:"Central Asia" OR TITLE_ABS:"Inner Asia" OR TITLE_ABS:Turan OR TITLE_ABS:"Eurasian steppe" OR TITLE_ABS:BMAC OR TITLE_ABS:"Bactria-Margiana" OR TITLE_ABS:Andronovo OR TITLE_ABS:Saka OR TITLE_ABS:Scyth* OR TITLE_ABS:Wusun OR TITLE_ABS:Usun OR TITLE_ABS:Kangju OR TITLE_ABS:Xiongnu OR TITLE_ABS:Turkic OR TITLE_ABS:"Golden Horde" OR TITLE_ABS:Botai OR TITLE_ABS:Berel OR TITLE_ABS:Gonur OR TITLE_ABS:Sarazm OR TITLE_ABS:Koken OR TITLE_ABS:"Boz-Barmak") AND FIRST_PDATE:[1800-01-01 TO 2026-07-25]
```

Exact E2:

```text
(TITLE_ABS:mtDNA OR TITLE_ABS:mitochondrial OR TITLE_ABS:mitogenom* OR TITLE_ABS:"hypervariable region" OR TITLE_ABS:"HVR-I" OR TITLE_ABS:"HVS-I" OR TITLE_ABS:"Y-SNP" OR TITLE_ABS:"Y-STR" OR TITLE_ABS:"Y chromosome" OR TITLE_ABS:uniparental) AND (TITLE_ABS:ancient OR TITLE_ABS:archaeological OR TITLE_ABS:burial OR TITLE_ABS:cemetery OR TITLE_ABS:kurgan OR TITLE_ABS:remains) AND (TITLE_ABS:Kazakhstan OR TITLE_ABS:Kazakh* OR TITLE_ABS:Kyrgyzstan OR TITLE_ABS:Kyrgyz* OR TITLE_ABS:Kirgiz* OR TITLE_ABS:Uzbekistan OR TITLE_ABS:Uzbek* OR TITLE_ABS:Tajikistan OR TITLE_ABS:Tajik* OR TITLE_ABS:Turkmenistan OR TITLE_ABS:Turkmen* OR TITLE_ABS:"Central Asia" OR TITLE_ABS:"Inner Asia" OR TITLE_ABS:Turan OR TITLE_ABS:"Eurasian steppe" OR TITLE_ABS:BMAC OR TITLE_ABS:"Bactria-Margiana" OR TITLE_ABS:Andronovo OR TITLE_ABS:Saka OR TITLE_ABS:Scyth* OR TITLE_ABS:Wusun OR TITLE_ABS:Usun OR TITLE_ABS:Kangju OR TITLE_ABS:Xiongnu OR TITLE_ABS:Turkic OR TITLE_ABS:"Golden Horde" OR TITLE_ABS:Botai OR TITLE_ABS:Berel OR TITLE_ABS:Gonur OR TITLE_ABS:Sarazm OR TITLE_ABS:Koken OR TITLE_ABS:"Boz-Barmak") AND FIRST_PDATE:[1800-01-01 TO 2026-07-25]
```

### C1/C2 — Crossref DOI discovery

Crossref `query.bibliographic` is relevance ranked and does not execute these as strict Boolean expressions. Its large totals must not enter a PRISMA flow.

```text
C1: query.bibliographic="ancient DNA" Kazakhstan Kyrgyzstan Uzbekistan Tajikistan Turkmenistan "Central Asia"
C2: query.bibliographic=ancient archaeological mtDNA mitochondrial "Y chromosome" Kazakhstan Kyrgyzstan Uzbekistan Tajikistan Turkmenistan
filter=from-pub-date:1800-01-01,until-pub-date:2026-07-25
```

### G1–G4 — Google Scholar-style discovery/citation chasing

```text
G1: "ancient DNA" (Kazakhstan OR Kyrgyzstan OR Uzbekistan OR Tajikistan OR Turkmenistan)
G2: (mtDNA OR mitochondrial OR "Y chromosome" OR Y-STR OR Y-SNP) ancient (Kazakhstan OR Kyrgyzstan OR Uzbekistan OR Tajikistan OR Turkmenistan)
G3: ("ancient genome" OR archaeogenomics) ("Central Asia" OR "Eurasian steppe" OR BMAC OR Saka OR Scythian)
G4: ("древняя ДНК" OR палеогенетика OR археогенетика OR мтДНК OR "Y-хромосома") (Казахстан OR Кыргызстан OR Узбекистан OR Таджикистан OR Туркменистан)
```

### R1/R2 — Russian-language indexes

```text
R1: ("древняя ДНК" OR палеогенетик* OR археогенетик*) AND (Казахстан OR Кыргызстан OR Киргиз* OR Узбекистан OR Таджикистан OR Туркменистан OR "Центральная Азия" OR "Средняя Азия")
R2: (митохондриальн* OR мтДНК OR "Y-хромосом*" OR гаплогрупп*) AND (древн* OR археолог* OR погребен* OR могильник OR курган) AND (Казахстан OR Кыргызстан OR Киргиз* OR Узбекистан OR Таджикистан OR Туркменистан OR БМАК OR Бактрийско-Маргиан* OR Андронов* OR сак* OR скиф* OR усун* OR кангюй* OR хунн* OR тюрк* OR "Золотая Орда")
```

## Inclusion and exclusion

Include human archaeological remains recovered in the five republics when an mtDNA/Y result or a genome-wide dataset with an extractable uniparental call is available and person-level site/date/provenance can be recovered. Include peer-reviewed articles; retain preprints in a separate evidence stream. Human host data from pathogen studies are eligible if reusable.

Exclude modern-only studies, animal/pathogen-only data without reusable human host genotypes, methods papers without regional persons, commercial genealogy, secondary claims without extractable primary records, and records with unresolvable geography/provenance. Comparator data from Xinjiang, Altai–Sayan, Mongolia, the southern Urals, northern Afghanistan and northeastern Iran form a predefined secondary stratum and must not be pooled into five-republic estimates.

## Deduplication and screening

1. Deduplicate citations by DOI, PMID and normalized title/author/year.
2. Dual-screen titles/abstracts and full texts; adjudicate disagreements; record a single explicit reason for each full-text exclusion.
3. Deduplicate people using AADR Individual ID plus skeletal/sample aliases, site, grave/context, date, molecular sex, publication and sequence accession.
4. Retain the best-quality genomic representation for analysis, while preserving every alias/source in a provenance table.
5. Do not merge people on site/period/haplogroup alone.
6. Keep HVR/RFLP and Y-STR-only evidence in a separate legacy tier; do not treat it as equivalent to complete mitogenomes or sequence-based Y calls.
7. PRISMA reporting must separate: bibliographic records, database/sample records, study-level duplicates, person-level duplicates, QC exclusions, kinship/site-balanced analytic sets.

## Safe gap statement

> Within the sources actually searched through 2026-07-25 and the targeted
> update through 2026-08-15, we did not identify an up-to-date all-five-republic
> synthesis accompanied by a reproducible, deduplicated person-level catalogue
> of ancient mtDNA and Y-chromosome lineages that explicitly handles duplicate
> genomic representations, kinship/cemetery clustering, call quality, uneven
> country-period coverage, legacy marker-only evidence and studies released
> after the latest AADR freeze.

This is a bounded finding from an incomplete targeted search, not proof that no
relevant publication exists and not a PRISMA-complete systematic-review claim.
