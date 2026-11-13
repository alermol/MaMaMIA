# MaMaMIA — Russian Terminology Decision Record

- **Status:** Terminology rulings + copy decisions (no code written).
- **Companion files:** `russian-localization.md` (implementation spec) and `russian-copy.md` (the EN→RU string table this glossary backs).
- **Audience:** the engineer implementing the i18n dictionary, and the reviewer checking biological correctness.
- **Convention:** each ruling gives the chosen Russian, grammatical notes, rejected alternatives + reason, and whether the term must stay in Latin script.

---

## 0. Global rulings

1. **Latin-script terms (do NOT translate):** `BAM`, `BEDPE`, `FDR`, `GC`, `MAPQ`, `TSV`, `CSV`, `HTML`, `SVG`, `ZIP`, `CSI`, `BAI`, `CIGAR`, `FLAG`, `RNAME`, `POS`, `QUAL`, `SEQ`, `@RG`, `@SQ`, `@HD`, `@PG`, `@CO`, `SM`, `ID`, `bamhusk`, `MaMaMIA`. Rationale: these are file formats, BAM-field names, tools or statistical acronyms with no established Russian rendering; forcing Cyrillic would break searches, exports and byte-asserted tests. `GC` stays Latin when it means *GC content* ("GC-состав", "GC-индекс"); "GC" is never spelled "ГЦ".
2. **Action buttons: infinitive, not imperative.** Use "Загрузить", "Выбрать", "Запустить анализ", "Скачать" — never "Загрузите"/"Выберите". Rationale: Russian UI convention for neutral, non-gendered action verbs; imperative ("Загрузите") is the polite second-person form and reads as an instruction to a person, not a control label. This is applied uniformly across `sidebar.*`, `params.*`, `results.*` and modals.
3. **Capitalisation: sentence case.** Russian headings/labels capitalise only the first letter of the sentence/standalone label ("Параметры анализа", "Мин. MAPQ"). No Title Case ("Параметры Анализа" is wrong). Acronyms keep their native case (MAPQ, FDR, BAM).
4. **Ё/ё is used** everywhere ("учёт", "ещё", "подтверждённый", "размещённый"). No "е" replacement.
5. **Number formatting is done in R, not in copy.** Russian copy must not contain hardcoded numerals with separators. R currently formats numbers with English conventions (`format(big.mark = ",", …)` → "1,000"; `formatC(…, digits = 2)` → "1.23" dot-decimal). This leaks into the Russian UI and is flagged in §4 of this record and in `russian-copy.md`; it needs a locale-aware formatter (space thousands, comma decimal) as a follow-up, not in the dictionary.
6. **Identifiers, column names, export headers, tool names, file formats** are never translated (out of scope).

---

## 1. Required terms (glossary)

| # | English | Chosen Russian | Grammar | Rejected alternatives | Latin? |
|---|---|---|---|---|---|
| 1 | introgression | **интрогрессия** | fem., "интрогрессия/интрогрессии/интрогрессий" | "интрогрессия" is the established term (no calque needed). "привнос" is wrong register | no |
| 2 | breadth-aware mapping | **картирование по широте покрытия** (noun phrase); header note: "Картирование интрогрессий по широте покрытия" | — | **"с учётом покрытия" is rejected**: it overstates — "покрытие" (coverage) in Russian usually means read depth, but "breadth" here measures the *fraction of a window covered by ≥1 read* (see §3.4). "с учётом глубины" is factually wrong | no |
| 3 | karyotype | **кариотип** | masc. | "кариотип" is standard; no alternative | no |
| 4 | call / calls (primary noun = a detected introgression segment) | **находка** ("находка/находки/находок") | fem. | **"вызов" rejected**: calque; collides with "вызов"="challenge/summons". "событие" is acceptable in prose but too generic for a table/tab. "сегмент" is too narrow (loses "detected, calibrated event"). "находка" = "finding", used in RU bioinformatics for detected events | no |
| 5 | no-call / unassigned | **нет находки** (no call at a position); **без размещения** / **не размещено** (unassigned recipient) | — | "без вызова" (calque) rejected. "неназначенный" is understandable but "не размещено" matches the app's "placement" vocabulary | no |
| 6 | donor / recipient (of introgression) | **донор / реципиент** | masc. both; "донорный"/"реципиентный" as adjectives | standard; "акцептор" rejected (wrong register for genomics) | no |
| 7 | unplaced donor | **неразмещённый донорный контиг** (chromosome state); short: "донор (не размещён)" | — | see §3.1: it is a *donor contig drawn as a neutral bar*, not a call. "неразмещённый донор" (without "контиг") is acceptable as a short label but must not be read as "a donor call" | no |
| 8 | deleted recipient | **делетированная хромосома реципиента** (chromosome state); short: "реципиент с делецией" | — | see §3.2: a *recipient chromosome with no/negligible coverage*, not a donor segment. "удалённый реципиент" rejected (reads as "removed from the list") | no |
| 9 | component evidence | **(n/a — see §3.3)** | — | the app has no literal "component evidence"; the real fields are `class_evidence` (which side drove the class) and `component_status` (consistency across components). See §3.3 | — |
| 10 | segment | **сегмент** | masc. | standard | no |
| 11 | window | **окно** | neut.; "окно/окна/окон" | standard (sliding window) | no |
| 12 | bin | **бин** (masc., indeclinable-ish; "бины") | — | "корзина" rejected (e-commerce calque); in genomics RU "бин" is used for fixed-size bins | no (transliterated) |
| 13 | resolution | **разрешение** | neut. | standard (windows per chromosome) | no |
| 14 | coverage / read depth | **покрытие** (umbrella) / **глубина покрытия** or **число прочтений** (read depth specifically) | neut. | keep them distinct from "breadth" (§3.4) | no |
| 15 | breadth | **широта покрытия** (fraction of window covered); gloss "доля покрытых оснований" | fem. | see §3.4. "ширина" rejected (wrong), "охват" too vague | no |
| 16 | GC content | **GC-состав** | masc. | "ГЦ-состав" rejected (GC stays Latin) | GC only |
| 17 | masking / masked region | **маскирование / маскированный участок** | — | "маскировка" (camouflage) rejected | no |
| 18 | calibration | **калибровка** | fem. | standard | no |
| 19 | FDR | **FDR** | — | keep Latin; occasionally gloss "уровень ложных находок" but do not replace the acronym | yes |
| 20 | classification / class | **классификация / класс** | — | standard | no |
| 21 | confidence | **достоверность** (class confidence) | fem. | "уверенность" is a literal calque and less idiomatic for a classifier; "достоверность класса" preferred. Filter label "Confidence" → "Достоверность" | no |
| 22 | score | **оценка** (a score value); the scoring *method* stays "score/скор" as a metric label "Score" → "Оценка" | fem. | "скор" (loanword) is used colloquially but "оценка" is the neutral, standard term; "балл" rejected (exam register) | no |
| 23 | mapping quality | **MAPQ** (keep); gloss "качество картирования" | — | do not translate MAPQ; "качество выравнивания" is wrong (MAPQ is a mapping confidence, not alignment identity) | MAPQ only |
| 24 | alignment | **выравнивание** | neut. | standard | no |
| 25 | assembly | **сборка** (genome assembly) | fem. | standard | no |
| 26 | contig | **контиг** | masc. | standard transliteration | no |
| 27 | chromosome | **хромосома** | fem. | standard | no |
| 28 | unplaced/unlocalised contig | **неразмещённый контиг** (see #7) | — | "нелокализованный" is the term for *unlocalised* (placed on chr but position unknown) in some assemblies; MaMaMIA uses "unplaced" for "no recipient host", so use "неразмещённый" | no |
| 29 | sample | **образец** | masc. | standard | no |
| 30 | cohort | **когорта** | fem. | standard | no |
| 31 | bundle | **архив результатов** / **пакет** (in "results bundle"); keep "bundle" transliterated as "бандл" only if needed — **rejected**, use "архив" | — | "бандл" is a calque; "пакет" collides with R packages; "архив" (ZIP) is the accurate word | no |
| 32 | BEDPE | **BEDPE** | — | keep Latin | yes |
| 33 | BAM | **BAM** | — | keep Latin | yes |
| 34 | mosaic / mosaicism | **мозаицизм / мозаичный** | masc. | standard | no |
| 35 | subclonal | **субклональный** | — | standard | no |
| 36 | recipient background | **реципиентный фон** | — | "фон реципиента" preferred | no |
| 37 | filtering | **фильтрация** | fem. | standard | no |
| 38 | threshold | **порог** | masc. | standard | no |
| 39 | parameter | **параметр** | masc. | standard | no |
| 40 | scoring | **подсчёт оценок / скоринг** | — | "скоринг" is accepted colloquially but "подсчёт оценок" is safer for headings; keep "scoring" as a method name only | no |
| 41 | evidence | **свидетельства / данные** (as in "Evidence" section → "Свидетельства") | — | "доказательства" is too strong (these are supportive signals, not proof); "свидетельства" is the standard softer term | no |

---

## 2. Specific copy rulings

### 2.1 The name and subtitle

- **`MaMaMIA` stays untranslated** — it is a proper name and the acronym.
- **Subtitle "Max Mapping for Meticulous Introgression Analysis" stays in English.** This is the *expansion of the acronym* (Ma-x Ma-pping … M-eticulous I-ntrogression A-nalysis): translating it ("Максимальное картирование для тщательного анализа интрогрессии") destroys the backronym, and a Russian subtitle would read as a slogan that no longer explains the name. It is confirmed as the official expansion in `uvr.toml:4` (`description = "Max Mapping for Meticulous Introgression Analysis"`) and `README.md:8`.
- **Header note "Breadth-aware introgression mapping" IS translated** → "Картирование интрогрессий по широте покрытия" (this is descriptive copy, not part of the acronym).
- The illustrative Russian subtitle in `russian-localization.md:61` ("Максимальное картирование для тщательного анализа интрогрессии") is **superseded** by this ruling: keep English. If product later wants a Russian slogan, a *separate*, non-backronym line can be added; it must not replace the subtitle.

### 2.2 Latin/English retention rule (confirmed)

`BAM`, `BEDPE`, `FDR`, `GC`, `TSV`, `CSV`, `HTML`, `SVG`, `ZIP`, `MAPQ` all remain Latin. Rule: a term stays Latin when it is (a) a file format/extension, (b) a BAM field/header tag, (c) a statistical acronym with an established Latin form, or (d) a tool name. When such a term is *qualified* in Russian, the qualifier is declined normally and the Latin token stays fixed ("GC-состав", "BAM-файл", "FDR-порог").

### 2.3 Plural agreement (1 / 2-4 / 5+)

No plural machinery in Phase 1. High-visibility strings that join a bare number to a noun are rephrased to "Label: N" form (number last) so no agreement is needed. The specific strings are enumerated in `russian-copy.md` §"Plural-avoidance rewrite list", including:

- `"%d segment%s omitted"` → "Опущено сегментов: %d"
- `"%d unplaced hidden"/"%d unplaced shown"` → "Скрыто неразмещённых: %d" / "Показано неразмещённых: %d"
- `"… call(s) selected for BEDPE export"` → "Выбрано для BEDPE-экспорта: %d"
- `"%d BAM file%s"` (footer + preflight log) → "BAM-файлов: %d"
- `"Analysis finished: %d %s"` (calibrated call/calls) → "Завершено — откалиброванных находок: %d"
- `"%d target(s) not in assembly '%s'"` → "Целей вне сборки '%s': %d"
- `"Too many BAM files: %d (limit %d)"` → "Слишком много BAM-файлов: %d (предел %d)"

The English templates keep their `%s`/`%d` order so `sprintf` call sites are unchanged; only the Russian template is rephrased.

### 2.4 Number formatting

UI copy presents numbers **untouched** (the app formats them in R). No Russian string in the dictionary contains a hardcoded "1,234.5". Flagged leaks (out of scope for the dictionary, require an R-side formatter follow-up): `format(big.mark = ",", scientific = FALSE)` produces English thousands separators ("1,000"); `formatC(…, format="f", digits=…)` produces dot decimals ("1.23"); `call_tooltip_signed` and `call_tooltip_bp` in `R/call-tooltip.R` do the same. Recommendation: introduce `mamamia_format_number(value, lang)` returning "1 000" and "1,23" for `ru` and routing through it at the ~20 `format(...)` call sites, rather than putting any number in copy.

### 2.5 Action buttons

Infinitive, consistently (see §0.2). "Run analysis" → "Запустить анализ" (not "Запустите анализ"); "Choose BAM file" → "Выбрать BAM-файл"; "Download results bundle (.zip)" → "Скачать архив результатов (.zip)".

### 2.6 Capitalisation and ё

Sentence case everywhere; use "ё". See §0.3–0.4.

---

## 3. Data/domain correctness checks (spec-vs-source)

### 3.1 "unplaced donor" is a contig/chromosome state, NOT a call

- Source: `www/htmlwidgets/mamamiaKaryotypeCard.js:406-417` builds the no-call card title from the *chromosome* flags: `chromosome.unplaced ? "unplaced donor" : chromosome.deleted ? "deleted recipient" : chromosome.added ? "added donor" : "recipient chromosome"`. The State row is "unplaced donor chromosome" (`:427-429`). The legend label is "Unplaced donor" (`R/karyotype-scene.R:93`).
- Meaning: a **donor contig drawn as a neutral dashed bar** because it carries calls with no recipient host, and "Hide unplaced" is off (`R/karyotype.R:1085-1154`).
- **Ruling:** translate as **"неразмещённый донорный контиг"** (or short "донор (не размещён)"). Do **not** translate as a donor-call concept ("неразмещённый донорский вызов"). The spec's literal "unplaced donor" is correct but ambiguous; the glossary pins it to the contig sense.

### 3.2 "deleted recipient" is a recipient-chromosome state, not a "deleted recipient segment"

- Source: a chromosome is "deleted" when its mean recipient breadth is below `deletion_breadth` and below `deletion_relative × median` (`R/karyotype.R:73-109`); it is drawn dashed. Legend "Deleted recipient" (`karyotype-scene.R:91`); no-call card State "deleted recipient (no coverage)" (`mamamiaKaryotypeCard.js:430-433`).
- There is **no literal "deleted recipient segment"** in the source. The spec's `"deleted recipient …"` ellipsis corresponds to these two literals, both about a *chromosome with no coverage*, not a "segment".
- The *call class* "Deletion" is a different concept (a recipient-loss segment; `R/dashboard.R:73-78`, `call-schema.R` `class` = deletion). Distinguish: class "Deletion" → **"Делеция"**; chromosome state "Deleted recipient" → **"Делетированная хромосома реципиента"** / short "Реципиент с делецией".

### 3.3 "component evidence" — the spec mislabels the concept

- The app has **no "component evidence" field**. Two distinct concepts exist:
  1. `class_evidence` (`call-schema.R:19,64,192`): "Evidence that supported the call: donor, recipient, donor+recipient, total imbalance or none." This is *which side drove the class*, not "components".
  2. `component_status` (`call-schema.R:21,67,197`): "Component consistency: consistent or heterogeneous", derived from the spread of donor/recipient residuals across **components** (`R/call-status.R:18-36`).
- **"component" in MaMaMIA = a subgenome / copy-set of a polyploid genome** (e.g., wheat A/B/D, or Th. intermedium J/St/V — see `R/profiles.R:84` "J/St/V resolution" and README "any number of donor and recipient components").
- **Ruling:** translate `class_evidence`-related copy as **"свидетельства класса"** (which side: донор / реципиент / донор+реципиент / нет), and `component` as **"компонент"** (masc.; gloss "субгеном" where it is specifically a subgenome). The least-misleading Russian for the chip "components heterogeneous" is **"компоненты неоднородны"**. Flag to the spec author: the term "component evidence" should be replaced by "class evidence" + "component status" in the spec §1.5 `payload$labels` list, or it will produce a key (`tooltip.component`) whose meaning the translator must not conflate with `tooltip.evidence`.

### 3.4 "Breadth-aware introgression mapping" — breadth ≠ depth

- Source: `R/app-bam.R:1675-1677` parameter help: "Breadth uses covered bases, count uses read depth, and combined uses both." The metric values are `breadth`/`count`/`combined` (`R/app-bam.R:2841-2845`); scoring uses the window metric (`R/scoring.R` consumes `metric`).
- **Breadth = fraction of window covered by ≥1 read** (covered bases / window), a 0–1 quantity; "count" = read depth (reads). Russian "покрытие" alone is ambiguous and usually means depth, so **"с учётом покрытия" overstates/confuses it**. Ruling: header note → **"Картирование интрогрессий по широте покрытия"**; metric label "Breadth" → **"Широта покрытия"**, "Count" → **"Число прочтений"** (gloss "глубина"), "Combined" → **"Комбинированный"**.

### 3.5 Subtitle expansion confirmed

`uvr.toml:4` and `README.md:8` both give "Max Mapping for Meticulous Introgression Analysis". No conflicting expansion exists. Ruling §2.1 applies.

### 3.6 Other source-vs-spec notes (minor)

- The spec §1.5 lists JS literals `"Component"`, `"State"`, `"Length"`, `"Position"`, `"Score"`, `"Chromosome"` (no-call card rows, `mamamiaKaryotypeCard.js:424,427,440,461,465`) and `"Click to filter this group"` (`:472-473`), `"group "` (`:400`), `"No assigned donor segment at this position."` (`:447`). These are real; see `russian-copy.md` `tooltip.*` / `karyotype.*`.
- The spec's `"unplaced donor"` / `"deleted recipient …"` are the two no-call card *titles*; the third title is `"added donor"` and the neutral is `"recipient chromosome"` — both are also user-facing and must be added to the dictionary (the spec omitted them).
- Matched-call tooltip text is fully server-side (`R/call-tooltip.R`), so its labels ("Placement", "Evidence", "Geometry", "Relation", "Assignment", "Group class", "Alternative", "Donor width"/"Call width", "Recipient run", "Overlap", "Width ratio", statuses "confident"/"ambiguous"/"needs_review", chips) are localized in R, not JS — consistent with the spec.
- `exploreIntrogressions()` (`R/dashboard.R:1590`) renders the same Karyotype/Calls tabs plus the full "Tooltip guide" modal and export cards; if it is a live entrypoint it needs the same labels (spec open decision #3).

---

## 4. Terminology rulings most likely to cause rework

1. **"call" → "находка"** (not "вызов"). This changes the Calls tab, the call note, counters, notifications and the export-bundle strings, and any test that greps "Call"/"Calls" in Russian must expect "Находка"/"Находки".
2. **Subtitle stays English** — do not wire the Russian slogan from the spec's illustrative dictionary; the `header.subtitle` ru value must equal the English literal (or be removed so fallback returns English).
3. **"unplaced donor" / "deleted recipient" are chromosome states**, and "component evidence" is a mislabel (§3.1–3.3) — the `payload$labels` keys and their Russian values must be corrected before the engineer generates them.
4. **"breadth" → "широта покрытия"**, not "покрытие"/"глубина" (§3.4).
5. **Number formatting** is R-side, not copy; the "1,000"/"1.23" leaks need a formatter follow-up (§2.4).

---

## 5. Terms that must stay in English/Latin (exhaustive for the in-scope UI)

`MaMaMIA`, `BAM`, `BEDPE`, `FDR`, `GC`, `MAPQ`, `SVG`, `ZIP`, `CSI`, `BAI`, `bamhusk`, `Breadth`/`Count`/`Combined` **only when they are the raw `metric` enum values in logs** (the human labels are translated per §3.4), `paired_log_ratio`/`donor_only`, `max`/`mean`/`sum` (enum values in logs only), `Heatmap`/`Calls` **only when they are the raw `mode` enum** (the radio labels are translated), `EN`/`RU` (switcher), and every BAM field tag (`FLAG`, `RNAME`, `POS`, `CIGAR`, `QUAL`, `SEQ`, `@HD/@SQ/@RG/@PG/@CO`, `SM`, `ID`).

---

## 6. Glossary of classified terms (class / copy state / status / placement)

These are the enumerated `class`/`copy_state`/`class_status`/`recipient_assignment` values rendered by the dashboard (`R/dashboard.R:1370-1404`), the filters (`R/dashboard.R:71-101`) and the tooltip chips. Each is a *value label*, translated once and reused.

| English | Russian | Notes |
|---|---|---|
| Substitution (class) | Замещение | donor replaces recipient segment |
| Addition (class) | Добавление | donor gain, no recipient loss |
| Deletion (class) | Делеция | recipient loss, no donor support |
| Complex (class) | Сложный | conflicting/incompatible spans |
| Other (class) | Прочее | scored but unclassified |
| Balanced (copy state) | Сбалансированное | total depth ~neutral |
| Gain (copy state) | Прирост | total depth up |
| Loss (copy state) | Потеря | total depth down |
| Compound (copy state) | Составное | mixed |
| Neutral (copy state) | Нейтральное | no direction |
| Confident (class status) | Уверенный | — |
| Ambiguous (class status) | Неоднозначный | — |
| Needs review (class status) | Требует проверки | — |
| Assigned (placement) | Размещён | placed on a recipient chr |
| Cross-group (placement) | Межгрупповой | translocation candidate |
| No host (placement) | Нет реципиента | no recipient loss to anchor |
| Group only (placement) | Только группа | fallback placement |
| Passing only (filter) | Только прошедшие | pass calibration |
