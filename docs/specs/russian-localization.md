# MaMaMIA Russian Localization — Implementation Specification

- **Status:** Decision record + spec (no code written).
- **Scope:** Add runtime-switchable Russian (ru) and English (en) UI, defaulting to English.
- **Target:** `bam_launcher_ui()` / `bam_launcher_server()` and the modules they render, plus the
  export-bundle/karyotype text paths described below.
- **Audience:** an implementing engineer who has not read this codebase. Do not guess; everything
  named here is pinned to current line numbers (verified against the working tree).

---

## 0. Goal and non-goals

**Goal:** a per-session language (`en` default) with a header switcher; switching at runtime
localizes every visible label in the main workflow (input, parameters, run, karyotype, calls,
export) **without recomputing the analysis, losing the uploaded BAM, resetting the active tab, or
destroying the karyotype widget**. English output stays byte-identical for the existing test suite
wherever feasible.

**Non-goals (deferred, documented in §6):** localized result bundles (HTML/CSV/README) shipped as
downloads; localized demo-data *content*; Russian plurals; persistence across sessions; browser
`navigator.language` auto-detection.

---

## 1. Decision record

### 1.1 Translation engine — hand-rolled `tr()`, no new dependency

**Chosen:** a hand-rolled `mamamia_tr(key, lang, ...)` backed by nested named lists in a new
`R/i18n.R`. **No `shiny.i18n`.**

**Why not `shiny.i18n`:**

- Adding `shiny.i18n` means editing `uvr.toml`, regenerating `uvr.lock`, and rebuilding the
  deployment image (`deployment/Dockerfile`, `deployment/start.R`). The current lock is already
  consistent; a new runtime dependency is a container-supply-chain change with no offsetting
  benefit for a *static, two-language* dictionary.
- The app is source-based (`mamamia_load()` sources every `R/*.R`, `R/utils.R` first; `R/load.R`
  excluded). A new `R/i18n.R` is picked up automatically with **zero** packaging changes.
- Missing-key detection is a first-class requirement: hand-rolled dictionaries are introspectable
  from R, so we can test key parity en↔ru, empty values, and fallback directly in testthat.

**Dictionary schema** (nested named lists, flattened to dot-paths for lookup):

```r
# R/i18n.R
mamamia_i18n_dicts <- list(
  en = list(
    header = list(
      subtitle  = "Max Mapping for Meticulous Introgression Analysis",
      note      = "Breadth-aware introgression mapping",
      language  = "Language"
    ),
    tabs = list(karyotype = "Karyotype", calls = "Calls", export = "Export"),
    sidebar = list(input = "Input", example_data = "Example data")
    # ... hundreds of leaves
  ),
  ru = list(
    header = list(
      subtitle  = "Максимальное картирование для тщательного анализа интрогрессии",
      note      = "Интрогрессионное картирование с учётом покрытия",
      language  = "Язык"
    ),
    tabs = list(karyotype = "Кариотип", calls = "Вызовы", export = "Экспорт"),
    sidebar = list(input = "Входные данные", example_data = "Пример данных")
    # ... mirrors en exactly (same key tree)
  )
)
```

Russian copy above is an **illustrative placeholder**; final copy is a translator's job and is out
of scope for this spec (see §8 open decisions).

**Key naming convention:** lowercase dot-path leaf keys grouped by UI region —
`header.*`, `tabs.*`, `sidebar.*`, `params.*`, `results.*`, `errors.*`, `tooltip.*`,
`karyotype.*`, `export.*`. Keys are the same in both languages.

**Interpolation:** `{name}` placeholders substituted from named `...` arguments, HTML-escaped by
default. For printf-style strings the *format template* is localized and `sprintf` is still applied
at the call site (the en and ru templates both use `%s`/`%d`):

```r
# {name} form (preferred for new/changed strings):
mamamia_tr("sidebar.example_note", lang, n = 8, recipient = "A", donor = "B")
# printf form (keep existing sprintf call sites):
sprintf(mamamia_tr("errors.windows", lang), n, max)
```

**Fallback behaviour:**

1. `lang` normalized by `mamamia_i18n_valid_lang()`; anything not `"ru"` → `"en"`.
2. Key missing in `ru` → return the **English** value (silent; caught statically by the parity test).
3. Key missing in *both* → `warning("i18n: missing key 'x.y.z'")` and return the key verbatim so
   the rendered page shows the gap and a test can `expect_warning`.

**Public API (all in `R/i18n.R`, `mamamia_`-prefixed to avoid globalenv collisions):**

```r
mamamia_tr(key, lang = "en", ..., .escape = TRUE)   # -> character(1)
mamamia_i18n_langs()                                  # -> c("en", "ru")
mamamia_i18n_valid_lang(lang)                         # -> "en" | "ru"
mamamia_i18n_keys(lang = "en")                        # -> flattened dot-paths (for tests)
mamamia_language_default()                            # -> "en"
```

`mamamia_tr` substitution must HTML-escape values with `htmltools::htmlEscape` unless
`.escape = FALSE` (for values that are themselves already-rendered HTML or trusted).

### 1.2 Runtime switch mechanics — hybrid: server `renderUI` for text-only subtrees + a JS message for labels inside stateful containers

Three options were compared:

- **(a) Whole-page `renderUI`/`uiOutput` keyed on a language reactiveVal.** Rejected: rebuilding the
  page destroys the file input, the `mamamiaKaryotypeCard` htmlwidget, and the active tab. The
  static `state` also carries non-reactive analysis data that a naive re-render would not restore.
- **(b) Client-side JS dictionary + `data-i18n`.** Rejected as the *source of truth*: it duplicates
  the dictionary, bypasses `sprintf`/`paste` server messages, and can't localize server-generated
  `renderUI`/notification text.
- **(c) `htmltools::tagQuery` + `shiny::replaceUI`.** Too invasive across a 5,200-line static UI
  and no better than (a) for state.

**Chosen: hybrid (a-prime) + (b-prime), server as the single source of truth.**

1. **Language state** is a per-session `shiny::reactiveVal("en")` in `bam_launcher_server()`.
2. **Text-only subtrees** are wrapped in `shiny::uiOutput` / `shiny::renderUI` whose expression
   reads `language()` and calls `mamamia_tr()`. Because these subtrees contain **no stateful
   inputs and no widgets**, rebuilding them on a language change is safe. The subtrees are:
   - Header subtitle + header note → new `uiOutput("mamamia_header_text")`.
   - Any currently-static section labels, intro text, help text (sidebar intro, parameter headings,
     call-note paragraph, empty-state copy) that are not already inside a `uiOutput`.
   - **Parameter labels are already `renderUI`** (`output[[paste0("param_label_", id)]]`,
     `output$resolution_param_label`, `output$resolution_ui` — see `R/app-bam.R:4363-4425`). Add
     `language()` as a dependency and wrap the label string in `mamamia_tr()`; no re-architecture.
3. **Labels inside stateful containers that must not be rebuilt** (the `bslib::navset_card_tab`
   titles, which share a container with the DT table and the karyotype widget) are updated by a
   `session$sendCustomMessage("mamamia-set-language", list(lang = language(), labels = ...))`
   handler in the existing inline `<script>` (see 1.7). The handler also sets
   `document.documentElement.lang`.

**What a language change does NOT do:** it never re-reads `input$bams`, never re-runs the analysis
(`result()` is untouched), never rebuilds the file input, the DT table, or the widget element, and
never changes the active tab. Scroll/selection survive because the DOM nodes carrying them are not
replaced.

**Switcher propagation:**

- Switcher input id: **`mamamia_language`**, built once (static) in the header next to the title and
  version badge. `shinyWidgets::radioGroupButtons` is already a dependency, so:

```r
shiny::div(
  class = "mamamia-lang-switch", role = "group",
  `aria-label` = "Language / Язык",   # static bilingual label; buttons show EN/RU
  shinyWidgets::radioGroupButtons(
    inputId = "mamamia_language",
    label = NULL,
    choices = c("EN" = "en", "RU" = "ru"),
    selected = "en", size = "sm", individual = TRUE, direction = "horizontal"
  )
)
```

- Server observer (in `bam_launcher_server`, next to the existing `observeEvent`s):

```r
language <- shiny::reactiveVal("en")
shiny::observeEvent(input$mamamia_language, {
    language(mamamia_i18n_valid_lang(input$mamamia_language))
}, ignoreNULL = TRUE, ignoreInit = TRUE)

shiny::observeEvent(language(), {
    session$sendCustomMessage("mamamia-set-language", list(
        lang = language(),
        tabs = list(
            karyotype = mamamia_tr("tabs.karyotype", language()),
            calls     = mamamia_tr("tabs.calls",     language()),
            export    = mamamia_tr("tabs.export",    language())
        )
    ))
})
```

Initial value is `"en"`; `selected = "en"` matches. No reliance on `navigator.language`.

### 1.3 Where language state lives and its lifetime

- Per-session `shiny::reactiveVal("en")` created at the top of `bam_launcher_server()` (after
  `session_dir`, before `result <- reactiveVal(NULL)` at `R/app-bam.R:3641`).
- Lifetime: the session. **Persistence across sessions is out of scope** (would need a cookie or
  `bookmarking`; not worth it here).
- Default is hard-coded `"en"` because the test suite asserts English output (see §5). Browser
  locale is deliberately **not** consulted.

### 1.4 Translation engine placement

`R/i18n.R` is sourced automatically by `mamamia_load()` (alphabetical after `R/utils.R`; it does
not need to precede `R/app-bam.R` because the dictionaries are data and `mamamia_tr` is only
*called* at UI/server build time, which happens after all files are sourced). No change to
`R/load.R`, `uvr.toml`, `uvr.lock`, or the deployment image is required.

### 1.5 Language-coupled front-end logic must become language-independent

Two places compare/derive behaviour from English text and will break under Russian:

1. **Tab enable/disable JS** (`R/app-bam.R:3047-3073`). `mamamiaResultTabs()` matches
   `link.textContent.trim()` against `'Calls'`, `'Export'`, `'Karyotype'`. Fix by giving every
   `nav_panel` an explicit `value` and matching on the `data-value` attribute instead of text:
   - In `bam_results_ui()` (`R/app-bam.R:2943-2977`): add `value = "karyotype"`, `"calls"`,
     `"export"` to the three `nav_panel(...)` calls.
   - In `explorer_results_ui()` (`R/dashboard.R:1038-1057`): add `value = "karyotype"`, `"calls"`.
   - Rewrite `mamamiaResultTabs` to select `.nav-link` by `closest('[data-value]')` and compare
     `data-value` to `"calls"`/`"export"` (enable/disable) and `"karyotype"` (active-tab fallback).
     Keep the `#explorer-tabs` selector.

2. **Karyotype widget** (`www/htmlwidgets/mamamiaKaryotypeCard.js`). Hardcoded user-facing strings:
   `aria-label` `"Karyotype in assembly order"` (line ~751), empty-state fallback
   `"No chromosomes in this result."` (~661), legacy/no-call tooltip literals `"call"`, `"no call"`,
   `"Component"`, `"State"`, `"Length"`, `"Position"`, `"Score"`, `"Chromosome"`, `"Click to filter
   this group"`, `"group "`, `"unplaced donor"`, `"deleted recipient …"`, `"No assigned donor
   segment at this position."`. Add a server-resolved `labels` block to the payload and have the JS
   prefer it, **keeping the English literals as fallback defaults** so
   `tests/testthat/test-karyotype-scene.R:637` (`grepl("unplaced donor", js)`) and
   `test-call-tooltip.R` (asserts the matched-call card only reads `matched.tooltip`) stay green:

```r
# in karyotype_payload() (R/karyotype.R) — add:
payload$labels <- list(
  aria_label = mamamia_tr("karyotype.aria_label", lang),
  empty      = mamamia_tr("karyotype.empty", lang),
  call       = mamamia_tr("tooltip.call", lang),
  no_call    = mamamia_tr("tooltip.no_call", lang),
  component  = mamamia_tr("tooltip.component", lang),
  state      = mamamia_tr("tooltip.state", lang),
  length     = mamamia_tr("tooltip.length", lang),
  position   = mamamia_tr("tooltip.position", lang),
  score      = mamamia_tr("tooltip.score", lang),
  chromosome = mamamia_tr("tooltip.chromosome", lang),
  click_hint = mamamia_tr("tooltip.click_to_filter", lang)
  # ... and every literal the JS currently emits
)
```

   JS changes: `var L = (payload && payload.labels) || {};` then `L.aria_label || "Karyotype in
   assembly order"`, etc. **Matched-call tooltip text is already fully server-side** (the tooltip
   model built in `R/call-tooltip.R` / `R/karyotype.R`), so its localization is pure R — no JS
   string surgery for that path.

### 1.6 `<html lang>` attribute

`bslib::page_fluid` defaults the document to `lang="en"`. Add `lang = "en"` explicitly to the
`page_fluid(...)` call in `bam_launcher_ui()` (`R/app-bam.R:3026`) and have the
`mamamia-set-language` handler set `document.documentElement.lang = x.lang`. This is an
accessibility correctness fix (screen readers announce the language) and costs one JS line.

### 1.7 New JS custom message handler

Append to the existing inline `<script>` in `bam_launcher_ui()` (`R/app-bam.R:3033-3490`) a handler:

```js
Shiny.addCustomMessageHandler('mamamia-set-language', function (x) {
  if (x && x.lang) { document.documentElement.lang = x.lang; }
  if (x && x.tabs) {
    var map = { karyotype: x.tabs.karyotype, calls: x.tabs.calls, export: x.tabs.export };
    document.querySelectorAll('#explorer-tabs .nav-link').forEach(function (link) {
      var item = link.closest('[data-value]');
      var value = item ? item.getAttribute('data-value') : null;
      if (value && map[value]) {
        // preserve the icon/HTML structure: replace only the trailing text node
        link.childNodes.forEach(function (n) {
          if (n.nodeType === Node.TEXT_NODE && n.textContent.trim() !== '') {
            n.textContent = map[value];
          }
        });
      }
    });
  }
});
```

(Exact text-node update strategy is an implementation detail; the contract is: **never change tab
labels by rebuilding the nav, never match on English text, keep `data-value` as the key**.)

---

## 2. Scope matrix

| Area | In scope (Phase 1) | Rationale |
|---|---|---|
| Header chrome (subtitle, note, switcher label) | **Yes** | First thing a user sees; cheap (`uiOutput` + `tr`) |
| Sidebar labels, intro, help, buttons, titles | **Yes** | Core workflow entry |
| Parameter panel labels/headings/legend (Resolution, Scoring, Masking, Calibration, `param_label_*`, `resolution_param_label`, `resolution_ui`) | **Yes** | Already `renderUI`-based; `tr()` is a one-line change per label |
| Run / progress stage labels, log messages, notifications, error text (`R/app-bam.R`, `R/preflight.R`) | **Yes** | Server messages are `sprintf`/`paste`; wrap templates in `tr()` |
| Results UI (tabs, call note, filter labels, selection controls, DT column headers via `bam_call_table_header`) | **Yes** | Core workflow; tab titles via JS handler, rest via `tr()` |
| Karyotype widget legend, empty-state, aria-label, no-call/legacy tooltip | **Yes** | Via `payload$labels` (§1.5) + `tr()` in `R/karyotype.R`/`R/call-tooltip.R` |
| `<html lang>` update + switcher aria | **Yes** | Accessibility |
| Plot/axis/legend text generated by **server-side R** (`ggplot2` titles in `R/export-bundle-plots.R`, `R/dashboard.R`) | **Yes for in-app plots** | In-app chart text is user-facing; localized via `tr()` at build time |
| Result bundles (HTML/CSV/README from `R/export-bundle*.R`, `R/cohort-export.R`) | **Deferred (documented follow-up)** | High blast radius: exported artifacts are byte-asserted by `test-export-bundle.R` (`"Calls by class"`), `test-dashboard.R` (`"Calls that pass calibration filters."`), `test-example-data.R` (`"Example dataset"`); they are a *different product surface* (a report the user downloads), and shipping untranslated-but-consistent English bundles is acceptable in Phase 1. See §8 for the explicit user decision. |
| Demo-data *content* (the bundled example dataset's biological labels) | **Deferred** | The demo data is scientific content, not UI copy; translating its inline labels would change fixtures and downstream assertions. |
| Russian plurals (1/2-4/5+ categories) | **Deferred — simple forms accepted** | See §6. Strings are written to avoid bare number↔noun agreement where possible. |

**Phased extraction order** (each phase is independently shippable):

1. **P0 — engine + switcher + header + tabs** (the mechanism, `R/i18n.R`, switcher, `mamamia-set-language`, `data-value` fix, `<html lang>`).
2. **P1 — sidebar + parameter panel + run/progress/errors** (input→parameters→run path).
3. **P2 — results UI + karyotype/tooltip labels + in-app plot text** (calls/export/karyotype path).
4. **P3 — (follow-up) exported bundles + demo content + plurals.**

The header switcher and **every visible label in the main workflow** must be covered by end of P2.

---

## 3. Cyrillic font

`www/mamamia-font.css` declares a single **Latin** `@font-face` for Inter
(`www/fonts/inter-latin.woff2`, unicode-range U+0000-00FF + Latin-1 supplement etc.). Russian text
would silently fall back to a system font. Required:

1. Add `www/fonts/inter-cyrillic.woff2` (Inter Cyrillic subset, SIL OFL 1.1 — the same license the
   repo already ships as `www/fonts/OFL.txt`; the OFL permits redistribution and modification and
   already covers the family, so no new license file is strictly required, but the new file's
   provenance must be recorded in a `www/fonts/README.md` noting the upstream Inter release and
   version).
2. Add a second `@font-face` in `www/mamamia-font.css` with `font-family: "Inter"`, the same
   `font-weight: 100 900`, and a Cyrillic `unicode-range` (`U+0400-045F, U+0490-0491, U+04B0-04B1,
   U+2116, U+2010-2027` — Cyrillic block + Ё/ё and common punctuation) pointing at the new file.
   Keeping two `@font-face` blocks with disjoint `unicode-range` makes the browser download only the
   subset it needs and preserves the existing Latin rendering byte-for-byte.
3. Add an explicit Cyrillic-capable fallback stack to `bam_dashboard_theme()`'s `base_font`
   (`R/app-theme.R:60`) and to the widget's `plotFont` in `mamamiaKaryotypeCard.js:7` — e.g.
   `"Inter", "Inter Variable", system-ui, -apple-system, "Segoe UI", Roboto, "Noto Sans",
   "Helvetica Neue", Arial, sans-serif` — so any not-yet-covered glyph degrades to a font that has
   Cyrillic rather than a tofu box.

**Assumption (not verifiable here):** the exact upstream Inter release/version and the checksum of
`inter-cyrillic.woff2` must be pinned by the implementer at download time (record version + sha256
in `www/fonts/README.md`). The Latin subset in-tree is 48,256 bytes; the Cyrillic subset should be
sourced from the same Inter release so the two subsets are metric-compatible.

---

## 4. Migration strategy for ~hundreds of inline strings

Goal: extract strings to `R/i18n.R` **without changing English output**, in commits ordered so the
suite stays green at every step.

1. **Commit 1 — engine only, no call sites changed.** Add `R/i18n.R` (dictionaries seeded with a
   handful of keys), add `tests/testthat/test-i18n.R` (dictionary tests). English strings are
   copied **verbatim** from the current literals into `mamamia_i18n_dicts$en`, so `tr()` returns
   the identical bytes.
2. **Commit 2 — mechanism wiring, no copy change.** Add `language()` reactiveVal, the switcher,
   `mamamia-set-language` handler, `<html lang>` update, and the `data-value` attribute + JS rewrite
   (§1.5). At this point `language` is always `"en"`, so no visible change; the JS change must be
   behaviour-preserving (verify `mamamiaResultTabs` still enables/disables the same tabs).
3. **Commit 3+ — extract region by region** (header, sidebar, params, run/errors, results, karyotype).
   Each extraction is: copy the literal into `en`, translate into `ru`, replace the literal at the
   call site with `mamamia_tr(...)`. English stays byte-identical; Russian keys are added in the
   same commit as their call sites (so the parity test never fails in CI).
4. **Rule for `sprintf`/`paste` call sites:** the localized *format template* goes through
   `mamamia_tr(key, lang)` and `sprintf`/`paste` is applied after. The English template is the
   current literal; the Russian template keeps the same `%s`/`%d` placeholders in the same order.
   For `paste(...)` used only to wrap long lines, collapse into a single template string first.
5. **Which tests must change and why:** none of the four named assertions may change their
   *expected string* if the strings stay English by default. `test-karyotype.R:505`
   (`expect_match(html, "Karyotype")`) and `:509` (`"Hide deletions"`) exercise
   `explorer_results_ui()`/`explorer_karyotype_ui()` with the default English language, so they
   keep passing **iff** the default language is `en` and `tr()` is pure (no side effects at build
   time). The only tests that change in this project are: (a) new `test-i18n.R`; (b) a new JS
   regression test asserting `mamamiaResultTabs`/the handler keys on `data-value` (§5); (c) **no**
   existing assertion should need its literal updated unless a string is *deliberately* changed —
   flag any such case in the PR rather than silently editing a test.

---

## 5. Testing plan

New test files (all under `tests/testthat/`, following the existing `testthat` >= 3.0 convention):

### `tests/testthat/test-i18n.R` — dictionary unit tests

- `test_that("en and ru dictionaries have identical key sets")` — `expect_identical(
  sort(mamamia_i18n_keys("en")), sort(mamamia_i18n_keys("ru")) )`.
- `test_that("no dictionary value is empty")` — for every leaf in both languages,
  `expect_true(nzchar(value))`.
- `test_that("mamamia_tr falls back to English for a missing Russian key")` — assert
  `mamamia_tr("header.subtitle", "ru")` equals `mamamia_tr("header.subtitle", "en")` when the ru
  leaf is absent, **and** that the static parity test is what surfaces it (i.e. a ru key never
  silently differs).
- `test_that("a key missing in both languages warns and returns the key")` — `expect_warning(
  x <- mamamia_tr("no.such.key", "en") )` and `expect_identical(x, "no.such.key")`.
- `test_that("mamamia_tr interpolates and HTML-escapes")` — `mamamia_tr("sidebar.example_note",
  "en", n = 8)` substitutes `{n}`; a value containing `<` is escaped when `.escape = TRUE` and left
  raw when `.escape = FALSE`.
- `test_that("invalid language normalizes to en")` — `mamamia_i18n_valid_lang("fr")` → `"en"`.
- `test_that("default language is en")` — `expect_identical(mamamia_language_default(), "en")`.

### `tests/testthat/test-i18n-shiny.R` — behaviour (uses `shiny::testServer`)

- `test_that("toggling the language input changes a label and preserves analysis state")` — start
  `bam_launcher_server` with a valid `state` (reuse `helper-load.R`/`helper-bam.R` fixtures),
  set `session$setInputs(mamamia_language = "ru")`, then assert (i) a representative label rendered
  via `tr()` (e.g. the resolved `param_label_metric` output) changed to the Russian string, and (ii)
  the `result()` reactive is unchanged (compare identity before/after — no recompute).
- `test_that("language defaults to en and observer is inert on init")` — with no input,
  `expect_identical(language(), "en")` and no custom message was sent (spy on
  `session$sendCustomMessage`).

### `tests/testthat/test-i18n-js.R` — regression for the tab-enable logic

- `test_that("tab enable/disable keys on data-value, not English text")` — read
  `R/app-bam.R` (or the rendered `<script>`) and assert it **does not** compare against the
  literals `'Calls'`, `'Export'`, `'Karyotype'` for enable/disable, and **does** reference
  `data-value`. (Mirrors the style of `test-karyotype-scene.R:630` which greps the widget JS.)
- `test_that("nav panels expose explicit values")` — `as.character(bam_results_ui())` contains
  `data-value="karyotype"`, `data-value="calls"`, `data-value="export"`.

### Existing tests that must stay green (regression guard)

- `test-example-data.R:92` (`"Example dataset"`), `test-dashboard.R:831/844/876`
  (`"Calls that pass calibration filters."`), `test-export-bundle.R:543` (`"Calls by class"`),
  `test-karyotype.R:505` (`"Karyotype"`), `:509/512/518` (`"Hide deletions"`, `"Hide unplaced"`,
  `"Tooltip help"`) — all default-English, must pass with **no literal edits**.
- `test-karyotype-scene.R:637` (`grepl("unplaced donor", js)`) — preserved by keeping JS literals
  as fallback.

**Run command:** `Rscript -e 'testthat::test_dir("tests/testthat")'` (or the repo's existing test
runner); the acceptance gate is the full suite plus the three new files.

---

## 6. Risks / edge cases

- **Cyrillic font** — without the subset (§3) the UI renders inconsistently. Mitigated by the
  second `@font-face` + fallback stack; visual check required.
- **Text overflow in fixed header / narrow sidebar** — Russian strings are longer than English
  (e.g. the subtitle, `header-note`, sidebar buttons). `.mamamia-header` is 150px-tall grid row and
  `.mamamia-header-note` uses `white-space: nowrap` + uppercase. Add CSS: allow `white-space:
  normal` / smaller letter-spacing for `ru` (a `.lang-ru` class on `<html>` or body toggled by the
  JS handler), and verify the sidebar at narrow widths. This is the main visual risk.
- **`sprintf`/plural agreement** — Russian has 3 plural categories. **Decision: accept simple
  forms in Phase 1**; write Russian copy to avoid bare `число + существительное` agreement (use
  "Окон: N" style) or accept the imprecision. No `ngettext`-style helper yet; if a string genuinely
  needs it later, add `mamamia_tr_plural(key, lang, n, ...)` behind the same dictionary.
- **HTML escaping in interpolations** — `mamamia_tr` escapes `...` values by default; call sites
  that interpolate already-safe HTML must pass `.escape = FALSE`. Sample names from user uploads
  are the highest-risk input (a `<` in a filename must never become a tag).
- **`lang` attribute** — must update both `<html lang>` (JS) and any per-element `lang` on
  server-rendered report fragments (`R/export-bundle*.R` hardcode `<html lang="en">`; deferred, but
  note it so the follow-up changes it consistently).
- **Switcher accessibility** — `radioGroupButtons` renders real buttons; provide a
  `role="group"` + static `aria-label` (bilingual, short) and rely on EN/RU button labels being
  language-neutral. Screen-reader users get the language announced via `<html lang>`.
- **Karyotype re-render on language change** — the payload strings are resolved by
  `karyotype_payload()`; the karyotype `renderX` output must take `language()` as a dependency so it
  re-renders (calling the widget's `renderValue` **in place**, not rebuilding the element) with the
  new labels, and must re-read `selected_group()` so selection is preserved. Re-rendering must not
  re-run the analysis.
- **Two tab UIs** — `bam_results_ui()` (launcher) and `explorer_results_ui()` (`R/dashboard.R`,
  used by `exploreIntrogressions()`) both render Karyotype/Calls tabs; both need `data-value`, and
  the JS selector `#explorer-tabs` only matches the launcher's nav. If `exploreIntrogressions()` is
  a live entrypoint, it needs the same label handling; otherwise note it as out of scope.

---

## 7. Files to create / modify

### Create

| File | Purpose |
|---|---|
| `R/i18n.R` | Translation engine: `mamamia_tr()`, dictionaries `mamamia_i18n_dicts`, key/validation helpers. |
| `www/fonts/inter-cyrillic.woff2` | Inter Cyrillic subset (OFL 1.1); downloaded, version + sha256 recorded. |
| `www/fonts/README.md` | Provenance of the new subset (upstream Inter release, version, checksum). |
| `tests/testthat/test-i18n.R` | Dictionary unit tests (parity, emptiness, fallback, interpolation, escaping). |
| `tests/testthat/test-i18n-shiny.R` | `testServer` tests for the switcher + state preservation. |
| `tests/testthat/test-i18n-js.R` | Regression test for `data-value`-keyed tab logic + nav values. |

### Modify

| File | Change |
|---|---|
| `R/app-bam.R` | Switcher in header (~3497-3540); `uiOutput("mamamia_header_text")` for subtitle/note; `page_fluid(..., lang="en")`; JS handler + `data-value` fix (§1.5); wrap `tr()` around every extracted string (sidebar, params, run/errors, results); `language()` reactiveVal + observers in `bam_launcher_server()`; make karyotype output depend on `language()`. |
| `R/app-theme.R` | Add Cyrillic-capable fallback to `base_font`; add `.lang-ru` CSS + header/sidebar overflow rules (§6). |
| `R/dashboard.R` | `value=` on `explorer_results_ui()` nav panels; wrap dashboard labels/table headers/plot titles in `tr()`. |
| `R/karyotype.R`, `R/karyotype-scene.R` | Add `payload$labels` block; wrap tooltip/legend/empty-state strings in `tr()`. |
| `R/call-tooltip.R` | Wrap tooltip model labels in `tr()`. |
| `R/example-data.R` | Wrap `bam_example_note()` text in `tr()` (keep English byte-identical). |
| `R/preflight.R` | Wrap preflight/log/error message templates in `tr()`. |
| `R/export-bundle*.R`, `R/cohort-export.R` | **Deferred** (see §2); do not touch in Phase 1. |
| `www/mamamia-font.css` | Add Cyrillic `@font-face` with correct `unicode-range`. |
| `www/htmlwidgets/mamamiaKaryotypeCard.js` | Prefer `payload.labels`, keep English literals as fallback; add Cyrillic to `plotFont`. |

---

## 8. Ordered implementation steps

1. **Engine.** Create `R/i18n.R` with the API in §1.1; seed `mamamia_i18n_dicts$en` with a minimal
   key set. Write `test-i18n.R`; run it green.
2. **Mechanism.** Add `language()` + switcher + `mamamia-set-language` handler + `<html lang>`
   update + `data-value` on nav panels + JS rewrite. Write `test-i18n-js.R`. Run full suite green
   (still English). Add Cyrillic font files + `@font-face` + fallback stack.
3. **Header + tabs.** Extract header strings; wire `mamamia_header_text`; localize tab titles via
   the handler. Manual visual check of RU header/overflow.
4. **Sidebar + parameters + run/errors.** Extract in commit-sized regions; labels already in
   `renderUI` get `language()` as a dependency + `tr()`. Add `test-i18n-shiny.R` for the switcher.
5. **Results + karyotype.** `payload$labels` in `R/karyotype.R`; `tr()` in tooltip/scene strings;
   karyotype output depends on `language()` and re-reads `selected_group()`. Wrap dashboard labels.
6. **Regression + acceptance.** Run the full suite (including the four named English assertions)
   and the manual RU pass; record screenshots; open the deferred-bundle follow-up.

---

## 9. Acceptance criteria (reviewer checklist)

- App boots in `en`; existing tests pass with **no English literal edits** (the four named
  assertions plus `test-karyotype-scene.R:637`).
- Header switcher (id `mamamia_language`) appears next to the title/version badge; clicking `RU`
  localizes header, sidebar, parameter labels, tabs, run/progress/errors, and karyotype
  tooltip/legend/empty-state **without** re-uploading the BAM, recomputing `result()`, resetting the
  active tab, or losing scroll/selection.
- `<html lang>` flips `en`↔`ru`; switcher has an accessible name.
- `R CMD`/testthat: `test-i18n.R`, `test-i18n-shiny.R`, `test-i18n-js.R` pass.
- Russian renders in Inter Cyrillic (no tofu); header/sidebar do not overflow at 980px and below.
- `data-value`-based tab logic: after switching to RU, the Calls/Export tabs are still correctly
  disabled until a result exists, and the Karyotype tab is still the active-tab fallback.
- No changes to `uvr.toml`, `uvr.lock`, `deployment/*`, or `R/load.R`.

---

## 10. Open decisions (user must resolve)

1. **Result bundles** (HTML/CSV/README from `R/export-bundle*.R`, `R/cohort-export.R`): ship
   Russian bundles in Phase 1, or defer (recommended) and file as follow-up? This is a product
   decision, not a technical one.
2. **Final Russian copy**: who owns/approves translations? The mechanism is spec'd here; the words
   are not.
3. **`exploreIntrogressions()`** (`R/dashboard.R`) — is it a live entrypoint that must localize, or
   a legacy/internal path we can leave English?
4. **Plurals**: accept simple forms (recommended for Phase 1) or invest in `mamamia_tr_plural` now?
5. **Persistence**: is per-session-only acceptable, or must the choice survive a reload (cookie)?

---

## 11. Weakest link

**The karyotype/tooltip localization path** (`R/karyotype.R` → `payload$labels` →
`mamamiaKaryotypeCard.js`), specifically the requirement that a language change re-render the
widget **in place** via `renderValue` without rebuilding the element and while preserving
`selected_group()`. It is the only place where localization and a stateful interactive widget
overlap, where the English fallback literals must survive (to keep `test-karyotype-scene.R:637`
green), and where the two-tab-UI ambiguity (§6) can silently produce a half-localized screen. Every
other decision degrades gracefully; this one fails visibly.
