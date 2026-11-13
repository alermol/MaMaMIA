## Dictionary and translation-engine contracts.

test_that("en and ru dictionaries have identical key sets", {
    expect_identical(
        sort(mamamia_i18n_keys("en"), method = "radix"),
        sort(mamamia_i18n_keys("ru"), method = "radix")
    )
})

test_that("no dictionary value is empty", {
    for (lang in mamamia_i18n_langs()) {
        values <- mamamia_i18n_dicts[[lang]]
        expect_true(length(values) > 0L)
        expect_true(all(nzchar(values)))
        expect_false(anyNA(values))
    }
})

test_that("the full tooltip-guide block is present", {
    keys <- mamamia_i18n_keys("en")
    expect_identical(sum(grepl("^tooltip_guide\\.", keys)), 137L)
    expect_true("tooltip_guide.caveats.heatmap_score.desc" %in% keys)
})

test_that("the subtitle stays English in both languages", {
    expect_identical(
        mamamia_tr("header.subtitle", "ru"),
        mamamia_tr("header.subtitle", "en")
    )
})

test_that("mamamia_tr falls back to English for a missing Russian key", {
    old <- mamamia_i18n_dicts
    on.exit(assign("mamamia_i18n_dicts", old, envir = globalenv()))
    broken <- old
    broken$ru <- as.list(broken$ru)
    broken$ru[["header.note"]] <- NULL
    assign("mamamia_i18n_dicts", broken, envir = globalenv())
    expect_identical(
        mamamia_tr("header.note", "ru"),
        mamamia_tr("header.note", "en")
    )
})

test_that("a key missing in both languages warns and returns the key", {
    expect_warning(
        value <- mamamia_tr("no.such.key", "en"),
        "missing key"
    )
    expect_identical(value, "no.such.key")
})

test_that("mamamia_tr interpolates placeholders", {
    expect_identical(
        mamamia_tr("sidebar.example_note", "en", n = 8, donor = "A"),
        "Example 8: A donor and {recipient} recipient."
    )
})

test_that("escaping policy escapes interpolated values exactly once", {
    value <- "A&B <C>"

    ## Plain text result: the value is not pre-escaped; the text sink
    ## (htmltools) performs the single escaping pass.
    plain <- mamamia_tr(
        "sidebar.example_note", "en",
        n = 1, donor = value, recipient = "R"
    )
    expect_match(plain, value, fixed = TRUE)
    rendered <- as.character(shiny::p(plain))
    expect_match(rendered, "A&amp;B &lt;C&gt;", fixed = TRUE)
    expect_false(grepl("&amp;amp;", rendered, fixed = TRUE))

    ## Raw-HTML sink: the value is escaped once by mamamia_tr, and the
    ## resulting string is handed to HTML() without a second pass.
    escaped <- mamamia_tr(
        "sidebar.example_note", "en",
        n = 1, donor = value, recipient = "R", .escape = TRUE
    )
    expect_match(escaped, "A&amp;B &lt;C&gt;", fixed = TRUE)
    expect_false(grepl("&amp;amp;", escaped, fixed = TRUE))
    raw <- as.character(shiny::HTML(escaped))
    expect_match(raw, "A&amp;B &lt;C&gt;", fixed = TRUE)
})

test_that("invalid language normalizes to en", {
    expect_identical(mamamia_i18n_valid_lang("fr"), "en")
    expect_identical(mamamia_i18n_valid_lang(NA_character_), "en")
    expect_identical(mamamia_i18n_valid_lang("RU"), "ru")
})

test_that("default language is en", {
    expect_identical(mamamia_language_default(), "en")
})
test_that("Russian error file labels contain no English File token", {
    one <- bam_error_file_label("sample.bam", lang = "ru")
    expect_identical(one, "Файл: sample.bam")
    expect_false(grepl("File", one, fixed = TRUE))

    many <- bam_error_file_label(
        paste0("sample", 1:7, ".bam"), lang = "ru"
    )
    expect_match(many, "^Файлы: ", perl = TRUE)
    expect_false(grepl("Files", many, fixed = TRUE))
    expect_match(many, "(ещё 2)", fixed = TRUE)
})

test_that("every call-table label/description key exists and is reachable", {
    keys_en <- mamamia_i18n_keys("en")
    keys_ru <- mamamia_i18n_keys("ru")
    columns <- bam_call_column_metadata()$column
    expect_length(columns, 64L)

    for (column in columns) {
        label_key <- paste0("call_table.label.", column)
        desc_key <- paste0("call_table.desc.", column)
        expect_true(label_key %in% keys_en)
        expect_true(label_key %in% keys_ru)
        expect_true(desc_key %in% keys_en)
        expect_true(desc_key %in% keys_ru)
    }

    header <- as.character(
        bam_call_table_header(c("group", "class", "score_mean"), "ru")
    )
    expect_match(header, "Группа", fixed = TRUE)
    expect_match(header, "Класс", fixed = TRUE)
    expect_match(header, "Оценка", fixed = TRUE)
    expect_false(grepl("Chromosome", header, fixed = TRUE))
    expect_false(grepl(">Class<", header, fixed = TRUE))
    expect_false(grepl(">Score<", header, fixed = TRUE))
})

test_that("oversized-upload and upload-help produce Russian copy", {
    modal <- as.character(bam_upload_too_large_modal(
        list(size = 1024, limit = 512, files = "sample.bam"),
        "ru"
    ))
    expect_match(modal, "Этот BAM имеет размер", fixed = TRUE)
    expect_match(modal, "данные BAM не обрабатывались", fixed = TRUE)
    expect_false(grepl("This BAM is", modal, fixed = TRUE))

    help <- as.character(bam_upload_help_ui(1024, "ru"))
    expect_match(help, "Всё, что читает MaMaMIA", fixed = TRUE)
    expect_match(help, "имена прочтений", fixed = TRUE)
    expect_match(help, "Длина прочтений", fixed = TRUE)
    expect_false(grepl("Everything MaMaMIA reads", help, fixed = TRUE))
})
test_that("the upload-help modal footer is localized", {
    en <- as.character(bam_upload_help_modal(1024, "en"))
    ru <- as.character(bam_upload_help_modal(1024, "ru"))
    expect_match(en, ">Close<", fixed = TRUE)
    expect_match(ru, "Закрыть", fixed = TRUE)
    expect_false(grepl(">Close<", ru, fixed = TRUE))
})

test_that("the BEDPE selection counter has no leading ellipsis", {
    expect_identical(
        paste(3, mamamia_tr("results.selection_status", "en")),
        "3 call(s) selected for BEDPE export"
    )
    expect_false(startsWith(
        mamamia_tr("results.selection_status", "en"), "\u2026"
    ))
})

test_that("cached-result messages are localized without English residue", {
    for (key in c("log.restored_cached", "log.using_cached_notify")) {
        ru <- mamamia_tr(key, "ru")
        expect_false(grepl("cached analysis result", ru, fixed = TRUE))
        expect_true(nzchar(ru))
    }
    expect_identical(
        mamamia_tr("log.restored_cached", "en"),
        "Restored a cached analysis result for these inputs; no recalculation was needed."
    )
    expect_identical(
        mamamia_tr("log.using_cached_notify", "en"),
        "Using the cached analysis result for these inputs."
    )
})
