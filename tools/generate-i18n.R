#!/usr/bin/env Rscript

## One-off generator for R/i18n.R.
##
## Parses the EN/RU tables in docs/specs/russian-copy.md and writes the
## dictionary as plain data. The generated file is committed; this script is
## not used at runtime. Re-run it only when the copy file changes:
##
##   Rscript tools/generate-i18n.R
##
## The parser is deliberately narrow: it reads the primary tables only (the
## plural-avoidance rewrite table repeats keys that already appear above it)
## and applies the few key-level corrections documented in
## docs/specs/russian-glossary.md and the task brief.

root <- normalizePath(
    Sys.getenv("MAMAMIA_ROOT", unset = getwd()),
    winslash = "/", mustWork = TRUE
)
source_path <- file.path(root, "docs", "specs", "russian-copy.md")
target_path <- Sys.getenv(
    "MAMAMIA_I18N_OUT",
    unset = file.path(root, "R", "i18n.R")
)

## Keep generation byte-identical when the caller runs under LC_ALL=C.
if (!grepl("UTF-8|utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
    for (candidate in c("C.UTF-8", "C.utf8", "en_US.UTF-8", "en_US.utf8")) {
        if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", candidate)))) {
            break
        }
    }
}

lines <- readLines(source_path, warn = FALSE, encoding = "UTF-8")

split_row <- function(line) {
    body <- sub("^\\|", "", sub("\\|\\s*$", "", line))
    cells <- strsplit(body, "(?<!\\\\)\\|", perl = TRUE)[[1L]]
    trimws(gsub("\\\\\\|", "|", cells))
}

## Some copy cells quote a value to make leading/trailing spaces visible
## (" from "). Drop only those outer quotes; the spaces are the copy.
strip_visible_quotes <- function(value) {
    if (grepl('^".*"$', value)) {
        substr(value, 2L, nchar(value) - 1L)
    } else {
        value
    }
}

is_key <- function(x) {
    grepl("^[a-z][a-z0-9_]*(\\.[a-z0-9_]+)+$", x)
}

rows <- list()
for (line_number in seq_along(lines)) {
    line <- lines[[line_number]]
    if (!startsWith(trimws(line), "|")) {
        next
    }
    cells <- split_row(trimws(line))
    if (length(cells) < 3L) {
        next
    }
    key <- gsub("`", "", cells[[1L]])
    if (!is_key(key)) {
        next
    }
    en <- strip_visible_quotes(cells[[2L]])
    ru <- strip_visible_quotes(cells[[3L]])
    ## Strip translator notes that follow the copy in the RU cell.
    ru <- sub("\\s*\\*\\(.*\\)\\*\\s*$", "", ru)
    if (key %in% names(rows)) {
        previous <- rows[[key]]
        if (!identical(previous$en, en) || !identical(previous$ru, ru)) {
            stop(
                "Conflicting duplicate dictionary key '", key,
                "' at copy line ", line_number,
                ": '", previous$en, "' / '", en, "'",
                call. = FALSE
            )
        }
    } else {
        rows[[key]] <- list(en = en, ru = ru)
    }
}

## Keys whose RU cell is not usable copy.
rows[["header.language"]] <- list(en = "Language", ru = "Язык")
rows[["errors.short_contigs"]] <- list(
    en = "Contig %s is shorter than the requested %s windows; choose a lower resolution.",
    ru = "Контиг %s короче запрошенных %s окон; выберите меньшее разрешение."
)
rows[["errors.short_contigs_plural"]] <- list(
    en = "Contigs %s are shorter than the requested %s windows; choose a lower resolution.",
    ru = "Контиги %s короче запрошенных %s окон; выберите меньшее разрешение."
)

## The copy file flags these as missing from results.* (glossary note in
## russian-copy.md, "Coverage note"): confidence choices and the enum value
## labels used by the filter pickers and the DT table.
extras <- list(
    "results.unplaced_note" = c(
        paste(
            "Calls that cannot be placed on a recipient chromosome or an",
            "added donor bar appear as neutral unplaced donor bars when",
            "Hide unplaced is off. The table lists only calls placed in the",
            "current karyotype."
        ),
        paste(
            "Находки, которые нельзя разместить на хромосоме реципиента или",
            "на полосе добавленного донора, отображаются как нейтральные",
            "полосы неразмещённых доноров, когда «Скрыть неразмещённые»",
            "выключено. Таблица содержит только находки, размещённые в",
            "текущем кариотипе."
        )
    ),
    "results.confidence_high" = c("High", "Высокая"),
    "results.confidence_moderate" = c("Moderate", "Средняя"),
    "results.confidence_low" = c("Low", "Низкая"),
    "results.class_substitution" = c("Substitution", "Замещение"),
    "results.class_addition" = c("Addition", "Добавление"),
    "results.class_deletion" = c("Deletion", "Делеция"),
    "results.class_complex" = c("Complex", "Сложный"),
    "results.class_other" = c("Other", "Прочее"),
    "results.status_confident" = c("Confident", "Уверенный"),
    "results.status_ambiguous" = c("Ambiguous", "Неоднозначный"),
    "results.status_needs_review" = c("Needs review", "Требует проверки"),
    "results.copy_balanced" = c("Balanced", "Сбалансированное"),
    "results.copy_gain" = c("Gain", "Прирост"),
    "results.copy_loss" = c("Loss", "Потеря"),
    "results.copy_compound" = c("Compound", "Составное"),
    "results.copy_neutral" = c("Neutral", "Нейтральное"),
    "results.assignment_assigned" = c("Assigned", "Размещён"),
    "results.assignment_cross_group" = c("Cross-group", "Межгрупповой"),
    "results.assignment_no_recipient_loss" = c("No host", "Нет реципиента"),
    "common.yes" = c("Yes", "Да"),
    "common.no" = c("No", "Нет"),
    "sidebar.example_note" = c(
        "Example {n}: {donor} donor and {recipient} recipient.",
        "Пример {n}: {donor} — донор, {recipient} — реципиент."
    ),
    "tooltip.donor" = c("donor", "донор"),
    "tooltip.recipient" = c("recipient", "реципиент"),
    "tooltip.margin" = c("margin", "отрыв"),
    "tooltip.assigned" = c("assigned", "размещён"),
    "tooltip.cross_group_candidate" = c("cross-group candidate", "межгрупповой кандидат"),
    "tooltip.group_only" = c("group only", "только группа"),
    "tooltip.evidence_donor" = c("donor", "донор"),
    "tooltip.evidence_recipient" = c("recipient", "реципиент"),
    "tooltip.evidence_donor_recipient" = c("donor+recipient", "донор+реципиент"),
    "tooltip.evidence_none" = c("none", "нет"),
    "tooltip.confidence_high" = c("high", "высокая"),
    "tooltip.confidence_moderate" = c("moderate", "средняя"),
    "tooltip.confidence_low" = c("low", "низкая"),
    "tooltip.copy_balanced" = c("balanced", "сбалансированное"),
    "tooltip.copy_gain" = c("gain", "прирост"),
    "tooltip.copy_loss" = c("loss", "потеря"),
    "tooltip.copy_compound" = c("compound", "составное"),
    "tooltip.copy_neutral" = c("neutral", "нейтральное"),
    "tooltip.total_gain" = c("total gain", "суммарный прирост"),
    "tooltip.total_loss" = c("total loss", "суммарная потеря"),
    "errors.prefix" = c("[%s] %s", "[%s] %s"),
    "errors.suggestion.config" = c(
        "Check MAMAMIA_CONFIG and the YAML species/assembly registry; validate it with readMaMaMIAConfig().",
        "Проверьте MAMAMIA_CONFIG и YAML-реестр видов/сборок; проверьте его с помощью readMaMaMIAConfig()."
    ),
    "errors.suggestion.BAM" = c(
        "Check that the file is a coordinate-sorted BAM; re-upload it or choose another file.",
        "Проверьте, что файл — координатно-отсортированный BAM; загрузите его снова или выберите другой файл."
    ),
    "errors.suggestion.assembly" = c(
        "Choose a species pair or assembly override that matches the BAM header, or use a BAM that contains every configured target.",
        "Выберите пару видов или переопределение сборки, соответствующее заголовку BAM, либо используйте BAM, содержащий все настроенные цели."
    ),
    "errors.suggestion.index" = c(
        "Check that the session temporary directory is writable and has enough free space, then re-upload the BAM.",
        "Проверьте, что временный каталог сеанса доступен для записи и на нём достаточно свободного места, затем загрузите BAM снова."
    ),
    "errors.suggestion.windows" = c(
        "Choose a lower resolution; a contig is shorter than the requested window count.",
        "Выберите меньшее разрешение; контиг короче запрошенного числа окон."
    ),
    "errors.suggestion.GC" = c(
        "Check that every gcrank index in the configuration exists and has a supported format_version.",
        "Проверьте, что каждый gcrank-индекс в конфигурации существует и имеет поддерживаемую format_version."
    ),
    "errors.suggestion.metrics" = c(
        "Check MAPQ and the exclude flags; verify that the BAM is intact and indexed.",
        "Проверьте MAPQ и флаги исключения; убедитесь, что BAM цел и проиндексирован."
    ),
    "errors.suggestion.score" = c(
        "Check that donor and recipient windows share chromosome groups and that every group has at least one donor and one recipient.",
        "Проверьте, что окна донора и реципиента относятся к одним хромосомным группам и в каждой группе есть хотя бы один донор и один реципиент."
    ),
    "errors.suggestion.mask" = c(
        "Check that the score table and window grid use the same windows; cohort recurrence requires at least two samples.",
        "Проверьте, что таблица оценок и сетка окон используют одни и те же окна; рекуррентность когорты требует минимум два образца."
    ),
    "errors.suggestion.segmentation" = c(
        "Check that both the donor and the recipient side have enough finite windows.",
        "Проверьте, что и на стороне донора, и на стороне реципиента достаточно конечных окон."
    ),
    "errors.suggestion.calibration" = c(
        "Check the calibration length bins and the control segments; a very small cohort needs the configured single-sample floors.",
        "Проверьте интервалы длины калибровки и контрольные сегменты; очень малой когорте нужны настроенные пороги для одного образца."
    ),
    "errors.suggestion.limits" = c(
        "Reduce the number of BAMs or their total size; the limits come from server.max_bams and server.max_total_bytes.",
        "Уменьшите число BAM-файлов или их суммарный размер; лимиты задаются в server.max_bams и server.max_total_bytes."
    ),
    "errors.suggestion.upload" = c(
        "Reduce the BAM locally with bamhusk (%s), then upload the smaller file: %s",
        "Уменьшите BAM локально с помощью bamhusk (%s), затем загрузите меньший файл: %s"
    ),
    "errors.suggestion.analysis" = c(
        "Check the input mode, the registry defaults (MAPQ, exclude flags, fit limits) and the server log.",
        "Проверьте режим входных данных, значения реестра по умолчанию (MAPQ, флаги исключения, лимиты подгонки) и журнал сервера."
    ),
    "log.preflight_none" = c(
        "Preflight: no state available.",
        "Проверка: состояние недоступно."
    ),
    "log.preflight_waiting" = c(
        "Preflight: waiting for a BAM upload.",
        "Проверка: ожидание загрузки BAM."
    ),
    "log.preflight_failed" = c(
        "Preflight failed.",
        "Проверка не пройдена."
    ),
    "log.preflight_failed_stage" = c(
        "Preflight failed at stage '%s': %s",
        "Проверка не пройдена на этапе '%s': %s"
    ),
    "log.assembly_candidates" = c(
        "Assembly candidates: %s",
        "Кандидаты сборок: %s"
    ),
    "log.assembly_candidate_item" = c(
        "%s (%d/%d targets, %.1f%%)",
        "%s (%d/%d целей, %.1f%%)"
    ),
    "log.preflight_ok" = c(
        "Preflight OK: %d BAM file%s.",
        "Проверка пройдена: BAM-файлов: %d."
    ),
    "log.preflight_note" = c(
        "Preflight note: %s",
        "Примечание проверки: %s"
    ),
    "log.detected_assembly" = c(
        "Detected assembly: %s.",
        "Определена сборка: %s."
    ),
    "log.run_params_rejected" = c(
        "Run parameters rejected at stage '%s': %s",
        "Параметры запуска отклонены на этапе '%s': %s"
    ),
    "log.run_refused_stage" = c(
        "Run refused at stage '%s': %s",
        "Запуск отклонён на этапе '%s': %s"
    ),
    "log.run_failed_stage" = c(
        "Run failed at stage '%s': %s",
        "Запуск завершился ошибкой на этапе '%s': %s"
    ),
    "errors.file_label" = c("File", "Файл"),
    "errors.file_label_plural" = c("Files", "Файлы"),
    "errors.file_more" = c("(+%d more)", "(ещё %d)")
)
for (key in names(extras)) {
    rows[[key]] <- list(en = extras[[key]][[1L]], ru = extras[[key]][[2L]])
}

keys <- sort(names(rows), method = "radix")
en <- vapply(keys, function(key) rows[[key]]$en, character(1))
ru <- vapply(keys, function(key) rows[[key]]$ru, character(1))
empty <- keys[!nzchar(en) | !nzchar(ru)]
if (length(empty) > 0L) {
    stop("Empty dictionary entries: ", paste(empty, collapse = ", "))
}

quote_r <- function(x) {
    encodeString(enc2utf8(as.character(x)), quote = '"')
}
render_dict <- function(values) {
    paste0(
        "    ", quote_r(names(values)), " = ", quote_r(unname(values)),
        collapse = ",\n"
    )
}

out <- c(
    "## Russian/English UI dictionary and translation helpers.",
    "##",
    "## Generated from docs/specs/russian-copy.md by tools/generate-i18n.R.",
    "## The committed file is plain data; no parser runs at application start.",
    "## Values are keyed by dot-path and are identical in both languages.",
    "",
    "mamamia_i18n_dicts <- list(",
    "  en = c(",
    render_dict(en),
    "  ),",
    "  ru = c(",
    render_dict(ru),
    "  )",
    ")",
    "",
    "#' Languages supported by the UI",
    "#' @noRd",
    "mamamia_i18n_langs <- function() c(\"en\", \"ru\")",
    "",
    "#' Normalize a language tag; only \"ru\" is accepted, everything else is en",
    "#' @noRd",
    "mamamia_i18n_valid_lang <- function(lang) {",
    "    lang <- tolower(trimws(as.character(lang)[[1L]]))",
    "    if (identical(lang, \"ru\")) \"ru\" else \"en\"",
    "}",
    "",
    "#' Default interface language",
    "#' @noRd",
    "mamamia_language_default <- function() \"en\"",
    "",
    "#' Flattened dictionary keys for one language",
    "#' @noRd",
    "mamamia_i18n_keys <- function(lang = \"en\") {",
    "    names(mamamia_i18n_dicts[[mamamia_i18n_valid_lang(lang)]])",
    "}",
    "",
    "#' Translate a dot-path key",
    "#'",
    "#' The template itself is trusted copy and is returned verbatim. Named",
    "#' `{value}` placeholders are substituted from `...`; they are inserted as",
    "#' plain text by default so a text sink such as [shiny::p()] escapes them",
    "#' once. Pass `.escape = TRUE` when the result is inserted into a raw HTML",
    "#' sink (for example [shiny::HTML()] or the karyotype payload consumed by",
    "#' JavaScript) so each interpolated value is escaped exactly once.",
    "#'",
    "#' Missing Russian keys fall back to English. A key missing in both",
    "#' languages warns and returns the key itself.",
    "#' @noRd",
    "mamamia_tr <- function(key, lang = mamamia_language_default(), ...,",
    "                       .escape = FALSE) {",
    "    key <- as.character(key)[[1L]]",
    "    lang <- mamamia_i18n_valid_lang(lang)",
    "    value <- as.list(mamamia_i18n_dicts[[lang]])[[key]]",
    "    if (is.null(value) || length(value) == 0L || is.na(value) ||",
    "        !nzchar(value)) {",
    "        value <- as.list(mamamia_i18n_dicts$en)[[key]]",
    "    }",
    "    if (length(value) == 0L || is.na(value) || !nzchar(value)) {",
    "        warning(\"i18n: missing key '\", key, \"'\", call. = FALSE)",
    "        return(key)",
    "    }",
    "    substitutions <- list(...)",
    "    if (length(substitutions) > 0L) {",
    "        if (is.null(names(substitutions)) ||",
    "            any(!nzchar(names(substitutions)))) {",
    "            stop(\"i18n placeholders must be named\", call. = FALSE)",
    "        }",
    "        for (name in names(substitutions)) {",
    "            replacement <- paste(",
    "                as.character(substitutions[[name]]), collapse = \"\"",
    "            )",
    "            if (isTRUE(.escape)) {",
    "                replacement <- htmltools::htmlEscape(replacement)",
    "            }",
    "            value <- gsub(",
    "                paste0(\"{\", name, \"}\"), replacement, value,",
    "                fixed = TRUE",
    "            )",
    "        }",
    "    }",
    "    unname(value)",
    "}"
)
writeLines(out, target_path, useBytes = TRUE)
message("Wrote ", target_path, " with ", length(keys), " keys")
