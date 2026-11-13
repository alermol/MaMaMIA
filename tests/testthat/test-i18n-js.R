## Language logic in the rendered UI must not depend on English tab text.

test_that("nav panels expose explicit data-value keys", {
    bam <- as.character(bam_results_ui())
    expect_match(bam, 'data-value="karyotype"', fixed = TRUE)
    expect_match(bam, 'data-value="calls"', fixed = TRUE)
    expect_match(bam, 'data-value="export"', fixed = TRUE)

    explorer <- as.character(explorer_results_ui("explorer"))
    expect_match(explorer, 'data-value="karyotype"', fixed = TRUE)
    expect_match(explorer, 'data-value="calls"', fixed = TRUE)
})

test_that("mamamiaResultTabs keys on data-value, not English text", {
    source <- readLines(
        mamamia_path("R", "app-bam.R"), warn = FALSE
    )
    start <- grep("function mamamiaResultTabs", source)
    expect_length(start, 1L)
    finish <- start + which(
        grepl("mamamia-result-tabs", source[start:length(source)])
    )[[1L]] - 1L
    block <- paste(source[start:finish], collapse = "\n")

    expect_match(block, "mamamiaTabValue", fixed = TRUE)
    expect_match(paste(source, collapse = "\n"), "[data-value]", fixed = TRUE)
    expect_false(grepl("!== 'Calls'", block, fixed = TRUE))
    expect_false(grepl("!== 'Export'", block, fixed = TRUE))
    expect_false(grepl("!== 'Karyotype'", block, fixed = TRUE))
    expect_false(grepl("textContent", block, fixed = TRUE))
})

test_that("the language handler updates tabs and the document language", {
    source <- paste(
        readLines(mamamia_path("R", "app-bam.R"), warn = FALSE),
        collapse = "\n"
    )
    expect_match(source, "mamamia-set-language", fixed = TRUE)
    expect_match(source, "document.documentElement.lang", fixed = TRUE)
    expect_match(source, "x.tabs.karyotype", fixed = TRUE)
})
