## BEDPE export, writing the segments flagged as putative introgressions with
## both chromosome coordinates.

#' @export
#' @importFrom rlang .data
writeToBedpe <- function(ISA, file) {
    stopifnot("Input is not ISA object" = inherits(ISA, "ISA"))

    output_colnames <- c(
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "name",
        "score",
        "strand1",
        "strand2",
        "chrom_name1",
        "chrom_name2",
        "subgen_name1",
        "subgen_name2"
    )

    putative <- dplyr::filter(ISA$segments, .data$putative_introgression)

    if (nrow(putative) == 0L) {
        warning("No putative introgressions to write; writing header only",
            call. = FALSE
        )
        empty <- data.frame(matrix(
            ncol = length(output_colnames),
            nrow = 0L,
            dimnames = list(NULL, output_colnames)
        ))
        utils::write.table(
            empty,
            file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        message("Output was written in ", file)
        return(invisible(NULL))
    }

    annotate_segment_coords(putative, ISA$out) |>
        dplyr::select(dplyr::all_of(
            c(
                "chr_id.don",
                "start.don",
                "stop.don",
                "chr_id.rec",
                "start.rec",
                "stop.rec",
                ISA$meta$metrics,
                "chr_name.don",
                "chr_name.rec",
                "subgenome.don",
                "subgenome.rec"
            )
        )) |>
        dplyr::mutate(
            name = paste("PIS", seq_len(dplyr::n()), sep = ""),
            .after = "stop.rec"
        ) |>
        dplyr::mutate(
            strand1 = ".",
            strand2 = ".",
            .after = ISA$meta$metrics
        ) |>
        stats::setNames(output_colnames) |>
        utils::write.table(file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

    message("Output was written in ", file)
    invisible(NULL)
}
