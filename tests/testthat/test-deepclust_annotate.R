
protein_file  <- system.file("seqs/ortho_thal_aa.fasta", package = "orthologr")
protein_file2 <- system.file("seqs/ortho_lyra_aa.fasta", package = "orthologr")

expected_cols <- c("representative_id", "member_id", "file_name")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

# --- Output structure --------------------------------------------------------

test_that("deepclust_annotate() returns a tibble with 3 columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(input_file = c(protein_file, protein_file2))
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
})

test_that("deepclust_annotate() output is non-empty", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(input_file = c(protein_file, protein_file2))
        expect_gt(nrow(result), 0L)
})

test_that("deepclust_annotate() has no missing file_name values", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(input_file = c(protein_file, protein_file2))
        expect_false(anyNA(result[["file_name"]]))
})

# --- Character vector input: file_name uses basenames ------------------------

test_that("deepclust_annotate() file_name column contains only source basenames (character vector)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(input_file = c(protein_file, protein_file2))
        expected_names <- c(basename(protein_file), basename(protein_file2))
        expect_true(all(result[["file_name"]] %in% expected_names))
})

test_that("deepclust_annotate() file_name contains both source basenames (character vector)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(input_file = c(protein_file, protein_file2))
        expected_names <- c(basename(protein_file), basename(protein_file2))
        expect_true(all(expected_names %in% result[["file_name"]]))
})

# --- Named list input: file_name uses list names -----------------------------

test_that("deepclust_annotate() uses list names as file_name when input_file is a named list", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(
                input_file = list(
                        "thal" = protein_file,
                        "lyra" = protein_file2
                )
        )
        expect_true(all(result[["file_name"]] %in% c("thal", "lyra")))
})

test_that("deepclust_annotate() file_name contains all list names (named list)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(
                input_file = list(
                        "thal" = protein_file,
                        "lyra" = protein_file2
                )
        )
        expect_true(all(c("thal", "lyra") %in% result[["file_name"]]))
})

test_that("deepclust_annotate() named list does not produce basenames in file_name", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(
                input_file = list(
                        "thal" = protein_file,
                        "lyra" = protein_file2
                )
        )
        expect_false(any(result[["file_name"]] %in% c(basename(protein_file),
                                                       basename(protein_file2))))
})

test_that("deepclust_annotate() named list returns a tibble with 3 expected columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(
                input_file = list(
                        "thal" = protein_file,
                        "lyra" = protein_file2
                )
        )
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
})

# --- Unnamed list input: file_name falls back to basenames -------------------

test_that("deepclust_annotate() falls back to basenames when input_file is an unnamed list", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust_annotate(
                input_file = list(protein_file, protein_file2)
        )
        expected_names <- c(basename(protein_file), basename(protein_file2))
        expect_true(all(result[["file_name"]] %in% expected_names))
})

# --- Duplicate ID warning ----------------------------------------------------

test_that("deepclust_annotate() warns when duplicate sequence IDs exist across input files", {
        # Use cluster_table to bypass diamond deepclust — no DIAMOND2 needed
        first_id <- names(seqinr::read.fasta(protein_file, seqtype = "AA",
                                              as.string = TRUE))[1]

        tmp_dup <- tempfile(fileext = ".fasta")
        file.copy(protein_file, tmp_dup)
        on.exit(unlink(tmp_dup), add = TRUE)

        dummy_clusters <- tibble::tibble(
                representative_id = first_id,
                member_id         = first_id
        )

        expect_warning(
                deepclust_annotate(
                        input_file    = c(protein_file, tmp_dup),
                        cluster_table = dummy_clusters
                ),
                regexp = "Duplicate sequence IDs"
        )
})

test_that("deepclust_annotate() duplicate warning names affected IDs", {
        first_id <- names(seqinr::read.fasta(protein_file, seqtype = "AA",
                                              as.string = TRUE))[1]

        tmp_dup <- tempfile(fileext = ".fasta")
        file.copy(protein_file, tmp_dup)
        on.exit(unlink(tmp_dup), add = TRUE)

        dummy_clusters <- tibble::tibble(
                representative_id = first_id,
                member_id         = first_id
        )

        expect_warning(
                deepclust_annotate(
                        input_file    = c(protein_file, tmp_dup),
                        cluster_table = dummy_clusters
                ),
                regexp = first_id,
                fixed  = TRUE
        )
})

test_that("deepclust_annotate() does not warn about duplicate IDs when all IDs are unique", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")

        # thal and lyra have completely distinct IDs — no warning expected
        warned <- FALSE
        withCallingHandlers(
                deepclust_annotate(input_file = c(protein_file, protein_file2)),
                warning = function(w) {
                        if (grepl("Duplicate sequence IDs", conditionMessage(w)))
                                warned <<- TRUE
                        invokeRestart("muffleWarning")
                }
        )
        expect_false(warned)
})

# --- Error handling ----------------------------------------------------------

test_that("deepclust_annotate() errors informatively when a path in a named list does not exist", {
        expect_error(
                deepclust_annotate(
                        input_file = list(
                                "thal"    = protein_file,
                                "missing" = "nonexistent_file.fasta"
                        )
                ),
                regexp = "not found"
        )
})

test_that("deepclust_annotate() errors informatively when a path in a character vector does not exist", {
        expect_error(
                deepclust_annotate(
                        input_file = c(protein_file, "nonexistent_file.fasta")
                ),
                regexp = "not found"
        )
})
