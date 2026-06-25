
protein_file  <- system.file("seqs/ortho_thal_aa.fasta",  package = "orthologr")
protein_file2 <- system.file("seqs/ortho_lyra_aa.fasta",  package = "orthologr")
cds_file      <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")

expected_cols <- c("representative_id", "member_id")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

# --- Output structure --------------------------------------------------------

test_that("deepclust() returns a tibble", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_true(tibble::is_tibble(result))
})

test_that("deepclust() output has exactly 2 expected columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_equal(colnames(result), expected_cols)
})

test_that("deepclust() output columns are character type", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_type(result[["representative_id"]], "character")
        expect_type(result[["member_id"]],         "character")
})

test_that("deepclust() output is non-empty", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)
        expect_gt(nrow(result), 0L)
})

# --- Cluster membership invariants -------------------------------------------

test_that("every input sequence appears as a cluster member", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)

        input_ids <- names(seqinr::read.fasta(protein_file, seqtype = "AA",
                                               as.string = TRUE))
        expect_true(all(input_ids %in% result[["member_id"]]))
})

test_that("all representative_ids are themselves cluster members", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = protein_file)

        expect_true(all(result[["representative_id"]] %in% result[["member_id"]]))
})

# --- seq_type = "cds" --------------------------------------------------------

test_that("deepclust() runs with seq_type = 'cds'", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = cds_file, seq_type = "cds")
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
        expect_gt(nrow(result), 0L)
})

# --- Multi-file input --------------------------------------------------------

test_that("deepclust() accepts a vector of input files and returns a tibble", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = c(protein_file, protein_file2))
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
        expect_gt(nrow(result), 0L)
})

test_that("deepclust() with multiple files contains all sequences from all files", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- deepclust(input_file = c(protein_file, protein_file2))

        ids1 <- names(seqinr::read.fasta(protein_file,  seqtype = "AA", as.string = TRUE))
        ids2 <- names(seqinr::read.fasta(protein_file2, seqtype = "AA", as.string = TRUE))

        expect_true(all(c(ids1, ids2) %in% result[["member_id"]]))
})

test_that("deepclust() errors informatively when a file does not exist", {
        expect_error(
                deepclust(input_file = c(protein_file, "nonexistent_file.fasta")),
                regexp = "not found"
        )
})

# --- approx_id threshold -----------------------------------------------------

test_that("stricter approx_id yields at least as many rows as a lenient threshold", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")

        # Higher approx_id -> fewer sequences can cluster together -> more singleton rows
        result_lenient <- deepclust(input_file = protein_file, approx_id = 30)
        result_strict  <- deepclust(input_file = protein_file, approx_id = 90)
        expect_gte(nrow(result_strict), nrow(result_lenient))
})

# --- save.output -------------------------------------------------------------

test_that("deepclust() writes a TSV when save.output is specified", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        out_dir  <- tempdir()
        filename <- paste0("deepclust_", basename(protein_file), ".tsv")

        deepclust(input_file = protein_file, save.output = out_dir)

        expect_true(file.exists(file.path(out_dir, filename)))
})
