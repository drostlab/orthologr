
cds_file     <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
subject_cds  <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")
protein_file <- system.file("seqs/ortho_thal_aa.fasta",  package = "orthologr")
subject_prot <- system.file("seqs/ortho_lyra_aa.fasta",  package = "orthologr")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

test_that("diamond_best() returns a tibble with CDS input", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_best(query_file = cds_file, subject_file = subject_cds)
        expect_true(tibble::is_tibble(result))
})

test_that("diamond_best() output is non-empty for related sequences", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_best(query_file = cds_file, subject_file = subject_cds)
        expect_gt(nrow(result), 0L)
})

test_that("diamond_best() returns at most one hit per query_id", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_best(query_file = cds_file, subject_file = subject_cds)
        expect_equal(nrow(result), dplyr::n_distinct(result[["query_id"]]))
})

test_that("diamond_best() returns fewer or equal hits compared to diamond()", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        all_hits  <- diamond(query_file = cds_file, subject_file = subject_cds)
        best_hits <- diamond_best(query_file = cds_file, subject_file = subject_cds)
        expect_lte(nrow(best_hits), nrow(all_hits))
})

test_that("diamond_best() has query_id column of type character", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_best(query_file = cds_file, subject_file = subject_cds)
        expect_type(result[["query_id"]], "character")
})

test_that("diamond_best() runs with seq_type = 'protein'", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_best(
                query_file   = protein_file,
                subject_file = subject_prot,
                seq_type     = "protein"
        )
        expect_true(tibble::is_tibble(result))
        expect_gt(nrow(result), 0L)
        expect_equal(nrow(result), dplyr::n_distinct(result[["query_id"]]))
})
