
cds_file     <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
subject_cds  <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")
protein_file <- system.file("seqs/ortho_thal_aa.fasta",  package = "orthologr")
subject_prot <- system.file("seqs/ortho_lyra_aa.fasta",  package = "orthologr")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

test_that("diamond_rec() returns a tibble with CDS input", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_rec(query_file = cds_file, subject_file = subject_cds)
        expect_true(tibble::is_tibble(result))
})

test_that("diamond_rec() output is non-empty for related sequences", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_rec(query_file = cds_file, subject_file = subject_cds)
        expect_gt(nrow(result), 0L)
})

test_that("diamond_rec() returns fewer or equal hits than diamond_best()", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        bh  <- diamond_best(query_file = cds_file, subject_file = subject_cds)
        rbh <- diamond_rec(query_file  = cds_file, subject_file = subject_cds)
        expect_lte(nrow(rbh), nrow(bh))
})

test_that("diamond_rec() query_id and subject_id are unique pairs", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_rec(query_file = cds_file, subject_file = subject_cds)
        n_pairs <- dplyr::n_distinct(result[["query_id"]], result[["subject_id"]])
        expect_equal(nrow(result), n_pairs)
})

test_that("diamond_rec() all query_ids are unique (each gene has one RBH)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_rec(query_file = cds_file, subject_file = subject_cds)
        expect_equal(nrow(result), dplyr::n_distinct(result[["query_id"]]))
})

test_that("diamond_rec() runs with seq_type = 'protein'", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond_rec(
                query_file   = protein_file,
                subject_file = subject_prot,
                seq_type     = "protein"
        )
        expect_true(tibble::is_tibble(result))
        expect_gt(nrow(result), 0L)
})

test_that("diamond_rec() RBH pairs are a subset of forward BH pairs", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        bh  <- diamond_best(query_file = cds_file, subject_file = subject_cds)
        rbh <- diamond_rec(query_file  = cds_file, subject_file = subject_cds)
        # Every RBH query_id must appear in BH
        expect_true(all(rbh[["query_id"]] %in% bh[["query_id"]]))
})
