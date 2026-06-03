context("Test: blast_best()")

blast_available <- isTRUE(tryCatch(is_installed_blast(), error = function(e) FALSE))

test_that("blast_best() runs properly ...", {
        skip_if_not(blast_available, "BLAST is not installed or not on PATH")
        result <- blast_best(
                query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr')
        )
        expect_true(tibble::is_tibble(result))
        expect_equal(nrow(result), dplyr::n_distinct(result[["query_id"]]))
})
