context("Test: blast()")

blast_available <- isTRUE(tryCatch(is_installed_blast(), error = function(e) FALSE))

test_that("blast() runs properly ...", {
        skip_if_not(blast_available, "BLAST is not installed or not on PATH")
        test_blast <- blast(
                query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr')
        )
        expect_true(tibble::is_tibble(test_blast))
})
