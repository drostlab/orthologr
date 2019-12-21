context("Test: blast()")

test_that(
        "blast() runs properly ...",
        {
                test_blast <- blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                                    subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
                
                expect_true(tibble::is_tibble(test_blast))
        }
)
