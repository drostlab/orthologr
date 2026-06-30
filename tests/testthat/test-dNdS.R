
query_cds    <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
subject_cds  <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

test_that("dNdS() errors on invalid ortho_detection ...", {
        expect_error(
                dNdS(
                        query_file      = query_cds,
                        subject_file    = subject_cds,
                        ortho_detection = "INVALID"
                )
        )
})

test_that("dNdS() errors on invalid aa_aln_type ...", {
        expect_error(
                dNdS(
                        query_file   = query_cds,
                        subject_file = subject_cds,
                        aa_aln_type  = "INVALID"
                )
        )
})

test_that("dNdS() errors on invalid codon_aln_tool ...", {
        expect_error(
                dNdS(
                        query_file     = query_cds,
                        subject_file   = subject_cds,
                        codon_aln_tool = "INVALID"
                )
        )
})

test_that("dNdS() errors on invalid dnds_est.method ...", {
        expect_error(
                dNdS(
                        query_file      = query_cds,
                        subject_file    = subject_cds,
                        dnds_est.method = "INVALID"
                )
        )
})

test_that("dNdS() returns a data frame with RBH and Comeron ...\n
          and is non-empty with RBH and Comeron ...\n
          and has required columns with RBH and Comeron ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- dNdS(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                aa_aln_type     = "pairwise",
                aa_aln_tool     = "NW",
                codon_aln_tool  = "pal2nal",
                dnds_est.method = "Comeron",
                comp_cores      = 1,
                quiet           = TRUE,
                print_citation  = FALSE
        )
        expect_true(is.data.frame(result))
        expect_gt(nrow(result), 0L)
        expect_true(all(c("query_id", "subject_id", "dN", "dS", "dNdS") %in% colnames(result)))
})

test_that("dNdS() dNdS column is numeric with RBH and Comeron ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- dNdS(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                aa_aln_type     = "pairwise",
                aa_aln_tool     = "NW",
                codon_aln_tool  = "pal2nal",
                dnds_est.method = "Comeron",
                comp_cores      = 1,
                quiet           = TRUE,
                print_citation  = FALSE
        )
        expect_true(is.numeric(result[["dNdS"]]))
})

test_that("dNdS() returns valid results with BH method ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- dNdS(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "BH",
                aa_aln_type     = "pairwise",
                aa_aln_tool     = "NW",
                codon_aln_tool  = "pal2nal",
                dnds_est.method = "Comeron",
                comp_cores      = 1,
                quiet           = TRUE,
                print_citation  = FALSE
        )
        expect_true(is.data.frame(result))
        expect_gt(nrow(result), 0L)
})

test_that("dNdS() works with Li estimation method ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- dNdS(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                aa_aln_type     = "pairwise",
                aa_aln_tool     = "NW",
                codon_aln_tool  = "pal2nal",
                dnds_est.method = "Li",
                comp_cores      = 1,
                quiet           = TRUE,
                print_citation  = FALSE
        )
        expect_true(is.data.frame(result))
        expect_gt(nrow(result), 0L)
})
