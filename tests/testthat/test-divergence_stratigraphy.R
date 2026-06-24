
query_cds    <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
subject_cds  <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

test_that("divergence_stratigraphy() errors on invalid ortho_detection ...", {
        expect_error(
                divergence_stratigraphy(
                        query_file      = query_cds,
                        subject_file    = subject_cds,
                        ortho_detection = "INVALID"
                )
        )
})

test_that("divergence_stratigraphy() errors on invalid aligner ...", {
        expect_error(
                divergence_stratigraphy(
                        query_file   = query_cds,
                        subject_file = subject_cds,
                        aligner      = "INVALID"
                )
        )
})

test_that("divergence_stratigraphy() returns a data frame ...\n
          and the output is non-empty", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                comp_cores      = 1,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true(is.data.frame(result))
        expect_gt(nrow(result), 0L)
})

test_that("divergence_stratigraphy() returns DS column with ds.values = TRUE ...\n
          and DS column is numeric ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                ds.values       = TRUE,
                comp_cores      = 1,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true("DS" %in% colnames(result))
        expect_true(is.numeric(result[["DS"]]))
})

test_that("divergence_stratigraphy() returns dNdS column with ds.values = FALSE ...\n
          and dNdS column is numeric with ds.values = FALSE ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                ds.values       = FALSE,
                comp_cores      = 1,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true("dNdS" %in% colnames(result))
        expect_true(is.numeric(result[["dNdS"]]))
})

test_that("divergence_stratigraphy() always returns query_id ...\n
          and excludes subject_id by default ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                comp_cores      = 1,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true("query_id" %in% colnames(result))
        expect_false("subject_id" %in% colnames(result))
})

test_that("divergence_stratigraphy() includes subject_id with subject.id = TRUE ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                subject.id      = TRUE,
                comp_cores      = 1,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true("subject_id" %in% colnames(result))
})

test_that("divergence_stratigraphy() works with BH ortho_detection ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "BH",
                comp_cores      = 1,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true(is.data.frame(result))
        expect_gt(nrow(result), 0L)
})

test_that("divergence_stratigraphy() works with NG method ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                dnds_est.method = "NG",
                comp_cores      = 1,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true(is.data.frame(result))
        expect_true("DS" %in% colnames(result))
})

test_that("divergence_stratigraphy() works with comp_cores = 2 ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        # PSOCK workers need an installed package; skip when loaded via load_all()
        skip_if(
                isTRUE(tryCatch(pkgload::is_dev_package("orthologr"), error = function(e) FALSE)),
                "orthologr loaded via load_all(); PSOCK workers require an installed package"
        )
        result <- divergence_stratigraphy(
                query_file      = query_cds,
                subject_file    = subject_cds,
                ortho_detection = "RBH",
                comp_cores      = 2,
                quiet           = TRUE,
                clean_folders   = TRUE
        )
        expect_true(is.data.frame(result))
        expect_gt(nrow(result), 0L)
})

