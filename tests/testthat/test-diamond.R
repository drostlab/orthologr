context("Test: diamond()")

cds_file     <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
subject_cds  <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")
protein_file <- system.file("seqs/ortho_thal_aa.fasta",  package = "orthologr")
subject_prot <- system.file("seqs/ortho_lyra_aa.fasta",  package = "orthologr")

expected_cols <- c(
        "query_id", "subject_id", "perc_identity", "num_ident_matches",
        "alig_length", "mismatches", "gap_openings", "n_gaps", "pos_match",
        "ppos", "q_start", "q_end", "q_len", "qcov", "qcovhsp",
        "s_start", "s_end", "s_len", "evalue", "bit_score", "score_raw"
)

diamond_available <- tryCatch({
        is_installed_diamond()
        TRUE
}, error = function(e) FALSE)

# --- Input validation (no DIAMOND needed) ------------------------------------

test_that("diamond() errors on invalid diamond_algorithm", {
        expect_error(
                diamond(query_file = cds_file, subject_file = subject_cds,
                        diamond_algorithm = "blastn"),
                regexp = "valid DIAMOND mode"
        )
})

test_that("diamond() errors on invalid database_maker", {
        expect_error(
                diamond(query_file = cds_file, subject_file = subject_cds,
                        database_maker = "hmmer"),
                regexp = "Please choose either"
        )
})

# --- Output structure --------------------------------------------------------

test_that("diamond() returns a tibble with CDS input", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(query_file = cds_file, subject_file = subject_cds)
        expect_true(tibble::is_tibble(result))
})

test_that("diamond() output has exactly 21 expected columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(query_file = cds_file, subject_file = subject_cds)
        expect_equal(colnames(result), expected_cols)
})

test_that("diamond() output is non-empty for related sequences", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(query_file = cds_file, subject_file = subject_cds)
        expect_gt(nrow(result), 0L)
})

test_that("diamond() query_id and subject_id columns are character", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(query_file = cds_file, subject_file = subject_cds)
        expect_type(result[["query_id"]],   "character")
        expect_type(result[["subject_id"]], "character")
})

test_that("diamond() evalue column respects the e-value cutoff", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(query_file = cds_file, subject_file = subject_cds,
                          eval = "1E-5")
        expect_true(all(result[["evalue"]] <= 1e-5))
})

test_that("diamond() perc_identity is between 0 and 100", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(query_file = cds_file, subject_file = subject_cds)
        expect_true(all(result[["perc_identity"]] >= 0))
        expect_true(all(result[["perc_identity"]] <= 100))
})

# --- Protein sequence input --------------------------------------------------

test_that("diamond() runs with seq_type = 'protein'", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(
                query_file   = protein_file,
                subject_file = subject_prot,
                seq_type     = "protein"
        )
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
        expect_gt(nrow(result), 0L)
})

# --- Sensitivity modes -------------------------------------------------------

test_that("diamond() runs with sensitivity_mode = 'sensitive'", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- diamond(
                query_file       = cds_file,
                subject_file     = subject_cds,
                sensitivity_mode = "sensitive"
        )
        expect_true(tibble::is_tibble(result))
        expect_gt(nrow(result), 0L)
})

test_that("diamond() with 'more-sensitive' mode returns at least as many hits as 'fast' mode", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result_fast      <- diamond(query_file = cds_file, subject_file = subject_cds,
                                    sensitivity_mode = "fast")
        result_sensitive <- diamond(query_file = cds_file, subject_file = subject_cds,
                                    sensitivity_mode = "more-sensitive")
        expect_gte(nrow(result_sensitive), nrow(result_fast))
})
