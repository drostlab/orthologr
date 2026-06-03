context("Test: set_diamond()")

cds_file     <- system.file("seqs/ortho_thal_cds.fasta",  package = "orthologr")
protein_file <- system.file("seqs/ortho_thal_aa.fasta",   package = "orthologr")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

# --- Input validation (no DIAMOND needed) ------------------------------------

test_that("set_diamond() errors on non-existent file", {
        expect_error(
                set_diamond(file = "nonexistent_file.fasta"),
                regexp = "seems not to exist"
        )
})

test_that("set_diamond() errors on invalid seq_type", {
        expect_error(
                set_diamond(file = cds_file, seq_type = "rna"),
                regexp = "Please choose either"
        )
})

test_that("set_diamond() errors on invalid makedb_type", {
        expect_error(
                set_diamond(file = cds_file, makedb_type = "nucleotide"),
                regexp = "valid DIAMOND database type"
        )
})

# --- Output structure with seq_type = "cds" (no makedb) ----------------------

test_that("set_diamond() returns a list of length 2 for CDS input", {
        result <- set_diamond(file = cds_file, seq_type = "cds")
        expect_type(result, "list")
        expect_length(result, 2L)
})

test_that("set_diamond()[[1]] is a data.table with expected columns (CDS)", {
        result <- set_diamond(file = cds_file, seq_type = "cds")
        dt <- result[[1]]
        expect_true(data.table::is.data.table(dt))
        expect_true(all(c("geneids", "seqs", "aa") %in% colnames(dt)))
})

test_that("set_diamond() returns non-empty sequence data for CDS input", {
        result <- set_diamond(file = cds_file, seq_type = "cds")
        expect_gt(nrow(result[[1]]), 0L)
})

test_that("set_diamond() translated amino acids are character strings", {
        result <- set_diamond(file = cds_file, seq_type = "cds")
        expect_true(is.character(result[[1]][["aa"]]))
})

# --- Output structure with seq_type = "protein" (no makedb) ------------------

test_that("set_diamond() returns a list of length 2 for protein input", {
        result <- set_diamond(file = protein_file, seq_type = "protein")
        expect_type(result, "list")
        expect_length(result, 2L)
})

test_that("set_diamond()[[1]] has columns geneids and aa for protein input", {
        result <- set_diamond(file = protein_file, seq_type = "protein")
        dt <- result[[1]]
        expect_true(data.table::is.data.table(dt))
        expect_true(all(c("geneids", "aa") %in% colnames(dt)))
        expect_gt(nrow(dt), 0L)
})

# --- makedb = TRUE (requires DIAMOND) ----------------------------------------

test_that("set_diamond() with makedb = TRUE creates a database file", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        result <- set_diamond(file = protein_file, seq_type = "protein", makedb = TRUE)
        dbname <- result[[2]]
        expect_true(nchar(dbname) > 0)
        db_path <- file.path(tempdir(), "_blast_db", paste0(dbname, ".dmnd"))
        expect_true(file.exists(db_path))
})
