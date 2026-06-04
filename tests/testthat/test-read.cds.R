context("Test: read.cds()")

cds_file <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")

test_that("read.cds() returns a data.table", {
        result <- read.cds(file = cds_file, format = "fasta")
        expect_true(data.table::is.data.table(result))
})

test_that("read.cds() output has columns 'geneids' and 'seqs'", {
        result <- read.cds(file = cds_file, format = "fasta")
        expect_true("geneids" %in% colnames(result))
        expect_true("seqs"    %in% colnames(result))
})

test_that("read.cds() returns non-empty data for a valid file", {
        result <- read.cds(file = cds_file, format = "fasta")
        expect_gt(nrow(result), 0L)
})

test_that("read.cds() 'seqs' column contains lowercase nucleotide strings", {
        result <- read.cds(file = cds_file, format = "fasta")
        # All characters should be lowercase a/c/g/t
        expect_true(all(grepl("^[acgt]+$", result[["seqs"]])))
})

test_that("read.cds() geneids are unique", {
        result <- read.cds(file = cds_file, format = "fasta")
        expect_equal(length(unique(result[["geneids"]])), nrow(result))
})

test_that("read.cds() errors on a non-existent file", {
        expect_error(
                read.cds(file = "nonexistent.fasta", format = "fasta"),
                regexp = "seems not to exist|does not seem to exist|could not be read"
        )
})

test_that("read.cds() errors on an unsupported format", {
        expect_error(
                read.cds(file = cds_file, format = "genbank"),
                regexp = "supported"
        )
})

test_that("read.cds() with delete_corrupt_cds = TRUE returns only triplet-length sequences", {
        result <- read.cds(file = cds_file, format = "fasta",
                           delete_corrupt_cds = TRUE)
        # All retained sequences should have length divisible by 3
        expect_true(all(nchar(result[["seqs"]]) %% 3 == 0))
})
