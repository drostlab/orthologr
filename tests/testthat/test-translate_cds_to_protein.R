context("Test: translate_cds_to_protein()")

cds_file <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")

test_that("translate_cds_to_protein() errors on a non-existent input file", {
        expect_error(
                translate_cds_to_protein("nonexistent_cds.fasta", tempfile(fileext = ".fasta")),
                regexp = "seems not to exist"
        )
})

test_that("translate_cds_to_protein() creates an output file", {
        out <- tempfile(fileext = ".fasta")
        on.exit(unlink(out))
        translate_cds_to_protein(input_file = cds_file, output_file = out)
        expect_true(file.exists(out))
})

test_that("translate_cds_to_protein() output is a non-empty FASTA file", {
        out <- tempfile(fileext = ".fasta")
        on.exit(unlink(out))
        translate_cds_to_protein(input_file = cds_file, output_file = out)
        aa_seqs <- Biostrings::readAAStringSet(out)
        expect_gt(length(aa_seqs), 0L)
})

test_that("translate_cds_to_protein() output contains only amino acid characters", {
        out <- tempfile(fileext = ".fasta")
        on.exit(unlink(out))
        translate_cds_to_protein(input_file = cds_file, output_file = out)
        aa_seqs <- Biostrings::readAAStringSet(out)
        # All sequences should be of class AAStringSet (Biostrings validates this)
        expect_s4_class(aa_seqs, "AAStringSet")
})

test_that("translate_cds_to_protein() output has same number of sequences as input", {
        out <- tempfile(fileext = ".fasta")
        on.exit(unlink(out))
        translate_cds_to_protein(input_file = cds_file, output_file = out,
                                 delete_corrupt_cds = FALSE)
        cds_in <- Biostrings::readDNAStringSet(cds_file)
        aa_out <- Biostrings::readAAStringSet(out)
        expect_equal(length(aa_out), length(cds_in))
})

test_that("translate_cds_to_protein() output has 1/3 width as input (since 3 nucleotides corresponds to one amino acid)", {
        out <- tempfile(fileext = ".fasta")
        on.exit(unlink(out))
        translate_cds_to_protein(input_file = cds_file, output_file = out,
                                 delete_corrupt_cds = FALSE)
        cds_in <- Biostrings::readDNAStringSet(cds_file)
        aa_out <- Biostrings::readAAStringSet(out)
        expect_equal(width(aa_out), width(cds_in) / 3)
})