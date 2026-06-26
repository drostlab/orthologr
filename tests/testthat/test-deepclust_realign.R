
protein_file  <- system.file("seqs/ortho_thal_aa.fasta", package = "orthologr")
protein_file2 <- system.file("seqs/ortho_lyra_aa.fasta", package = "orthologr")

expected_cols <- c("representative_id", "member_id",
                   "approx_pident", "evalue", "bitscore",
                   "qcovhsp", "scovhsp")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

# Helper: run deepclust once and reuse the result across tests that need it.
# Wrapped in skip_if_not so it is only evaluated when DIAMOND2 is present.
get_clusters <- function() {
        deepclust(input_file = c(protein_file, protein_file2))
}

# =============================================================================
# deepclust_realign() -- tibble input
# =============================================================================

# --- Output structure ---------------------------------------------------------

test_that("deepclust_realign() returns a tibble (tibble clusters input)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_true(tibble::is_tibble(result))
})

test_that("deepclust_realign() output has exactly the 7 expected columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_equal(colnames(result), expected_cols)
})

test_that("deepclust_realign() output is non-empty", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_gt(nrow(result), 0L)
})

# --- Column types -------------------------------------------------------------

test_that("deepclust_realign() id columns are character type", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_type(result[["member_id"]],         "character")
        expect_type(result[["representative_id"]], "character")
})

test_that("deepclust_realign() numeric columns are double type", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_type(result[["approx_pident"]], "double")
        expect_type(result[["evalue"]],        "double")
        expect_type(result[["bitscore"]],      "double")
        expect_type(result[["qcovhsp"]],       "double")
        expect_type(result[["scovhsp"]],       "double")
})

# --- Value ranges -------------------------------------------------------------

test_that("approx_pident values are in [0, 100]", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_true(all(result[["approx_pident"]] >= 0 & result[["approx_pident"]] <= 100,
                        na.rm = TRUE))
})

test_that("evalue values are non-negative", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_true(all(result[["evalue"]] >= 0, na.rm = TRUE))
})

test_that("bitscore values are positive", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_true(all(result[["bitscore"]] > 0, na.rm = TRUE))
})

test_that("qcovhsp and scovhsp values are in [0, 100]", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_true(all(result[["qcovhsp"]] >= 0 & result[["qcovhsp"]] <= 100,
                        na.rm = TRUE))
        expect_true(all(result[["scovhsp"]] >= 0 & result[["scovhsp"]] <= 100,
                        na.rm = TRUE))
})

# --- Cluster membership consistency ------------------------------------------

test_that("all member_ids in deepclust_realign output are present in deepclust clusters", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_true(all(result[["member_id"]] %in% clusters[["member_id"]]))
})

test_that("all representative_ids in deepclust_realign output are present in deepclust clusters", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        expect_true(all(result[["representative_id"]] %in% clusters[["representative_id"]]))
})

# =============================================================================
# deepclust_realign() -- file path input for clusters
# =============================================================================

test_that("deepclust_realign() accepts a TSV file path as clusters", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        out_dir <- tempdir()

        # write the deepclust TSV to a known path
        deepclust(
                input_file  = c(protein_file, protein_file2),
                save.output = out_dir
        )
        tsv_path <- file.path(
                out_dir,
                paste0("deepclust_",
                       paste0("deepclust_merged_2_files.fasta"),
                       ".tsv")
        )
        skip_if_not(file.exists(tsv_path), "deepclust TSV not found for file-path test")

        result <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = tsv_path
        )
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
        expect_gt(nrow(result), 0L)
})

# --- save.output -------------------------------------------------------------

test_that("deepclust_realign() writes a TSV when save.output is specified", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        out_dir  <- tempdir()

        deepclust_realign(
                input_file  = c(protein_file, protein_file2),
                clusters    = clusters,
                save.output = out_dir
        )

        written_files <- list.files(out_dir, pattern = "^deepclust_realign_.*\\.tsv$")
        expect_gt(length(written_files), 0L)
})

# --- Single-file input -------------------------------------------------------

test_that("deepclust_realign() works with a single input file", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters_single <- deepclust(input_file = protein_file)
        result <- deepclust_realign(
                input_file = protein_file,
                clusters   = clusters_single
        )
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), expected_cols)
        expect_gt(nrow(result), 0L)
})

# =============================================================================
# deepclust_realign() -- error handling (no DIAMOND2 needed)
# =============================================================================

test_that("deepclust_realign() errors when input_file does not exist", {
        expect_error(
                deepclust_realign(
                        input_file = "nonexistent.fasta",
                        clusters   = data.frame(representative_id = "A",
                                                member_id         = "A")
                ),
                regexp = "not found"
        )
})

test_that("deepclust_realign() errors when clusters data frame is missing required columns", {
        expect_error(
                deepclust_realign(
                        input_file = protein_file,
                        clusters   = data.frame(rep_id = "A", mem_id = "A")
                ),
                regexp = "representative_id"
        )
})

test_that("deepclust_realign() errors when clusters is neither a file path nor a data frame", {
        expect_error(
                deepclust_realign(
                        input_file = protein_file,
                        clusters   = 42L
                ),
                regexp = "clusters"
        )
})

# =============================================================================
# deepclust_annotate() -- cluster_table parameter
# =============================================================================

test_that("deepclust_annotate() with cluster_table returns a tibble with 3 expected columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters
        )
        expect_true(tibble::is_tibble(result))
        expect_equal(colnames(result), c("representative_id", "member_id", "file_name"))
})

test_that("deepclust_annotate() cluster_table result is non-empty", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters
        )
        expect_gt(nrow(result), 0L)
})

test_that("deepclust_annotate() cluster_table has no missing file_name values", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        result <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters
        )
        expect_false(anyNA(result[["file_name"]]))
})

test_that("deepclust_annotate() cluster_table produces same rows as a fresh deepclust_annotate run", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters  <- get_clusters()
        from_table <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters
        )
        fresh <- deepclust_annotate(input_file = c(protein_file, protein_file2))
        expect_equal(nrow(from_table), nrow(fresh))
        expect_setequal(from_table[["member_id"]], fresh[["member_id"]])
})

test_that("deepclust_annotate() errors when cluster_table is missing required columns", {
        expect_error(
                deepclust_annotate(
                        input_file    = c(protein_file, protein_file2),
                        cluster_table = data.frame(rep_id = "A", mem_id = "A")
                ),
                regexp = "representative_id"
        )
})

# =============================================================================
# deepclust_annotate() -- realign parameter
# =============================================================================

test_that("deepclust_annotate() with realign returns a tibble with 8 expected columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        realign  <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        result <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters,
                realign       = realign
        )
        expect_true(tibble::is_tibble(result))
        expect_equal(
                colnames(result),
                c("representative_id", "member_id", "file_name",
                  "approx_pident", "evalue", "bitscore", "qcovhsp", "scovhsp")
        )
})

test_that("deepclust_annotate() realign result is non-empty", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        realign  <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        result <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters,
                realign       = realign
        )
        expect_gt(nrow(result), 0L)
})

test_that("deepclust_annotate() realign alignment columns are numeric", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        realign  <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        result <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters,
                realign       = realign
        )
        expect_type(result[["approx_pident"]], "double")
        expect_type(result[["evalue"]],        "double")
        expect_type(result[["bitscore"]],      "double")
        expect_type(result[["qcovhsp"]],       "double")
        expect_type(result[["scovhsp"]],       "double")
})

test_that("deepclust_annotate() realign has same row count as without realign", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        realign  <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        with_realign    <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters,
                realign       = realign
        )
        without_realign <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters
        )
        expect_equal(nrow(with_realign), nrow(without_realign))
})

test_that("deepclust_annotate() realign does not alter representative_id or member_id columns", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        clusters <- get_clusters()
        realign  <- deepclust_realign(
                input_file = c(protein_file, protein_file2),
                clusters   = clusters
        )
        with_realign    <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters,
                realign       = realign
        )
        without_realign <- deepclust_annotate(
                input_file    = c(protein_file, protein_file2),
                cluster_table = clusters
        )
        expect_equal(with_realign[["representative_id"]],
                     without_realign[["representative_id"]])
        expect_equal(with_realign[["member_id"]],
                     without_realign[["member_id"]])
})

test_that("deepclust_annotate() errors when realign is missing required columns", {
        expect_error(
                deepclust_annotate(
                        input_file = c(protein_file, protein_file2),
                        realign    = data.frame(member_id = "A",
                                                representative_id = "A")
                ),
                regexp = "approx_pident"
        )
})
