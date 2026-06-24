
# Build a minimal mock dNdS table matching the columns filter_dNdS() uses
mock_dnds <- tibble::tibble(
        query_id   = paste0("gene", 1:6),
        subject_id = paste0("subj", 1:6),
        dN         = c(0.01, 0.02, NA,   0.05, 0.10, 0.00),
        dS         = c(0.10, NA,   0.30, 0.20, 0.50, 0.00),
        dNdS       = c(0.10, NA,   NA,   0.25, 0.20, NA)
)

test_that("filter_dNdS() returns a data frame / tibble", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        expect_true(is.data.frame(result))
})

test_that("filter_dNdS() removes rows with NA in dN", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        expect_false(any(is.na(result[["dN"]])))
})

test_that("filter_dNdS() removes rows with NA in dS", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        expect_false(any(is.na(result[["dS"]])))
})

test_that("filter_dNdS() removes rows where dNdS > threshold", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 0.15)
        expect_true(all(result[["dNdS"]] <= 0.15))
})

test_that("filter_dNdS() retains rows that pass all filters", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 2)
        # gene1: dN=0.01, dS=0.10, dNdS=0.10 -> should be retained
        expect_true("gene1" %in% result[["query_id"]])
})

test_that("filter_dNdS() with a very low threshold removes high dNdS rows", {
        result <- filter_dNdS(mock_dnds, dnds.threshold = 0.05)
        expect_true(nrow(result) < nrow(mock_dnds))
})

test_that("filter_dNdS() with threshold = Inf retains all non-NA rows", {
        result_all_na_removed <- dplyr::filter(mock_dnds,
                                               !is.na(dN), !is.na(dS), !is.na(dNdS))
        result <- filter_dNdS(mock_dnds, dnds.threshold = Inf)
        expect_equal(nrow(result), nrow(result_all_na_removed))
})
