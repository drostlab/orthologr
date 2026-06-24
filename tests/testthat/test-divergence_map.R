
mock_dnds_tbl <- tibble::tibble(
        query_id   = paste0("AT", 1:10),
        subject_id = paste0("AL", 1:10),
        dNdS       = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95),
        dN         = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10),
        dS         = c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65)
)

test_that("divergence_map() returns a data frame ...", {
        result <- divergence_map(mock_dnds_tbl)
        expect_true(is.data.frame(result))
})

test_that("divergence_map() returns divergence_strata and query_id columns ...", {
        result <- divergence_map(mock_dnds_tbl)
        expect_true(all(c("divergence_strata", "query_id") %in% colnames(result)))
})

test_that("divergence_map() includes subject_id with subject.id = TRUE ...", {
        result <- divergence_map(mock_dnds_tbl, subject.id = TRUE)
        expect_true("subject_id" %in% colnames(result))
})

test_that("divergence_map() excludes subject_id with subject.id = FALSE ...", {
        result <- divergence_map(mock_dnds_tbl, subject.id = FALSE)
        expect_false("subject_id" %in% colnames(result))
})

test_that("divergence_map() preserves row count ...", {
        result <- divergence_map(mock_dnds_tbl, subject.id = FALSE)
        expect_equal(nrow(result), nrow(mock_dnds_tbl))
})

test_that("divergence_map() returns deciles by default ...", {
        result <- divergence_map(mock_dnds_tbl, n_quantile = 10)
        n_unique <- length(unique(result[["divergence_strata"]]))
        expect_equal(n_unique, 10)
})

test_that("divergence_map() returns quintiles with n_quantile = 5 ...", {
        result <- divergence_map(mock_dnds_tbl, n_quantile = 5)
        n_unique <- length(unique(result[["divergence_strata"]]))
        expect_equal(n_unique, 5)
})

test_that("divergence_map() strata are numeric ...", {
        result <- divergence_map(mock_dnds_tbl)
        expect_true(is.numeric(result[["divergence_strata"]]))
})

test_that("divergence_map() strata form sequential integers ...", {
        result <- divergence_map(mock_dnds_tbl, n_quantile = 10)
        strata <- sort(unique(result[["divergence_strata"]]))
        expect_equal(strata, 1:length(strata))
})

test_that("divergence_map() preserves query_ids ...", {
        result <- divergence_map(mock_dnds_tbl, subject.id = FALSE)
        expect_equal(sort(unique(result[["query_id"]])), sort(unique(mock_dnds_tbl[["query_id"]])))
})

test_that("divergence_map() preserves subject_ids with subject.id = TRUE ...", {
        result <- divergence_map(mock_dnds_tbl, subject.id = TRUE)
        expect_equal(sort(unique(result[["subject_id"]])), sort(unique(mock_dnds_tbl[["subject_id"]])))
})

test_that("divergence_map() errors on low variance dNdS values ...", {
        low_var <- tibble::tibble(
                query_id   = paste0("G", 1:5),
                subject_id = paste0("S", 1:5),
                dNdS       = rep(0.5, 5),
                dN         = rep(0.01, 5),
                dS         = rep(0.02, 5)
        )
        expect_error(divergence_map(low_var, n_quantile = 10))
})

test_that("divergence_map() works with data.table input ...", {
        dt_input <- data.table::as.data.table(mock_dnds_tbl)
        result <- divergence_map(dt_input)
        expect_true(is.data.frame(result))
        expect_equal(nrow(result), nrow(mock_dnds_tbl))
})

test_that("divergence_map() assigns lower strata to lower dNdS values ...", {
        result <- divergence_map(mock_dnds_tbl, n_quantile = 10)
        min_idx <- which.min(mock_dnds_tbl[["dNdS"]])
        max_idx <- which.max(mock_dnds_tbl[["dNdS"]])
        min_strata <- result[min_idx, "divergence_strata", drop = TRUE]
        max_strata <- result[max_idx, "divergence_strata", drop = TRUE]
        expect_lt(min_strata, max_strata)
})
