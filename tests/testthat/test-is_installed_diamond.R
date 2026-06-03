context("Test: is_installed_diamond()")

diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))

test_that("is_installed_diamond() returns TRUE when DIAMOND2 is on PATH", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        expect_true(is_installed_diamond())
})

test_that("is_installed_diamond() returns TRUE when a valid exec path is provided", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        # Resolve the actual path diamond lives on so we can pass it explicitly
        diamond_path <- dirname(Sys.which("diamond"))
        skip_if(nchar(diamond_path) == 0, "Could not resolve diamond PATH")
        expect_true(is_installed_diamond(diamond_exec_path = diamond_path))
})
