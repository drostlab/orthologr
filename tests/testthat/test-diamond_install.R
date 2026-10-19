# ---------------------------------------------------------------------------
# Offline tests: use a fake cache directory (R_USER_CACHE_DIR) and a fake
# 'diamond' file, so no download happens. Safe for GitHub Actions / CRAN.
# ---------------------------------------------------------------------------

# run `code` with a fresh, empty user cache and (optionally) an empty PATH
with_fake_cache <- function(code, empty_path = FALSE) {
        old <- Sys.getenv(c("R_USER_CACHE_DIR", "PATH"), unset = NA)
        on.exit({
                for (v in names(old)) {
                        if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
                }
        }, add = TRUE)
        fake_cache <- tempfile("fake_cache_")
        dir.create(fake_cache)
        on.exit(unlink(fake_cache, recursive = TRUE), add = TRUE)
        Sys.setenv(R_USER_CACHE_DIR = fake_cache)
        if (empty_path) Sys.setenv(PATH = fake_cache)
        force(code)
}

place_fake_diamond <- function() {
        dir <- orthologr:::diamond_cache_dir()
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        file.create(file.path(dir, orthologr:::diamond_exe_name()))
        dir
}

test_that("diamond_cache_dir() lives in the R user cache", {
        with_fake_cache({
                dir <- orthologr:::diamond_cache_dir()
                dir.create(dir, recursive = TRUE)
                expect_true(startsWith(normalizePath(dir),
                                       normalizePath(Sys.getenv("R_USER_CACHE_DIR"))))
                expect_match(dir, "orthologr")
        })
})

test_that("diamond_install() does nothing if DIAMOND is already cached", {
        with_fake_cache({
                dir <- place_fake_diamond()
                expect_message(res <- diamond_install(), "already installed")
                expect_equal(res, dir)
                expect_silent(diamond_install(quiet = TRUE))
        })
})

test_that("diamond_resolve_path() respects a user supplied path", {
        with_fake_cache({
                place_fake_diamond()
                expect_equal(orthologr:::diamond_resolve_path("/some/path"), "/some/path")
        }, empty_path = TRUE)
})

test_that("diamond_resolve_path() falls back to the cached copy", {
        with_fake_cache({
                expect_null(orthologr:::diamond_resolve_path(NULL))
                dir <- place_fake_diamond()
                expect_equal(orthologr:::diamond_resolve_path(NULL), dir)
        }, empty_path = TRUE)
})

test_that("diamond_path() finds the cached copy or errors", {
        with_fake_cache({
                expect_error(diamond_path(), "diamond_install")
                dir <- place_fake_diamond()
                expect_equal(diamond_path(), dir)
        }, empty_path = TRUE)
})

# ---------------------------------------------------------------------------
# Online tests: download the real binary from GitHub. They only run
# interactively / on a local machine (devtools::test(), devtools::check()),
# never on CRAN or CI, because they need network access and a supported OS.
# ---------------------------------------------------------------------------

skip_if_no_download <- function() {
        skip_on_cran()
        skip_on_ci()
        skip_if_offline("github.com")
        info <- Sys.info()
        supported <- info[["sysname"]] %in% c("Darwin", "Windows") ||
                (info[["sysname"]] == "Linux" && info[["machine"]] %in% c("x86_64", "amd64"))
        skip_if_not(supported, "No pre-compiled DIAMOND binary for this platform")
}

test_that("diamond_install() downloads a working DIAMOND binary", {
        skip_if_no_download()
        with_fake_cache({
                dir <- suppressMessages(diamond_install(quiet = TRUE))
                exe <- file.path(dir, orthologr:::diamond_exe_name())
                expect_true(file.exists(exe))
                out <- system2(exe, "version", stdout = TRUE, stderr = TRUE)
                expect_match(out[1], "diamond version")
                # the DIAMOND wrappers now find it even without diamond on the PATH
                expect_true(is_installed_diamond(diamond_exec_path = dir))
        })
})

test_that("diamond_install(force = TRUE) re-downloads", {
        skip_if_no_download()
        with_fake_cache({
                dir <- place_fake_diamond()
                exe <- file.path(dir, orthologr:::diamond_exe_name())
                expect_equal(file.size(exe), 0)
                suppressMessages(diamond_install(force = TRUE, quiet = TRUE))
                expect_gt(file.size(exe), 0)
        })
})

test_that("diamond_install() fails clearly for a non-existing version", {
        skip_if_no_download()
        with_fake_cache({
                expect_error(suppressWarnings(diamond_install(version = "v0.0.0", quiet = TRUE)),
                             "Could not download DIAMOND")
        })
})
