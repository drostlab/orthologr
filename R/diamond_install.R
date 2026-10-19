#' @title Install the DIAMOND aligner into the user cache
#' @description Downloads a pre-compiled DIAMOND2 binary from the official
#' GitHub releases (\url{https://github.com/bbuchfink/diamond/releases}) that
#' matches the operating system of the user and stores it in a user-level cache
#' directory (see \code{\link[tools]{R_user_dir}}). No administrator rights are required.
#' Once installed, \code{\link{diamond}}, \code{\link{diamond_best}},
#' \code{\link{diamond_rec}}, \code{\link{set_diamond}}, \code{\link{deepclust}} and
#' \code{\link{deepclust_realign}} find the cached binary automatically whenever
#' \code{path = NULL} and no \code{diamond} is available on the \code{PATH}.
#' @param version a DIAMOND release tag, e.g. \code{"v2.1.9"}. The default \code{"latest"}
#' uses the most recent release.
#' @param force logical; re-download even if DIAMOND is already cached.
#' @param quiet logical; suppress progress messages.
#' @details Pre-compiled binaries are provided upstream for Linux (x86_64), macOS and Windows (x86_64).
#' On other platforms (e.g. Linux ARM) please install DIAMOND via conda:
#' \code{conda install -c bioconda diamond}.
#' @return Invisibly, the directory containing the \code{diamond} executable. This can be passed
#' to the \code{path} argument of the DIAMOND functions.
#' @examples
#' \dontrun{
#' diamond_install()
#' diamond_path()
#' }
#' @seealso \code{\link{diamond_path}}, \code{\link{diamond}}
#' @export
diamond_install <- function(version = "latest", force = FALSE, quiet = FALSE) {
        
        sysname <- Sys.info()[["sysname"]]
        machine <- Sys.info()[["machine"]]
        
        asset <- switch(
                sysname,
                Linux   = if (machine %in% c("x86_64", "amd64")) "diamond-linux64.tar.gz" else NA_character_,
                Darwin  = "diamond-macos.tar.gz",
                Windows = "diamond-windows.zip",
                NA_character_
        )
        if (is.na(asset))
                stop("No pre-compiled DIAMOND binary is available for ", sysname, " (", machine, "). ",
                     "Please install it via conda: 'conda install -c bioconda diamond'.",
                     call. = FALSE)
        
        exe_name <- diamond_exe_name()
        dest     <- diamond_cache_dir()
        exe      <- file.path(dest, exe_name)
        
        if (file.exists(exe) && !force) {
                if (!quiet) message("DIAMOND is already installed in ", dest, " (use force = TRUE to reinstall).")
                return(invisible(dest))
        }
        
        if (identical(version, "latest")) {
                url <- paste0("https://github.com/bbuchfink/diamond/releases/latest/download/", asset)
        } else {
                version <- sub("^v?", "v", version)
                url <- paste0("https://github.com/bbuchfink/diamond/releases/download/", version, "/", asset)
        }
        
        dir.create(dest, recursive = TRUE, showWarnings = FALSE)
        archive <- file.path(tempdir(), asset)
        on.exit(unlink(archive), add = TRUE)
        
        if (!quiet) message("Downloading ", url, " ...")
        old_timeout <- options(timeout = max(600, getOption("timeout")))
        on.exit(options(old_timeout), add = TRUE)
        
        tryCatch(
                utils::download.file(url, archive, mode = "wb", quiet = quiet),
                error = function(e)
                        stop("Could not download DIAMOND from ", url, ": ", conditionMessage(e), call. = FALSE),
                warning = function(w)
                        stop("Could not download DIAMOND from ", url, ": ", conditionMessage(w), call. = FALSE)
        )
        
        if (grepl("\\.zip$", asset)) {
                utils::unzip(archive, exdir = dest)
        } else {
                utils::untar(archive, exdir = dest)
        }
        
        if (!file.exists(exe))
                stop("The downloaded archive did not contain '", exe_name, "'.", call. = FALSE)
        
        Sys.chmod(exe, mode = "0755")
        
        # verify that the binary runs
        out <- tryCatch(
                suppressWarnings(system2(exe, "version", stdout = TRUE, stderr = TRUE)),
                error = function(e) NULL
        )
        if (is.null(out) || !any(grepl("diamond", out, ignore.case = TRUE)))
                stop("DIAMOND was downloaded to ", dest, " but could not be executed. ",
                     "Please install it via conda instead: 'conda install -c bioconda diamond'.",
                     call. = FALSE)
        
        if (!quiet) message("Installed ", out[1], " in ", dest)
        invisible(dest)
}

#' @title Locate the DIAMOND executable used by orthologr
#' @description Returns the directory holding the \code{diamond} executable that
#' \code{orthologr} will use. A \code{diamond} on the \code{PATH} takes precedence
#' (returned as \code{NULL}, i.e. no \code{path} needed); otherwise the copy installed with
#' \code{\link{diamond_install}} is returned.
#' @return \code{NULL} if \code{diamond} is on the \code{PATH}, the cache directory if only the
#' cached copy exists, and an error otherwise.
#' @seealso \code{\link{diamond_install}}
#' @export
diamond_path <- function() {
        if (nzchar(Sys.which("diamond"))) return(NULL)
        dir <- diamond_cache_dir()
        if (file.exists(file.path(dir, diamond_exe_name()))) return(dir)
        stop("DIAMOND was not found. Run diamond_install() or install it via conda: ",
             "'conda install -c bioconda diamond'.", call. = FALSE)
}

diamond_exe_name <- function() {
        if (.Platform$OS.type == "windows") "diamond.exe" else "diamond"
}

diamond_cache_dir <- function() {
        file.path(tools::R_user_dir("orthologr", which = "cache"), "diamond")
}

# internal: if the user gave no path and diamond is not on the PATH, fall back
# to the copy installed by diamond_install()
diamond_resolve_path <- function(path) {
        if (!is.null(path) || nzchar(Sys.which("diamond"))) return(path)
        dir <- diamond_cache_dir()
        if (file.exists(file.path(dir, diamond_exe_name()))) return(dir)
        path
}
