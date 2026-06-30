get_filename <- function(file_path) {
        tools::file_path_sans_ext(basename(file_path))
}
