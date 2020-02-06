check_ram <- function() {
        mem <- benchmarkme::get_ram()
        if (is.na(mem)) {
                return (NA)
        } else {
                res <- as.numeric(mem) / 1000000000
                names(res) <- "ram_in_gb"
                return (res)
        }
}