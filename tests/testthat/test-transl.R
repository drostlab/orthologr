
test_that("transl() returns a character string", {
        result <- transl("ATGATG")
        expect_type(result, "character")
        expect_length(result, 1L)
})

test_that("transl() output length equals input length / 3", {
        seq <- "ATGATGATG"  # 9 nt -> 3 aa
        result <- transl(seq)
        expect_equal(nchar(result), nchar(seq) / 3)
})

test_that("transl() correctly translates ATG (Met start codon)", {
        # ATG encodes Methionine (M)
        result <- transl("ATG")
        expect_equal(toupper(result), "M")
})

test_that("transl() correctly translates a known dipeptide", {
        # ATG = M (Met), AAA = K (Lys)
        result <- transl("ATGAAA")
        expect_equal(toupper(result), "MK")
})

test_that("transl() translates a stop codon to '*'", {
        # TAA is a canonical stop codon
        result <- transl("TAA")
        expect_equal(result, "*")
})

test_that("transl() handles a full ORF (start + coding + stop)", {
        # ATG (M) + GGT (G) + TAA (*)
        result <- transl("ATGGGTTAA")
        expect_equal(toupper(result), "MG*")
})

test_that("transl() output matches seqinr::translate() directly", {
        seq <- "ATGAAAGGG"
        expected <- seqinr::c2s(seqinr::translate(seqinr::s2c(seq)))
        expect_equal(transl(seq), expected)
})
