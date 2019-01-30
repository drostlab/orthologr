## `orthologr` version  0.2.0

### New Features

- new function `retrieve_longest_isoforms()` which enables retrieval of the longest isoforms from a proteome file and save results as fasta file for downstream analyses
- new function `OF2CoreOrthologs()` to retrieve core orthologs across multiple species from Orthofinder2 output
- new function `extract_features()`: Helper function to extract gene loci and splice variant IDs from GFF files
- new function `filter_best_hits()`: Helper function to select best BLAST hit based on minimum evalue
- new function `generate_ortholog_tables()`: Generate ortholog tables by gene locus and splice varaint

## `orthologr` version  0.1.0

### New Features

- `read.cds()` now trimms corrupted CDS (= CDS not divisible by 3) when `delete_corrupt_cds = FALSE` is specified
- the new default value of argument `delete_corrupt_cds` in `dNdS()` is now `FALSE`.
Thus, given the new trimming feature in `read.cds()`, corrupted CDS equences will be trimmed before being translated.


## `orthologr` version  0.0.5

### Bug fixes

- The default setting of the `BLAST` argument `max_target_seqs 1` was removed from `blast_best()` and `blast_rec()` due to
the misunderstood functionality of the `BLAST` argument (See details [here](https://watermark.silverchair.com/bty833.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAokwggKFBgkqhkiG9w0BBwagggJ2MIICcgIBADCCAmsGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMvrRsd_dqJkkrO5AiAgEQgIICPNCvHwDI6PCzjXjP4FtOO2O0HPaNbVmBBOW7I0YBwK6an33k32Zm3d24U1pB5TsoPAmqK1TBBz8y9IEEynMl7IQKmupg3MXTiUqUe0oLrMDIrL4Szunm255RdiIlvq5RDR5V3B7ejVCxTbwbmvbggBHyiMCcvEEufaNROT8z7XXaUiHx3OyD1HaRkXIM4MC5x_b734rkKFIploMWjLSX5KZS1i2prd_Shn2seyV4E2cj3RVw9YzcBEZvRHpbhzWIOWSjXR3ZqYTlDCxzk5CO0d15J57u3C64sM1rQRo2Atny5KR5VL4oLxy8bBBGsQl-k6lQavisoshkBZPfLMr5kEHMrCmCtw2Q2A6pzgYC-jJ7udunpA7DK2GZsSQc-ApYXHFdEOovrtFYNHFrqlE8IsooVB7U89re0lMzZM4gYoMzcDuvLFQ1_msIl51QMmPi8yIY6E5DA2h1Ho93wUxekwTBdxB3ovAtpz4_J93w8cn4z2EK9ygHoxLRQ90kfd7nawSMTxFwBxDB0iK3_w0YtnObzcyd_JAWZM9fobNHDc1CY9j4DWKclIxhqHeCHKMtw9nB9CSa7Bccbf1EBHw4uf1-qYl30IdUsKQUpP5YOtKA7nUDTpTJiKKEf1V3MmW6rPreBayi45-QPk7GkhlrgJgaXdU4LMSGxjaDdcXhOX3W35Y95VSAZEf3rlueplgDTRnRxu8NNZ8S1c5ccA1vP0gFth80NTwzUzZy2kSITIndok_0MWSeuavaRQ7g) and #9 ; Many thanks to @armish)

## `orthologr` version  0.0.4

### New Features 

- Users can now control via the new `delete_corrupt_cds` argument in `dNdS()` and related downstream functions
whether or not corrupted input coding sequences shall be removed prior to dN/dS inference.
In case corrupted CDS exist, the `dNdS()` function will now generate a separate fasta file which
stores all corrupted CDS so that they can be investigated. See issue #8 for details.

#### Function updates

- `dNdS()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `read.cds()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `cds2aa()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `set_blast()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast_best()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast_rec()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

## `orthologr` version  0.0.3

- Fixing internal path bug that caused that wrong pal2nal paths were generated when using multiple sequence aligners -> see issue https://github.com/HajkD/orthologr/issues/5 (Many thanks to Dr. Mario López-Pérez)

## `orthologr` version  0.0.2

- Fixing a major bug that caused KaKs_Calculator to not be able to correctly 
parse the kaks computation output (Many thanks to [Hongyi Li](https://github.com/lihongyi123) who spotted the bug and found a solution).
