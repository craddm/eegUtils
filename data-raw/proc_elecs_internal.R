biosemi64 <- import_chans("data-raw/biosemi64.txt")
biosemi64alpha <- import_chans("data-raw/biosemi64alpha.txt")
biosemi128 <- import_chans("data-raw/biosemi128.txt")
biosemi256 <- import_chans("data-raw/biosemi256.txt")
biosemi16 <- import_chans("data-raw/biosemi16.txt")
biosemi32 <- import_chans("data-raw/biosemi32.txt")
electrodeLocs <- import_chans("data-raw/standard_1005.elc")
usethis::use_data(biosemi64,
                  biosemi64alpha,
                  biosemi128,
                  biosemi256,
                  biosemi16,
                  biosemi32,
                  electrodeLocs,
                  internal = TRUE,
                  overwrite = TRUE)
