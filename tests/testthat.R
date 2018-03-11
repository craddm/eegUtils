library(testthat)
library(eegUtils)
Sys.setenv("R_TESTS" = "")
test_check("eegUtils")
