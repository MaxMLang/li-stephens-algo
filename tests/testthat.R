library(testthat)
lapply(list.files("./src", full.names = T), source)

test_dir(test_path())

