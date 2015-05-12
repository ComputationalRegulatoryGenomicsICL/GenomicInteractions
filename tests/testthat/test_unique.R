## test uniq/duplicated
test_that("unique works", {
  data(hic_example_data)

  expect_equal(unique(hic_example_data[c(1:4, 2:5)]), hic_example_data[1:5])
})


test_that("duplicated works", {
  data(hic_example_data)
  
  expect_equal(which(duplicated(hic_example_data[c(1:4, 2:5)])),5:7)
})
