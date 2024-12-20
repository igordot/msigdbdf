test_that("msigdb_sqlite()", {
  expect_error(msigdb_sqlite())
  expect_error(msigdb_sqlite("x"))
})
