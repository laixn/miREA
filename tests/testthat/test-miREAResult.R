test_that("miREAResult constructor and accessors work", {
  res <- miREAResult(
    result = list(TG_ORA = data.frame(pathway = "PW1", padj = 0.01)),
    time = data.frame(method = "TG_ORA", time = 0.1),
    valid_methods = "TG_ORA",
    parameter = list(n_pathway = 1)
  )

  expect_s4_class(res, "miREAResult")
  expect_equal(res$valid_methods, "TG_ORA")
  expect_equal(res[["valid_methods"]], "TG_ORA")
  expect_equal(res$time$time, 0.1)
  expect_equal(res[["parameter"]]$n_pathway, 1)
  expect_output(print(res), "miREAResult object")

  df <- as.data.frame(res)
  expect_equal(df$method, "TG_ORA")
  expect_equal(nrow(df), 1)
  expect_equal(res$failed_methods, list())
})

test_that("miREAResult setValidity() rejects inconsistent slots", {
  expect_error(
    miREAResult(
      result = list(TG_ORA = data.frame(pathway = "PW1", padj = 0.01)),
      time = data.frame(method = "TG_ORA", time = 0.1),
      valid_methods = c("TG_ORA", "MiR_ORA"),
      parameter = list(n_pathway = 1)
    ),
    "valid_methods must be a subset of names\\(result\\)"
  )

  expect_error(
    miREAResult(
      result = list(TG_ORA = data.frame(pathway = "PW1", padj = 0.01)),
      time = data.frame(method = "MiR_ORA", time = 0.1),
      valid_methods = "TG_ORA",
      parameter = list(n_pathway = 1)
    ),
    "time\\$method must be a subset of valid_methods"
  )
})

test_that("miREAResult $ and [[ error clearly on unknown slot names", {
  res <- miREAResult(
    result = list(TG_ORA = data.frame(pathway = "PW1", padj = 0.01)),
    time = data.frame(method = "TG_ORA", time = 0.1),
    valid_methods = "TG_ORA",
    parameter = list(n_pathway = 1)
  )

  expect_error(res$not_a_slot, "is not a valid miREAResult slot")
  expect_error(res[["not_a_slot"]], "is not a valid miREAResult slot")
})

test_that("miREAResult show() reports failed methods when present", {
  res <- miREAResult(
    result = list(TG_ORA = data.frame(pathway = "PW1", padj = 0.01)),
    time = data.frame(method = "TG_ORA", time = 0.1),
    valid_methods = "TG_ORA",
    failed_methods = list(MiR_ORA = "boom"),
    parameter = list(n_pathway = 1)
  )
  expect_equal(res$failed_methods, list(MiR_ORA = "boom"))
  expect_output(print(res), "Failed methods: MiR_ORA")
})
