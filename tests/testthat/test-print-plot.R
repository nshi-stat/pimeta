d <- get_data("sbp")

test_that("print methods", {
  for (m in c("HTS", "HK", "SJ", "KR", "APX")) {
    r <- pima(d$y, d$se, method = m)
    expect_output(print(r), "Prediction & Confidence Intervals")
    expect_output(print(r, trans = "exp"), "exponential")
  }
  for (m in c("DL", "HK", "SJ", "KR", "APX", "PL", "BC")) {
    expect_output(print(cima(d$y, d$se, method = m)), "Confidence Interval")
  }
})

test_that("plot methods", {
  grDevices::pdf(NULL)
  expect_error(plot(pima(d$y, d$se, method = "HTS")), NA)
  expect_error(plot(pima(d$y, d$se, method = "HTS"), trans = "exp", ntick = 4), NA)
  expect_error(plot(cima(d$y, d$se, method = "DL")), NA)
  expect_error(plot(pima(d$y, d$se, method = "HTS"), studylabel = 1:3),
               "same length")
  grDevices::dev.off()
})
