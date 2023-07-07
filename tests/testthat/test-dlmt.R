
test_that("bad names fail as expected", {
  df <- sim_data(ath_per_game = 5, num_games = 5, tmax = 3) %>%
    dplyr::rename(pt_est = y)

  expect_error({
    res <- dlmt(
      df,

      outcome = pt_est,

      centering = c("median"),

      optimize = T,
      optimize_method = c("optim"),

      smooth = T,
      details = T)
  }, "Error in column names: please rename pt_est column to a different name.")
})

test_that("multiplayer model works", {
  df <- sim_data(ath_per_game = 5, num_games = 5, tmax = 3) %>%
    dplyr::rename(ti = t, ga = g, at = a, ou = y)
  res <- dlmt(
    df,

    outcome = ou,
    time = ti,
    event = ga,
    unit = at,

    centering = c("median"),

    optimize = T,
    optimize_method = c("optim"),

    smooth = T,
    details = T)

  expect_equal(nrow(res), 3)
  expect_equal(nrow(attr(res, "predictions")), 5*5*3)
})

test_that("head-to-head model works", {
  df <- sim_data(ath_per_game = 2, num_games = 5, tmax = 3) %>%
    dplyr::rename(ti = t, ga = g, at = a, ou = y)
  res <- dlmt(
    df,

    outcome = ou,
    time = ti,
    event = ga,
    unit = at,

    centering = "none",

    optimize = T,
    optimize_method = c("optim"),

    smooth = T,
    details = T)

  expect_equal(nrow(res), 3)
  expect_equal(nrow(attr(res, "predictions")), 5*3)
})


test_that("model runs without transformations", {
  df <- sim_data(ath_per_game = 5, num_games = 5, tmax = 3)
  res <- dlmt(
    df,
    optimize = F,
    smooth = T,
    details = T)

  expect_equal(nrow(res), 3)
  expect_equal(nrow(attr(res, "predictions")), 5*5*3)
})

test_that("can specify number of knots", {
  df <- sim_data(ath_per_game = 5, num_games = 5, tmax = 3)
  res <- dlmt(
    df,
    optimize = T,
    smooth = F,
    details = T)

  mp <- attr(res, "model_params")
  lam <- attr(res, "lambda")

  res2 <- dlmt(
    df,
    optimize = F,
    lam = lam,
    model_params = mp,
    smooth = T,
    details = T
  )

  expect_equal(attr(res2, "model_params"), mp)
  expect_equal(attr(res2, "lambda"), lam)
})

test_that("h2h model runs without transformations", {
  df <- sim_data(ath_per_game = 2, num_games = 10, tmax = 3)

  expect_warning({
    res <- dlmt(
      df,
      optimize = F,
      smooth = T,
      details = T)})
  expect_equal(nrow(attr(res, "predictions")), nrow(df)/2)
})

test_that("dummy test", {
  df <- sim_data(ath_per_game = 5, num_games = 5, tmax = 3)
  preprocess_data(df)
  expect_equal(1,1)
})



