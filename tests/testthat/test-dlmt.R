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

    model_params = NULL,
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

    model_params = NULL,
    optimize = T,
    optimize_method = c("optim"),

    smooth = T,
    details = T)

  expect_equal(nrow(res), 3)
  expect_equal(nrow(attr(res, "predictions")), 5*3)
})
