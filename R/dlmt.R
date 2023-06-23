
# main function for package; runs DLM with transformations

# TODO: optional inputs
#  - knot locations
#  - optimization method (MAP, MCMC)
# output: posterior somethings

# TODO: note that this assumes minimum t is 1! I think it's baked into the model...?
dlmt <- function(
    df,

    outcome = y,
    time = t,
    event = g,
    unit = a,

    centering = c("none", "median", "mean"),
    training_periods = NULL,
    model_class = NULL,
    method = c("optim", "mcmc")) {

  centering <- match.arg(centering)
  method <- match.arg(method)
  p <- length(unique(dplyr::pull(df, {{unit}})))

  # browser()

  # set number of training periods
  if (is.null(training_periods)) {
    tmax <- max(dplyr::pull(df, {{time}}))
    training_periods <- tmax - round(tmax / 3)
  }
  stopifnot(training_periods < tmax)

  # automatically determine if model should be h2h or multiplayer
  if (is.null(model_class)) {
    max_players <- df %>%
      dplyr::group_by({{event}}) %>%
      dplyr::summarize(num_units = length(unique({{unit}}))) %>%
      dplyr::summarize(max_players = max(num_units)) %>%
      dplyr::pull(max_players)

    if (max_players == 2) {
      model_class <- "h2h"
      df <- to_h2h(df)
    } else {
      model_class <- "multi"
    }
  }

  # preprocess dataset
  df <- df %>%
    preprocess_data(event_center = centering,
                    model_class = model_class)

  # make training dataset
  df_train <- df %>%
    dplyr::filter({{time}} <= training_periods)
  tmax_train <- max(dplyr::pull(df_train, {{time}}))

  # prepare monotone spline transformation
  training_knots <- get_knots(df_train, outcome = {{outcome}}, num_knots = 3)
  i_basis <- get_ispline_basis(
    dplyr::pull(df_train, {{outcome}}),
    knots = training_knots,
    degree = 3)
  m_basis <- get_mspline_basis(
    dplyr::pull(df_train, {{outcome}}),
    knots = training_knots,
    degree = 3)
  init_t_params <- get_ispline_coefs(
    dplyr::pull(df_train, {{outcome}}),
    basis = i_basis,
    nonneg = T)

  ### prepare objects for stan

  # build X matrices for model
  X <- build_X(df_train,
               time={{time}}, event={{event}}, unit={{unit}},
               tmax = tmax_train, p = p,
               model_class = model_class)
  # record which units participate in each time period
  participants <- build_participants(df_train,
                                     time={{time}}, unit={{unit}},
                                     tmax = tmax_train,
                                     model_class = model_class)
  # record number of (unique) units in each time period
  nT <- build_nT(df_train, time={{time}})
  nT_unique <- build_nT_unique(df_train,
                               time={{time}}, unit={{unit}},
                               tmax = tmax_train,
                               model_class = model_class)
  stan_dat <- list(
    p = p,
    tmax = tmax_train,

    nT = nT,
    nT_unique = nT_unique,
    sumnT = sum(nT),
    sumnT_unique = sum(nT_unique),

    B = ncol(i_basis),
    sumB = sum(init_t_params[-1]),
    init_lam = init_t_params,
    i_basis = i_basis,
    m_basis = m_basis,
    alpha = sd(dplyr::pull(df_train, {{outcome}})),

    y = dplyr::pull(df_train, {{outcome}}),
    Xbar = X,
    ptcps = participants
  )

  ### run stan model!

  if (method == "optim") {
    res <- rstan::optimizing(
      stanmodels$model_unc,
      data = stan_dat,
      hessian = F,
      # hessian = T,
      # draws = 1000,
      # importance_resampling=T,   # get approximate posterior draws?
      # init = list(lam_scale = sum(init_t_params[-1])),   # initialize lam_scale at the right spot so things don't break
      verbose = T
    )
  } else {
    res <- rstan::sampling(
      stanmodels$model_unc,
      data = stan_dat,
      cores = 4,              # number of cores (could use one per chain)
      refresh = 1             # no progress shown
    )
  }

  res
}

if (F) {
  df <- sim_data(ath_per_game = 5)
  res <- dlmt(df, method="optim")

  df <- sim_data(ath_per_game = 2)
  res <- dlmt(df, method="optim")
}



