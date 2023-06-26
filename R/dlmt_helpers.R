
# helper functions for dlmt()



# # transform outcomes ------------------------------------------------------
#
# # Monotone spline fit
# mono_spline <- function(lambda, spline_basis) {
#   as.vector(spline_basis %*% lambda[-1:0] + lambda[1])
# }
#
# transform_data <- function(df, outcome=y, i_basis, lambda) {
#   df %>%
#     dplyr::mutate({{outcome}} := mono_spline(lambda, i_basis))
# }



# run dlm with transformations --------------------------------------------

#' Run no-intercept DLM
#'
#' @param df data to run model on
#' @param model_params model hyperparameters w, v0, a0, b0
#' @param cap_var cap maximum variance for players?
#'
#' @return many lists of model results
#'
run_dlmt <- function(
    df,
    outcome = y,
    time = t,
    event = g,
    unit = a,
    model_params,
    model_class = c("multi", "h2h"),
    clear_diags=F,
    cap_var=T) {

  model_class <- match.arg(model_class)
  unit_str <- rlang::as_string(rlang::ensym(unit))

  # some useful constants
  w <- model_params[1]
  v0 <- model_params[2]
  p <- attr(df, "p")

  # define & initialize model values
  as <- vector()
  bs <- vector()
  ms <- list()
  Vinvs <- list()
  Vs <- list()
  Xs <- list()
  ptcps <- list()

  as[1] <- model_params[3]
  bs[1] <- model_params[4]
  ms[[1]] <- rep(0, p)
  Vinvs[[1]] <- 1/v0 * Matrix::Diagonal(p)   # store V-inverses
  Vs[[1]] <- v0 * Matrix::Diagonal(p)        # store Vs

  mint <- min(dplyr::pull(df, {{time}}))
  maxt <- max(dplyr::pull(df, {{time}}))

  # iteratively run model on each time period
  for (tp in mint:maxt) {

    ### filter df and event-size key down to time t
    df_time <- df %>%
      dplyr::filter({{time}} == tp)

    if (model_class == "multi") {
      event_size_key_time <- attr(df, "event_size_key") %>%
        dplyr::filter({{time}} == tp)
    }

    ### compute useful quantities

    # total number of competitors in time t
    n_t <- nrow(df_time)

    if (n_t == 0) {
      browser()
      stop("Haven't implemented empty time periods yet...")
    }

    # X matrix: n_t by p matrix indicating which athlete the scores belong to
    X <- build_Xt(df_time, event={{event}}, unit=unit_str, p=p, model_class=model_class)
    if (model_class == "multi") {
      # centered X matrix: (I - GDG^T)X in our notation
      H <- build_Ht(event_size_key_time, n_t)
      X <- H %*% X
    }
    Xs[[tp]] <- X

    # p by 1 vector indicating which athletes competed at all in time period
    ptcps[[tp]] <- Matrix::colSums(X != 0) != 0

    # X^T X: pxp matrix
    # note: t(X) %*% Xbar is faster than t(Xbar) %*% Xbar,
    #  but code is neater this way
    XtX <- Matrix::t(X) %*% X

    # Xbarty: px1 vector with entries equal to the sum of scores for each athlete
    #  minus the average score for the game
    Xbarty <- as.vector(Matrix::t(X) %*% dplyr::pull(df_time, {{outcome}}))

    # yty is the (scalar) sum of squared scores
    yty <- df_time %>%
      dplyr::summarize(yty = sum({{outcome}}^2)) %>%
      dplyr::pull(yty)


    ### record priors of m, V, a, b

    m_pr <- ms[[tp]]
    a_pr <- as[tp]
    b_pr <- bs[tp]

    # Add innovation variance to prior of V
    if (tp > 1) {
      # Only add innovation variance to athletes who have played a game
      #  - hacky: uses m == 0 to indicate that athlete has never played
      innov <- Matrix::Diagonal(p, x=w*as.vector(!ms[[tp]] == 0))

      # Add innovation variance
      if (clear_diags) {
        V_pr <- Matrix::Diagonal(x=diag(Vs[[tp]])) + innov
      } else {
        # note: we don't add any innovation covariance to off-diagonals
        V_pr <- Vs[[tp]] + innov
      }
      if (cap_var) {
        # cap variance at initial variance v0
        V_pr[V_pr > v0] <- v0
      }

      Vinv_pr <- solve(V_pr)   # annoying that we have to do this: could maybe avoid with tricks
    } else {   # tp == 1
      Vinv_pr <- Vinvs[[tp]]
    }


    ### Update values of m, V, a, b

    Vinv_post <- Vinv_pr + XtX
    V_post <- solve(Vinv_post)   # VERY SLOW
    m_post <- V_post %*% (Vinv_pr %*% m_pr + Xbarty)
    a_post <- a_pr + n_t / 2
    mVm_prior <- Matrix::t(m_pr) %*% Vinv_pr %*% m_pr
    mVm_post <- Matrix::t(m_post) %*% Vinv_post %*% m_post
    b_post <- b_pr + (yty + mVm_prior - mVm_post) / 2


    ### save posterior results for this time period

    ms[[tp]] <- as.vector(m_post)
    Vinvs[[tp]] <- Vinv_post
    Vs[[tp]] <- V_post
    as[tp] <- a_post
    bs[tp] <- b_post

    ### Update priors for next time period
    if (tp < maxt){
      ms[[tp+1]] <- m_post
      Vinvs[[tp+1]] <- Vinv_post
      Vs[[tp+1]] <- V_post
      as[tp+1] <- a_post
      bs[tp+1] <- b_post
    }
  }

  # reduce size of store X matrices; only store columns for players who played in time t
  for (tp in mint:maxt) {
    Xs[[tp]] <- Xs[[tp]][,ptcps[[tp]]]
  }

  # Store results as a tibble
  res <- dplyr::tibble(
    {{time}} := mint:maxt,
    a = as,
    b = bs,
    m = ms,
    V = Vs,
    X = Xs,
    ptcps = ptcps)

  # Record result attributes
  attr(res, "model_params") <- model_params

  res
}


# run smoother

# input: results
# output: smoothed results
smooth_results <- function(res, time=t, w) {
  max_t <- max(dplyr::pull(res, {{time}}))

  # initialize storage lists
  m_smooth_list <- list()
  V_smooth_list <- list()

  # at time T, use final time-period results
  m_smooth_list[[max_t]] <- as.numeric(res$m[[max_t]])
  V_smooth_list[[max_t]] <- res$V[[max_t]]

  for (tp in (max_t-1):1) {
    res_t <- res %>% dplyr::filter({{time}} == tp)

    # only add w if player has started playing by given tp
    w_vec <- rep(w, length(res$m[[1]]))
    w_vec[res_t$m[[1]] == 0] <- 1
    w_mat <- diag(w_vec)
    Vplusw <- res_t$V[[1]] + w_mat

    G_t <- res_t$V[[1]] %*% solve(Vplusw)
    m_smooth_list[[tp]] <- as.numeric(res_t$m[[1]] + G_t %*% (m_smooth_list[[tp+1]] - res_t$m[[1]]))
    V_smooth_list[[tp]] <- res_t$V[[1]] + G_t %*% (V_smooth_list[[tp+1]] - Vplusw ) %*% t(G_t)
  }

  res$m <- m_smooth_list
  res$V <- V_smooth_list
  res$a <- res$a[[max_t]]
  res$b <- res$b[[max_t]]

  return(res)
}


# get predictions from model

#' Get DLM predictive distributions for each event
#'
#' @param df data
#' @param res fitted model
#'
#' @return df with model predictions for each time point
get_predictions <- function(df,
                            res,
                            outcome = y,
                            time = t,
                            event = g,
                            unit = a,
                            model_class = c("multi", "h2h")) {
  model_class <- match.arg(model_class)

  # useful values
  w <- attr(res, "model_params")[1]
  v0 <- attr(res, "model_params")[2]
  p <- attr(df, "p")

  # TODO: these are actually 'ath' values in the original code...
  #  - here we're using the "athleteID"s, which might lead to problems?
  if (model_class == "multi") {
    # gather observed outcomes in appropriate format (one row per time period)
    true_values <- df %>%
      dplyr::group_by({{time}}) %>%
      dplyr::summarize(
        {{event}} := list({{event}}),
        {{unit}} := list({{unit}}),
        {{outcome}} := list({{outcome}}),
        .groups="drop") %>%
      dplyr::mutate(X = res$X,
                    ptcps = res$ptcps)
  } else if (model_class == "h2h") {
    unit_str <- rlang::as_string(rlang::ensym(unit))
    true_values <- df %>%
      dplyr::group_by({{time}}) %>%
      dplyr::summarize(
        {{event}} := list({{event}}),
        "{{unit}}1" := list(.data[[paste0(unit_str, 1)]]),
        "{{unit}}2" := list(.data[[paste0(unit_str, 2)]]),
        {{outcome}} := list({{outcome}}),
        .groups="drop") %>%
      dplyr::mutate(X = res$X %>%
                      purrr::discard(is.null),    # account for empty time periods
                    ptcps = res$ptcps %>%
                      purrr::discard(is.null))
  }

  # make "predictive" values for time 1
  output1 <- dplyr::tibble(
    {{time}} := 1,
    a = attr(res, "model_params")[3],
    b = attr(res, "model_params")[4],
    m = list(rep(0, p)),
    V = list(v0 * Matrix::Diagonal(p)),
  )

  # link each time t+1 with predicted values from time t
  model_outputs <- res %>%
    dplyr::select(-X, -ptcps) %>%
    dplyr::mutate({{time}} := {{time}}+1) %>%
    dplyr::bind_rows(output1)

  # if only one athlete, make sure that X is a matrix
  if (is.null(nrow(true_values$X[[1]]))) {
    true_values <- true_values %>%
      dplyr::rowwise() %>%
      dplyr::mutate(X = list(as.matrix(X)))
  }

  # compute mean and variance estimates for each time period: SUPER slow
  true_values %>%
    dplyr::left_join(model_outputs,
                     by=rlang::as_string(rlang::ensym(time))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pt_est = list(as.vector(X %*% m[ptcps])),
      sigma = list(as.matrix(b/a * (Matrix::Diagonal(nrow(X)) + X %*%
                                      (V[ptcps,ptcps] + w*Matrix::Diagonal(sum(ptcps))) %*%
                                      Matrix::t(X)) )),
      var_theta = list(as.matrix(b/a * X %*% V[ptcps, ptcps] %*% Matrix::t(X)))) %>%
    # computes ONLY diagonal elements of sigma
    #  - faster, but we need full sigma matrix for accurate log-likelihood evaluations
    # sigma = list(as.matrix(b/a * (Diagonal(nrow(X)) +
    #              rowSums( (X %*% (V[ptcps,ptcps] + w*Diagonal(sum(ptcps)))) * X )) ))) %>%
    dplyr::ungroup()
}




# optimize dlm with transformations ---------------------------------------

# note this inputs preprocessed df
# TODO: note that this assumes minimum t is 1! I think it's baked into the model...?
optimize_dlmt <- function(
    df,

    outcome = y,
    time = t,
    event = g,
    unit = a,

    knots,
    i_basis,
    m_basis,

    model_class,
    training_periods = NULL,
    method = c("optim", "mcmc")) {
  method <- match.arg(method)
  p <- attr(df, "p")

  # get number of training periods
  if (is.null(training_periods)) {
    tmax <- max(dplyr::pull(df, {{time}}))
    training_periods <- tmax - round(tmax / 3)
    stopifnot(training_periods < tmax)
  }

  # set training dataset & training spline bases
  training_indices <- which(dplyr::pull(df, {{time}}) <= training_periods)
  df_train <- df[training_indices,]
  tmax_train <- max(dplyr::pull(df_train, {{time}}))
  i_basis_train <- i_basis[training_indices,]
  m_basis_train <- m_basis[training_indices,]

  ### prepare objects for stan

  # get initial parameters for transformation optimization
  init_t_params <- get_ispline_coefs(
    dplyr::pull(df_train, {{outcome}}),
    basis = i_basis_train,
    nonneg = T)

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

    B = ncol(i_basis_train),
    sumB = sum(init_t_params[-1]),
    init_lam = init_t_params,
    i_basis = i_basis_train,
    m_basis = m_basis_train,
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
      verbose = F
    )
  } else {
    browser()
    stop("Haven't implemented MCMC sampling version yet")
    res <- rstan::sampling(
      stanmodels$model_unc,
      data = stan_dat,
      cores = 4,              # number of cores (could use one per chain)
      refresh = 1             # no progress shown
    )
  }

  # add intercept term to output
  res$lam0 <- init_t_params[1]

  res
}




# construct useful objects for rstan optimization of transformation -------

#' Given event-size key for one time period, build H_t matrix
#'
#' @param esk event-size key
#' @param n_t number of players in games within time period
#'
#' @return H_t matrix
build_Ht <- function(esk, n_t) {
  matrix_list <- esk %>%
    dplyr::pull(event_size) %>%
    purrr::map(function(x) matrix(1/x, x, x))

  Matrix::Diagonal(n_t) - Matrix::bdiag(matrix_list)
}


build_X <- function(df, time=t, event=g, unit=a,
                    tmax, p,
                    model_class = c("multi", "h2h")) {
  model_class <- match.arg(model_class)
  unit_str <- rlang::as_string(rlang::ensym(unit))

  1:tmax %>%
    purrr::map(function(x) {
      df_time <- df %>%
        dplyr::filter({{time}} == x)
      n_t <- nrow(df_time)
      if (n_t == 0) { return(tibble()) }

      X <- build_Xt(df_time,
                    event = {{event}},
                    unit = unit_str,
                    p = p,
                    model_class = model_class)

      if (model_class == "multi") {
        event_size_key_time <- df_time %>%
          dplyr::group_by({{event}}) %>%
          dplyr::summarize(event_size = dplyr::n())
        H <- build_Ht(event_size_key_time, n_t)
        X <- H %*% X
      }

      as.matrix(X)
    }) %>%
    do.call(rbind, .)
}

#' Given data for one time period, build X_t matrix
#'
#' @param df data for single time period
#' @param p total number of athletes
#'
#' @return X_t matrix
build_Xt <- function(
    df,
    event=g,
    unit="a",
    p,
    model_class = c("multi", "h2h")) {
  model_class <- match.arg(model_class)

  if (model_class == "multi") {
    X_key <- df %>%
      dplyr::select({{event}}, paste0(unit, "_id")) %>%
      dplyr::arrange({{event}}) %>%
      dplyr::mutate(rowID = dplyr::row_number())
    Matrix::sparseMatrix(
      i = X_key$rowID,
      j = dplyr::pull(X_key, paste0(unit, "_id")),
      dims = c(nrow(df), p))
  } else {
    X <- matrix(0, nrow=nrow(df), ncol=p)
    for (i in 1:nrow(df)) {
      X[i, dplyr::pull(df, paste0(unit, "1_id"))[i]] <- 1
      X[i, dplyr::pull(df, paste0(unit, "2_id"))[i]] <- -1
    }
    X
  }

}

# column numbers of athletes who participate in each time period
build_participants <- function(df, time=t, unit=a,
                               tmax,
                               model_class = c("multi", "h2h")) {
  model_class <- match.arg(model_class)
  unit_str <- rlang::as_string(rlang::ensym(unit))

  if (model_class == "multi") {
    cols <- paste0(unit_str, "_id")
  } else {
    cols <- c(
      paste0(unit_str, "1_id"),
      paste0(unit_str, "2_id")
    )
  }

  1:tmax %>%
    purrr::map(function(x) {
      df %>%
        dplyr::filter({{time}} == x) %>%
        dplyr::select(all_of(cols)) %>%
        unlist(use.names=F) %>%
        unique()
    }) %>%
    unlist()
}

build_nT <- function(df, time=t) {
  df %>%
    dplyr::group_by({{time}}) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::pull(n)
}

build_nT_unique <- function(df, time=t, unit=a,
                            tmax,
                            model_class = c("multi", "h2h")) {
  model_class <- match.arg(model_class)
  unit_str <- rlang::as_string(rlang::ensym(unit))

  if (model_class == "multi") {
    cols <- paste0(unit_str, "_id")
  } else {
    cols <- c(
      paste0(unit_str, "1_id"),
      paste0(unit_str, "2_id")
    )
  }

  1:tmax %>%
    purrr::map_dbl(function(x) {
      df %>%
        dplyr::filter({{time}} == x) %>%
        dplyr::select(all_of(cols)) %>%
        unlist(use.names=F) %>%
        unique() %>%
        length()
    })
}
