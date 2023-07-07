
#' Run dynamic linear model with transformations (Che and Glickman, 2023+)
#'
#' @param df dataframe
#' @param outcome (unquoted) name of outcome column
#' @param time (unquoted) name of time column
#' @param event (unquoted) name of event ID column
#' @param unit (unquoted) name of unit ID column
#' @param centering centering to apply to outcomes from each event: "median", "mean", or "none"
#' @param optimize find an optimal outcome transformation and innovation variance parameter?
#' @param optimize_method method for estimating optimal parameters: "optim" for MAP estimation, "mcmc" for full MCMC estimation
#' @param lam (optional) vector of transformation parameters lambda; default to identity transformation. Note that the given lambda parameters will be used on the I-spline basis associated with the training data in the given df
#' @param model_params vector of model parameters (w, V0, a0, b0): (innovation variance, initial variance, shape and scale parameters for inverse-gamma prior on observation variance)
#' @param smooth smooth posterior latent ability estimates?
#' @param details include visualization of transformation, model predictions, and model estimates in output?
#' @param ... optional parameters, e.g., number of knots, number of training periods
#'
#' @return df with posterior model estimates at each time period, optionally with helpful attributes
#' @export
dlmt <- function(
    df,
    outcome = y,
    time = t,
    event = g,
    unit = a,
    centering = c("median", "mean", "none"),

    optimize = T,
    optimize_method = c("optim", "mcmc"),
    lam = NULL,
    model_params = NULL,

    smooth = F,
    details = F,
    ...) {

  args <- list(...)
  centering <- match.arg(centering)
  optimize_method <- match.arg(optimize_method)

  check_names({{outcome}})
  check_names({{time}})
  check_names({{event}})
  check_names({{unit}})

  # automatically determine if model should be h2h or multiplayer
  max_players <- df %>%
    dplyr::group_by({{event}}) %>%
    dplyr::summarize(num_units = length(unique({{unit}}))) %>%
    dplyr::summarize(max_players = max(num_units)) %>%
    dplyr::pull(max_players)
  if (max_players == 2) {
    cat("Using head-to-head model...")
    model_class <- "h2h"
    df <- to_h2h(df,
                 outcome = {{outcome}},
                 time = {{time}},
                 event = {{event}},
                 unit = {{unit}})
    if (centering != "none") {
      warning("For head-to-head data, centering should be \"none\": manually fixing")
      centering <- "none"
    }
  } else {
    cat("Using multiplayer model...")
    model_class <- "multi"
  }

  # preprocess dataset
  dfp <- df %>%
    preprocess_data(
      outcome = {{outcome}},
      time = {{time}},
      event = {{event}},
      unit = {{unit}},
      event_center = centering,
      model_class = model_class)

  # build spline basis (if optimizing or given a transformation)
  if (optimize | !is.null(lam)) {
    # prepare monotone spline transformation
    knots <- get_knots(
      dfp,
      outcome = {{outcome}},
      num_knots = args$num_knots)
    i_basis <- get_ispline_basis(
      dplyr::pull(dfp, {{outcome}}),
      knots = knots,
      degree = 3)
    m_basis <- get_mspline_basis(
      dplyr::pull(dfp, {{outcome}}),
      knots = knots,
      degree = 3)
  }

  # set default model parameters, if not specified
  if (is.null(model_params)) {
    model_params <- c(0.5, 10, 2, sd(dplyr::pull(df, {{outcome}})))
    using_default_mp <- T
  } else {
    using_default_mp <- F
  }

  # if optimizing, get optimized transformation and w
  if (optimize) {
    stopifnot(!is.null(optimize_method))
    if (!is.null(lam) | !using_default_mp) {
      warning("Using optimized lambda and w values instead of user-inputted lambda and w values")
    }
    opt <- optimize_dlmt(
      df = dfp,
      outcome = {{outcome}},
      time = {{time}},
      event = {{event}},
      unit = {{unit}},

      knots = knots,
      i_basis = i_basis,
      m_basis = m_basis,

      training_periods = args$training_periods,
      model_class = model_class,
      method = optimize_method,
      verbose = args$verbose
    )
    lam <- c(opt$lam0, opt$theta_tilde[-1])   # theta_tilde without w
    model_params[1] <- opt$theta_tilde[1]
  }

  # transform dataset
  if (optimize | !is.null(lam)) {
    dft <- dfp %>%
      dplyr::mutate({{outcome}} := as.vector(i_basis %*% lam[-1:0] + lam[1]))
  } else {
    dft <- dfp
  }

  # run dlm with transformations
  res <- run_dlmt(
    dft,
    outcome = {{outcome}},
    time = {{time}},
    event = {{event}},
    unit = {{unit}},
    model_params = model_params,
    model_class = model_class)
  if (optimize | !is.null(lam)) {
    attr(res, "lambda") <- lam
  }

  if (smooth) {
    res <- smooth_results(res, time={{time}}, w=model_params[1])
  }

  if (details) {
    # visualize transformation
    if (optimize | !is.null(lam)) {
      yfit <- as.vector(lam[1] + i_basis %*% lam[-1])
    } else {
      # visualize identity transformation
      yfit <- dplyr::pull(dfp, {{outcome}})
    }
    transformation_plot <- dplyr::tibble(
      x    = dplyr::pull(dfp, {{outcome}}),
      yfit = yfit
    ) %>%
      ggplot2::ggplot(ggplot2::aes(x=x, y=yfit)) +
      ggplot2::geom_line(linewidth=1) +
      ggplot2::geom_abline(lty = "dotted") +
      ggplot2::labs(y = "Transformed outcome",
                    x = "Outcome") +
      ggplot2::theme_classic()

    # store predictions for each game
    preds <- get_predictions(df = dft,
                             res = res,
                             outcome = {{outcome}},
                             time = {{time}},
                             event = {{event}},
                             unit = {{unit}},
                             model_class = model_class)

    # store posterior parameter estimates
    maxt <- max(dplyr::pull(df, {{time}}))
    m_df <- attr(dft, "unit_key") %>%
      cbind(purrr::reduce(res$m_posts, cbind)) %>%
      dplyr::as_tibble(.name_repair = "minimal")
    names(m_df) <- c("unit", "unit_id", paste0("t", 1:maxt))

    V_df <- attr(dft, "unit_key") %>%
      cbind(purrr::map(res$V_posts, diag) %>%
              purrr::reduce(cbind)) %>%
      dplyr::as_tibble(.name_repair = "minimal")
    names(V_df) <- c("unit", "unit_id", paste0("t", 1:maxt))

    sig_df <- tail(res, n=1) %>%
      dplyr::select(a_posts, b_posts) %>%
      dplyr::mutate(mean = b_posts/(a_posts-1),
                    mode = b_posts/(a_posts+1))

    attr(res, "transformation_plot") <- transformation_plot
    attr(res, "predictions") <- preds
    attr(res, "m_values") <- m_df
    attr(res, "V_values") <- V_df
    attr(res, "sig_values") <- sig_df
  }

  res
}


