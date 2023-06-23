
# function to simulate fake time-series data

#' Given true parameters, generate synthetic data according to DLM model
#'
#' @param W ratio of innovation variance to observation variance
#' @param V0 ratio of intial player ability variance to observation variance
#' @param sd observation standard deviation
#' @param inv_trans_fun (inverse of) transformation to apply to data
#' @param num_athletes total number of athletes
#' @param ath_per_game number of athletes that play each game
#' @param num_games number of games per time period
#' @param tmax number of time periods
#'
#' @return a synthetic dataset
#'
sim_data <- function(
    W = 0.5,
    V0 = 10,
    sd = 5,
    inv_trans_fun = function(x) {x},
    num_athletes = 50,
    ath_per_game = 10,
    num_games = 25,
    tmax = 10,
    HX = T) {   # generate scores around centered thetas?

  if (ath_per_game > num_athletes) {
    stop("Number of athletes per game (ath_per_game) must be less than the total number of athletes in the dataset (num_athletes)")
  }

  # Athlete IDs: 1 to p
  as <- 1:num_athletes

  # Generate initial athlete skills:
  #   rnorm(0, sd * V0)
  theta0 <- dplyr::tibble(
    a = as,
    theta = rnorm(num_athletes, mean=0, sd=sqrt(sd^2 * V0))
  )

  df_list <- list()
  theta_list <- list()
  theta_list[[1]] <- theta0

  for(t in 1:tmax) {

    # Update true athlete thetas via normal random walk
    if (t > 1) {
      theta_t <- rnorm(num_athletes,
                       mean=theta_list[[t-1]]$theta,
                       sd=sqrt(sd^2 * W))

      theta_t_df <- dplyr::tibble(
        a = as,
        theta = theta_t
      )
      theta_list[[t]] <- theta_t_df
    }

    # Generate observed scores for all games in time t
    game_list <- list()
    for (game in 1:num_games) {
      # randomly sample players to play in game g
      ath_playing <- sample(as,
                            ath_per_game,
                            replace=F)
      ath_thetas <- theta_list[[t]]$theta

      # generate scores for athletes
      if (HX) {
        # psi ~ N(Xbar * theta, sigma^2)
        scores <- rnorm(ath_per_game,
                        mean = ath_thetas[ath_playing] -
                          mean(ath_thetas[ath_playing]),
                        sd = sd)
      } else {
        scores <- rnorm(ath_per_game,
                        mean = ath_thetas[ath_playing],
                        sd = sd)
      }


      # store game outcome
      game_outcome <- dplyr::tibble(
        a = ath_playing,
        theta = ath_thetas[ath_playing],
        t_y = scores,
        y = inv_trans_fun(scores))
      game_list[[game]] <- game_outcome
    }


    # Record game
    df_list[[t]] <- game_list %>%
      dplyr::bind_rows(.id = "g") %>%
      dplyr::mutate(t = t,
                    g = as.numeric(g) + num_games*(t-1)) %>%
      dplyr::select(-theta)
  }

  # Store simulated outcome data, along with true athlete thetas at each time t
  thetas <- plyr::ldply(theta_list) %>%
    dplyr::mutate(t = rep(1:tmax, each=num_athletes))
  df <- df_list %>%
    dplyr::bind_rows() %>%
    dplyr::left_join(thetas, by=c("a", "t")) %>%
    dplyr::select(t, g, a, theta, t_y, y)

  # Store attributes of synthetic dataset
  synth_params <- c(W, V0, sd)
  names(synth_params) <- c("W", "V0", "sd")
  attr(df, "params") <- synth_params
  synth_sampsize <- c(num_athletes, ath_per_game, num_games, tmax)
  names(synth_sampsize) <- c("num_athletes", "ath_per_game", "num_games", "tmax")
  attr(df, "size") <- synth_sampsize
  attr(df, "thetas") <- thetas

  return(df)
}


# Plot how simulated data look
if (FALSE) {

  source("data_generating_process.R")

  W <- 0.5
  V0 <- 10
  sigma <- 10
  num_athletes <- 100
  ath_per_game <- 10
  num_games <- 5
  tmax <- 25

  synth_data <- sim_data(
    W=W, V0=V0, sd=sigma,
    num_athletes=num_athletes,
    ath_per_game=ath_per_game,
    num_games=num_games,
    tmax=tmax)

  # visualize thetas for a few athletes
  synth_data %>%
    dplyr::filter(a %in% sample(unique(a), 6)) %>%
    ggplot2::ggplot(ggplot2::aes(x=t, color=factor(a))) +
    ggplot2::geom_line(ggplot2::aes(y=theta)) +
    ggplot2::facet_wrap(~a) +
    ggplot2::labs(y="") +
    ggplot2::guides(color="none") +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = '#FDFAF5'),
      plot.background = ggplot2::element_rect(fill = '#FDFAF5'))
}
