
# helper functions for dlmt()

#####
# construct useful matrices for DLM
#####

#' Given data for one time period, build X_t matrix
#'
#' @param df data for single time period
#' @param p total number of athletes
#'
#' @return X_t matrix
build_Xt <- function(df, event=g, unit=a, p) {
  X_key <- df %>%
    dplyr::select({{event}}, {{unit}}) %>%
    dplyr::arrange({{event}}) %>%
    dplyr::mutate(rowID = dplyr::row_number())

  Matrix::sparseMatrix(
    i = X_key$rowID,
    j = dplyr::pull(X_key, {{unit}}),
    dims = c(nrow(df), p))
}


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

      if (model_class == "multi") {
        event_size_key_time <- df_time %>%
          dplyr::group_by({{event}}) %>%
          dplyr::summarize(event_size = dplyr::n())
        X <- build_Xt(df_time,
                      event={{event}},
                      unit=paste0(unit_str, "_id"),
                      p)
        H <- build_Ht(event_size_key_time, n_t)
        X <- H %*% X
      } else {
        if (n_t == 0) { return(tibble()) }
        X <- matrix(0, nrow=n_t, ncol=p)
        for (i in 1:n_t) {
          X[i, dplyr::pull(df_time, paste0(unit_str, "1_id"))[i]] <- 1
          X[i, dplyr::pull(df_time, paste0(unit_str, "2_id"))[i]] <- -1
          # X[i, (df_time$ath2)[i]] <- -1
        }
      }

      X %>%
        as.matrix()
    }) %>%
    do.call(rbind, .)
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
        dplyr::select(cols) %>%
        unlist(use.names=F)
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
        dplyr::select(cols) %>%
        unlist(use.names=F) %>%
        unique() %>%
        length()
    })
}
