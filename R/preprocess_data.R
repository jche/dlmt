

#' Generate useful keys and optionally center data
#'
#' @param df raw data: columns unitID, eventID, t, y
#' @param outcome outcome variable
#' @param time time variable
#' @param event event variable
#' @param unit unit variable
#' @param event_center name of event-centering method: "none", "median", or "mean"
#' @param model_class class of model, either "multi" for multiplayer model or "h2h" for head-to-head model
#'
#' @return preprocessed df with helpful attributes
#'
preprocess_data <- function(df,
                            outcome = y,
                            time = t,
                            event = g,
                            unit = a,
                            event_center = c("none", "median", "mean"),
                            model_class = c("multi", "h2h")) {
  event_center <- match.arg(event_center)
  model_class <- match.arg(model_class)

  # generate unit ID key
  unit_str <- rlang::as_string(rlang::ensym(unit))
  if (model_class == "multi") {
    unique_units <- df %>%
      dplyr::pull({{unit}}) %>%
      unique() %>%
      sort()
  } else if (model_class == "h2h") {
    # record names of two unit columns
    unit1 <- paste0(unit_str, 1)
    unit2 <- paste0(unit_str, 2)
    unique_units <- c(dplyr::pull(df, unit1), dplyr::pull(df, unit2)) %>%
      unique() %>%
      sort()
  }
  unit_key <- dplyr::tibble("{{unit}}" := unique_units) %>%
    dplyr::mutate("{{unit}}_id" := 1:dplyr::n())

  # add ID to dataset
  if (model_class == "multi") {
    df <- df %>%
      dplyr::left_join(unit_key, by=unit_str) %>%
      dplyr::relocate(paste0(unit_str, "_id"),
                      .after=dplyr::all_of(unit_str))
  } else if (model_class == "h2h") {
    names(unit_key) <- c(unit1, paste0(unit1, "_id"))
    df <- df %>%
      dplyr::left_join(unit_key, by=unit1)
    names(unit_key) <- c(unit2, paste0(unit2, "_id"))
    df <- df %>%
      dplyr::left_join(unit_key, by=unit2)
    names(unit_key) <- c(unit_str, paste0(unit_str, "_id"))
  }

  # event-center dataset
  if (model_class == "h2h") {
    stopifnot("For two-player events, event centering should be \"none\"" =
                event_center == "none")
  }
  center_fun <- function(x, type) {
    switch(type,
           none = x,
           mean = x - mean(x),
           median = x - median(x))
  }
  df <- df %>%
    dplyr::group_by({{event}}) %>%
    dplyr::mutate({{outcome}} := center_fun({{outcome}}, event_center)) %>%
    dplyr::ungroup()

  # add helpful attributes to df
  attr(df, "p") <- length(unique_units)
  attr(df, "unit_key") <- unit_key
  if (model_class == "multi") {
    event_size_key <- df %>%
      dplyr::group_by({{time}}, {{event}}) %>%
      dplyr::summarize(event_size = dplyr::n(), .groups="drop")
    attr(df, "event_size_key") <- event_size_key
  }

  df
}
