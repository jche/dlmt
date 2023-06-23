
#####
# map "multicompetitor" df with 2 players per event into a "h2h" df
#####

to_h2h <- function(df,
                   outcome = y,
                   time = t,
                   game = g,
                   athlete = a) {
  # check that data is indeed head-to-head
  temp <- df %>%
    dplyr::group_by({{game}}) %>%
    dplyr::summarize(n = dplyr::n())
  stopifnot("Data do not have exactly 2 competitors in some events" = max(temp$n) == 2)

  df %>%
    dplyr::mutate(rowID = 1:dplyr::n()) %>%   # so that draws can be preserved in data
    dplyr::group_by({{time}}, {{game}}) %>%
    dplyr::mutate(athlete_number = c("athlete1", "athlete2")) %>%
    tidyr::pivot_wider(
      names_from = athlete_number,
      # names_prefix = "{{athlete}}",
      values_from = {{athlete}}) %>%
    dplyr::summarize(
      {{outcome}} := {{outcome}}[1] - {{outcome}}[2],
      "{{athlete}}1" := athlete1[1],
      "{{athlete}}2" := athlete2[2]) %>%
    dplyr::ungroup()
}

