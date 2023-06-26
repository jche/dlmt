
#####
# map "multicompetitor" df with 2 players per event into a "h2h" df
#####

to_h2h <- function(df,
                   outcome = y,
                   time = t,
                   event = g,
                   unit = a) {
  # check that data is indeed head-to-head
  temp <- df %>%
    dplyr::group_by({{event}}) %>%
    dplyr::summarize(n = dplyr::n())
  stopifnot("Data do not have exactly 2 competitors in some events" = max(temp$n) == 2)

  df %>%
    dplyr::mutate(rowID = 1:dplyr::n()) %>%   # so that draws can be preserved in data
    dplyr::group_by({{time}}, {{event}}) %>%
    dplyr::mutate(unit_number = c("unit1", "unit2")) %>%
    tidyr::pivot_wider(
      names_from = unit_number,
      # names_prefix = "{{unit}}",
      values_from = {{unit}}) %>%
    dplyr::summarize(
      {{outcome}} := {{outcome}}[1] - {{outcome}}[2],
      "{{unit}}1" := unit1[1],
      "{{unit}}2" := unit2[2]) %>%
    dplyr::ungroup()
}

