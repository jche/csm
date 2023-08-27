
# utility functions

# compute treatment effects
get_att_ests <- function(matched_df) {
  matched_df %>%
    group_by(Z) %>%
    summarize(mn = sum(Y*weights) / sum(weights)) %>%
    summarize(est = last(mn) - first(mn)) %>%
    pull(est)
}


#' aggregate to tx/sc units
#'
#' @param wl list of tibbles with tx unit in first row, and associated weights
#'
#' @return df with tx and corresponding sc units
agg_sc_units <- function(wl) {
  if (!is.data.frame(wl)) {
    wl <- wl %>%
      purrr::map_dfr(~mutate(., subclass=id[1]))
  }

  wl %>%
    dplyr::group_by(subclass, Z) %>%
    dplyr::summarize(dplyr::across(dplyr::starts_with("X"),
                                   ~sum(.x * weights)),
                     .groups="drop_last") %>%
    dplyr::mutate(id = c(NA, subclass[1]), .before="subclass") %>%
    dplyr::mutate(weights = 1) %>%
    dplyr::ungroup()
}


#' aggregate to original tx and (weighted) co units
#'
#' @param wl list of tibbles with tx unit in first row, and associated weights
#'
#' @return df with original tx and sc units, with weights
agg_co_units <- function(wl) {
  if (!is.data.frame(wl)) {
    wl <- dplyr::bind_rows(wl)
  }

  wl %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(
      dplyr::across(-dplyr::contains("weights"), ~dplyr::first(.)),
      weights = sum(weights))
}

