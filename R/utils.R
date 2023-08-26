
# utility functions

# compute treatment effects
get_att_ests <- function(matched_df) {
  matched_df %>%
    group_by(Z) %>%
    summarize(mn = sum(Y*weights) / sum(weights)) %>%
    summarize(est = last(mn) - first(mn)) %>%
    pull(est)
}


#' aggregate scweights to tx/sc units
#'
#' @param scweights list of sc weights (tibbles with tx unit in first row)
#'
#' @return df with tx and corresponding sc units
agg_sc_units <- function(scweights) {
  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      map_dfr(~mutate(., subclass=id[1]))
  }

  scweights %>%
    group_by(subclass, Z) %>%
    summarize(across(starts_with("X"),
                     ~sum(.x * weights)),
              across(starts_with("Y"),
                     ~sum(.x * weights)),
              .groups="drop_last") %>%
    mutate(id = c(NA, subclass[1]), .before="subclass") %>%
    mutate(weights = 1) %>%
    ungroup()
}


#' aggregate sc weights to original tx and (weighted) co units
#'
#' @param scweights list of sc weights (tibbles with tx unit in first row)
#'
#' @return df with original tx and sc units, with weights
agg_co_units <- function(scweights) {
  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      bind_rows()
  }

  scweights %>%
    group_by(id) %>%
    summarize(across(-contains("weights"), ~first(.)),
              weights = sum(weights),
              subclass = NA,
              dist = NA)
}


agg_avg_units <- function(scweights) {
  if (!is.data.frame(scweights)) {
    scweights <- scweights %>%
      map_dfr(~mutate(., subclass=id[1]))
  }

  scweights %>%
    group_by(subclass, Z) %>%
    summarize(across(starts_with("X"),
                     ~sum(.x * 1/n())),
              across(starts_with("Y"),
                     ~sum(.x * 1/n()))) %>%
    mutate(id = c(NA, subclass[1]), .before="subclass") %>%
    mutate(weights = 1) %>%
    ungroup()
}

