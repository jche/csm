
# functions for conducting matching

# function 1:
#  1. find caliper matches
#  2. compute sc weights
#  3. aggregate appropriately (default tx, weighted co)

# function 2: produce tx effect estimate!

# idea: keep list of tx/co pairings as default
csm <- function(df,
                formula = NULL,
                metric = c("maximum", "euclidean", "manhattan"),
                caliper = 1,
                cal_method = c("adaptive", "fixed", "1nn"),
                est_method = c("scm", "average"),
                return = c("sc_units", "agg_co_units", "all"),
                dist_scaling = NULL,
                knn = NULL) {
  metric <- match.arg(metric)
  cal_method <- match.arg(cal_method)
  est_method <- match.arg(est_method)
  return <- match.arg(return)

  # rename columns:
  #  - treatment: Z
  #  - covs: X1 through Xp
  if (is.null(formula)) {
    df <- df %>%
      dplyr::select(Z, dplyr::starts_with("X"))
  } else {
    df <- df %>%
      dplyr::select(formula[[2]], all.vars(formula[[3]]))
  }
  df <- df %>%
    dplyr::mutate(id = 1:dplyr::n(), .before=Z)
  p <- ncol(df)-2
  names(df) <- c("id", "Z", paste0("X", 1:p))

  # if not given, build scaling vector for distance metric
  #  - used for caliper matching AND scm
  # NOTE: NA values indicate exact matches
  if (is.null(dist_scaling)) {
    dist_scaling <- df %>%
      dplyr::summarize(dplyr::across(dplyr::starts_with("X"),
                                     function(x) {
                                       if (is.numeric(x)) 1/sd(x)
                                       else NA
                                     })) %>%
      as.numeric()
  } else if (length(dist_scaling) == 1) {
    # if scaling is scalar c, scale all covariates by c
    dist_scaling <- rep(dist_scaling, p)
  }

  ### use cal_method: generate matches

  # get caliper matches
  cal_matches <- df %>%
    get_cal_matches(
      dist_scaling = dist_scaling,
      metric = metric,
      caliper = caliper,
      cal_method = cal_method,
      knn = knn)

  ### use est_method: scm or average

  weighted_matches <- cal_matches$matches %>%
    gen_weights(dist_scaling = dist_scaling,
                est_method = est_method,
                metric = metric)

  ### use return: sc units, aggregated weights per control, all weights
  res <- switch(return,
                sc_units     = agg_sc_units(weighted_matches),
                agg_co_units = agg_co_units(weighted_matches),
                all          = dplyr::bind_rows(weighted_matches))

  unmatched_units <- setdiff(df %>% filter(Z==1) %>% pull(id),
                             res %>% filter(Z==1) %>% pull(id))
  if (length(unmatched_units) > 0) {
    warning(glue::glue("Dropped the following treated units from data:
                        \t {paste(unmatched_units, collapse=\", \")}"))
  }

  browser()

  ### store some useful attributes

  # store information about scaling and adaptive calipers
  adacalipers_df <- dplyr::tibble(
    id = dplyr::filter(df, Z == 1)$id,
    adacal = cal_matches$adacalipers)
  attr(res, "scaling") <- dist_scaling
  attr(res, "adacalipers") <- adacalipers_df

  # store information about feasible units/subclasses
  feasible_units <- dplyr::filter(adacalipers_df, adacal <= caliper)$id
  attr(res, "unmatched_units") <- unmatched_units
  attr(res, "feasible_units")  <- feasible_units
  attr(res, "feasible_subclasses") <- dplyr::filter(res, id %in% feasible_units)$subclass

  # keep a bunch of data around, just in case
  attr(res, "weighted_matches") <- weighted_matches
  attr(res, "dm") <- cal_matches$dm
  attr(res, "dm_uncapped") <- cal_matches$dm_uncapped

  return(res)
}


# matching method ---------------------------------------------------------

# get caliper matches
get_cal_matches <- function(df,
                            dist_scaling,
                            metric,
                            caliper,
                            cal_method = c("adaptive", "1nn", "fixed"),
                            knn = NULL) {
  cal_method <- match.arg(cal_method)

  # store helpful constants
  ntx <- sum(df$Z)     # number of tx units
  p   <- ncol(df)-2    # number of matched covariates

  ### step 0: generate distance matrix

  dm <- gen_dm(df,
               dist_scaling=dist_scaling,
               metric=metric)
  dm_uncapped <- dm   # store uncapped distance matrix

  ### step 1: get caliper for each treated unit

  if (cal_method == "adaptive") {
    # default: p+1 nearest neighbors
    if (is.null(knn)) {
      knn <- p+1
    }

    min_dists <- 1:ntx
    for (i in 1:ntx) {
      temp <- dm[i,]
      temp_sorted <- sort(temp)

      # idea:
      #  - if 1nn farther than caliper: stretch to 1nn
      #  - if 1nn closer than caliper: shrink to knn, if possible
      if (temp_sorted[1] > caliper) {
        min_dists[i] <- temp_sorted[1]
      } else {
        min_dists[i] <- min(caliper, temp_sorted[knn])
      }
    }
  } else if (cal_method == "1nn") {
    min_dists <- apply(dm, 1, min)
  } else {
    min_dists <- rep(caliper, nrow(dm))
  }

  ### step 2: remove all control units farther than caliper

  for (i in 1:ntx) {
    temp <- dm[i,]
    temp[temp > min_dists[i]] <- NA
    dm[i,] <- temp

    # if min_dists is too big (i.e., no exact match on discrete covs),
    #  ensure that no matches are returned
    if (min_dists[i] >= 1000) {
      dm[i,] <- NA
    }
  }

  ### step 3: generate df of matched controls for each treated unit

  df_list <- purrr::map(1:ntx, function(x) {
    # record matches
    vec <- dm[x,]
    vec <- vec[!is.na(vec)]

    # if no matches, drop treated unit
    if (length(vec) == 0) {
      return(NULL)
    }

    # record id of tx/matched co units
    tx_id <- rownames(dm)[x]
    co_ids <- names(vec)

    df %>%
      dplyr::filter(id %in% c(tx_id, co_ids)) %>%
      dplyr::mutate(# dist = c(0, vec),
        subclass = x, .after="id")
  })

  # store number of co matches per tx unit
  num_matches <- rowSums(!is.na(dm))

  return(list(matches = df_list %>% purrr::discard(is.null),   # drop unmatched tx units
              adacalipers = min_dists,
              dm = dm,
              dm_uncapped = dm_uncapped,
              num_matches = num_matches
  ))
}




# estimate within matched sets --------------------------------------------

# Generate weights for list of matched sets
gen_weights <- function(matched_gps,
                        dist_scaling,
                        est_method = c("scm", "average"),
                        metric = c("maximum", "euclidean", "manhattan")) {
  est_method = match.arg(est_method)
  metric = match.arg(metric)

  if (est_method == "scm") {
    gp_weights <- purrr::map(matched_gps,
                             ~gen_sc_weights(.x, dist_scaling, metric),
                             .progress="Producing SCM units...")   # add progress bar
  } else if (est_method == "average") {
    gp_weights <- purrr::map(matched_gps,
                             function(x) {
                               x %>%
                                 dplyr::group_by(Z) %>%
                                 dplyr::mutate(weights = 1/dplyr::n()) %>%
                                 dplyr::ungroup()
                             })
  }
  return(gp_weights)
}
