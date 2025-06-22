source("src/common.R")
source("rbo/low.R")
source("rbo/hig.R")
source("src/rbo_vectorized.R")
library(future.apply)
library(dplyr)


for(year in .WEB_YEAR) {
  print(year)

  set.seed(0)

  # Read runs
  d <- read_inputs(glue("scratch/01-trec-download/web{year}/"), topk = 1000)
  meta <- rio::import(glue("data/web{year}.csv"))
  stopifnot(all(unique(meta$sys) %in% unique(d$sys)) && all(unique(d$sys) %in% unique(meta$sys)))#TODO
  topics <- sort(unique(d$topic))

  # For every pair of runs within the same group, compute RBO scores

  # First, precompute a data frame with all the comparisons' metadata
  # (will be easier to paralellize)
  r <- lapply(unique(meta$group), function(g) {
    group_syss <- meta |> filter(group == g) |> pull(sys) # runs from the same group
    if(length(group_syss) > 1) {
      lapply(utils::combn(group_syss, 2, simplify = FALSE), function(syss) {
        data.frame(sys1 = syss[1], sys2 = syss[2], topic = topics)
      }) |> bind_rows()
    }
  }) |> bind_rows()

  # And then proceed with actual calculations
  r <- future_lapply(1:nrow(r), function(i) {
    sys1 <- r$sys1[i]
    sys2 <- r$sys2[i]
    top <- r$topic[i]

    run1 <- d |> filter(sys == sys1, topic == top)
    run2 <- d |> filter(sys == sys2, topic == top)

    # Only calculate RBOs when having valid data (no missing or non-numeric scores)
    if(all(!is.na(run1$score)) && all(!is.na(run2$score)) &&
       all(!is.nan(run1$score)) && all(!is.nan(run2$score))) {

      run1ties <- extract_ranking(run1$doc, run1$score)
      run2ties <- extract_ranking(run2$doc, run2$score)

      m <- low_rankings(run1ties,run2ties)
      M <- hig_rankings(run1ties,run2ties)

      #tie-aware versions
      a <- rbo(run1ties, run2ties, .P, "a", c("max","min","ext"))

      #lower bound
      low <- rbo(m$x, m$y, .P, "a", c("max","min","ext"))

      #upper bound
      hig <- rbo(M$x, M$y, .P, "a", c("max","min","ext"))

      #average
      avg_max <- rbo_avg_max(run1ties, run2ties, .P)
      avg_min <- rbo_avg_min(run1ties, run2ties, .P)
      avg_ext <- rbo_avg_ext(run1ties, run2ties, .P)

      d <- data.frame(p = .P,
                      rbo_a_max = a$max, rbo_a_min = a$min, rbo_a_ext = a$ext,
                      rbo_low_max = low$max, rbo_low_min = low$min, rbo_low_ext = low$ext,
                      rbo_hig_max = hig$max, rbo_hig_min = hig$min, rbo_hig_ext = hig$ext,
                      rbo_avg_max = avg_max, rbo_avg_min = avg_min, rbo_avg_ext = avg_ext)

      d <- d |> mutate(sys1 = sys1, sys2 = sys2, topic = top)
      d

    }else{
      d <- data.frame(sys1, sys2, topic = top,
                 p = NA,
                 rbo_a_max = NA, rbo_a_min = NA, rbo_a_ext = NA,
                 rbo_low_max = NA, rbo_low_min = NA, rbo_low_ext = NA,
                 rbo_hig_max = NA, rbo_hig_min = NA, rbo_hig_ext = NA,
                 rbo_avg_max = NA, rbo_avg_min = NA, rbo_avg_ext = NA)
      d
    }
  }) |> bind_rows()

  # Save all
  path <- glue("output/rbo-trec")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  rio::export(r, glue("{path}/web{year}.csv"))
}
