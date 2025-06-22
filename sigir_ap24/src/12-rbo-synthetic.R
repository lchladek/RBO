source("src/common.R")
source("src/simulate.R")
source("rbo/low.R")
source("rbo/hig.R")
source("src/rbo_vectorized.R")

future_lapply(1:10, future.seed = 0, function(i) {
  len <- sample(10:100, 2, replace = TRUE)
  xy <- simulate_rankings(min(len), max(len),
                          n = 1000, tau = runif(1,.5,1),
                          frac_ties_x = runif(1, .1, 1), frac_ties_y = runif(1, .1, 1))
  x <- xy$x
  y <- xy$y

  n_x <- length(unlist(x))
  n_t_x <- sum(lengths(x)[lengths(x)>1])
  n_y <- length(unlist(y))
  n_t_y <- sum(lengths(y)[lengths(y)>1])


  m <- low_rankings(x,y)
  M <- hig_rankings(x,y)

  #tie-aware versions
  a <- rbo(x, y, .P, "a", c("max","min","ext"))
  #lower bound
  low <- rbo(m$x, m$y, .P, "a", c("max","min","ext"))

  #upper bound
  hig <- rbo(M$x, M$y, .P, "a", c("max","min","ext"))

  #average
  avg_max <- rbo_avg_max(x, y, .P)
  avg_min <- rbo_avg_min(x, y, .P)
  avg_ext <- rbo_avg_ext(x, y, .P)

  d <- data.frame(i, n_x, n_y, n_t_x, n_t_y,
                  p = .P,
                  rbo_a_max = a$max, rbo_a_min = a$min, rbo_a_ext = a$ext,
                  rbo_low_max = low$max, rbo_low_min = low$min, rbo_low_ext = low$ext,
                  rbo_hig_max = hig$max, rbo_hig_min = hig$min, rbo_hig_ext = hig$ext,
                  rbo_avg_max = avg_max, rbo_avg_min = avg_min, rbo_avg_ext = avg_ext)
  d

}) |> bind_rows() -> d

path <- glue("output/rbo-synthetic/")
dir.create(path, recursive = TRUE, showWarnings = FALSE)
rio::export(d, glue("{path}/rbo.csv"))
