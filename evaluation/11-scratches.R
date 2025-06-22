# run 11-synthetic-rankings first


qdata <- lapply(args$quantiles, function(q) {
  dist <- with(csv, get(paste("q", q, sep="")))
  est <- with(out, get(paste("est.q", q, sep="")))
  est - dist
})
names(qdata) <- args$quantiles

boxplot(qdata)

low_rbo_result <- low_rankings(pair$x, pair$y)
high_rbo_result <- hig_rankings(pair$x, pair$y)

rbo_low <- rbo(x=low_rbo_result$x, y=low_rbo_result$y, ties = "w", score='min', p=p)
rbo_high <- rbo(x=high_rbo_result$x, y=high_rbo_result$y, ties = "w", score='min', p=p)

qdata <- lapply(args$quantiles, function(q) {
  dist <- with(out, get(paste("q", q, sep="")))
  est <- with(out, get(paste("est.q", q, sep="")))
  (est - dist)/dist
})
names(qdata) <- args$quantiles
boxplot(qdata)

boxplot(list(a=out$q0 - out$q0.025, b=out$est.q0 - out$q0))
boxplot(out$q0 - out$q0.025)

sim_dist <- unlist(sapply(seq(nrow(est)), function(i) {rep(est$rbo[i], 10000*est$prob[i])}))


hist(dist)
hist(sim_dist)

hist(sim_dist2)

#=====

pair <- simulate_rankings(len_x=15, len_y=15, n=25)
print(lapply(pair, to_string))

all <- tie_distribution_all(pair, score="min")
pmf <- estimate_pmf(pair)
pmf_cull_twos <- estimate_pmf(pair, convolution_fn = cull_twos_convolution)
pmf_cull_all <- estimate_pmf(pair, convolution_fn = cull_all_convolution)
pmf_cull_all2 <- estimate_pmf(pair, convolution_fn = cull_all_convolution2)

#=====

pair <- lapply(list(
  x="a b c d e f g h i",
  y="(w x y z g f p q a b)"
), from_string)


all <- tie_distribution_all(pair, score="min")
hist(all)

low_rbo_result <- low_rankings(pair$x, pair$y)
high_rbo_result <- hig_rankings(pair$x, pair$y)

eval_budget <- function() {
  {
    n <- 1000
    set.seed(101)
    rankings <- lapply(seq(n), simulate_rankings, len_x=15, len_y=15, n=30)
  }
  
  {
    for (i in seq(n)) {
      pair <- rankings[[i]]
      tries <- 10
      df <- data.frame()
  
      for (trunc in seq(15)) {
        for (tr in seq(tries)) {
          a <- Sys.time()
          pmf <- estimate_pmf(pair, truncation_index = trunc)
          b <- Sys.time()
          qs <- pmf_quantile(pmf, args$quantiles)
          qdf <- as.data.frame(as.list(qs))
          names(qdf) <- args$quantiles
          
          df <- rbind(df, cbind(
            qdf, data.frame(time=b-a)
          ))
        }
      }
  
    }
  }
}


rbo_low <- rbo(x=low_rbo_result$x, y=low_rbo_result$y, ties = "w", score='min', p=p)
rbo_high <- rbo(x=high_rbo_result$x, y=high_rbo_result$y, ties = "w", score='min', p=p)


#====

bitstring_to_int <- function(s) {
  sum(3 ^ (length(s) - seq_along(s)) * s)
}

#===
# speed comparison of the full distribution evaluation

set.seed(50)
pairs <- list()
while(length(pairs) < 10) {
  pair <- try(simulate_rankings(10, 10, 20))
  if (tie_permutations(pair) < 100000) {
    pairs[[length(pairs) + 1]] <- pair
  }
}

times <- sapply(
  seq_along(pairs),
  function(i) { 
    message(i)
    pair <- pairs[[i]]
    o <- Sys.time()
    tie_distribution_all(pair, p=0.9, score="min")
    m <- Sys.time() 
    tie_distribution_all2(pair, p=0.9, score="min")
    n <- Sys.time()
    c(m-o, n-m)
  }
)
rownames(times) <- c("old", "new")

View(times)

print(paste(time_a, "should be a faster permutation algo than", time_b))

#====

bits_to_str <- function(bits) {
  paste(bits, collapse = '')
}

str_to_bits <- function(str) {
  sapply(strsplit(str, "")[[1]], as.integer, USE.NAMES=FALSE)
}

bitstring <- function(len, ones=c()) {
  bits <- rep(0L, len)
  bits[ones] <- 1L
  bits
}
