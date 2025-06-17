#!/usr/bin/env -S Rscript --vanilla

####################################################################################################
# Copyright 2025 Lukáš Chládek <l@chla.cz>                                             MIT LICENSE #
#                                                                                                  #
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software    #
# and associated documentation files (the "Software"), to deal in the Software without             #
# restriction, including without limitation the rights to use, copy, modifS, merge, publish,       #
# distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the    #
# Software is furnished to do so, subject to the following conditions:                             #
#                                                                                                  #
# The above copyright notice and this permission notice shall be included in all copies or         #
# substantial portions of the Software.                                                            #
#                                                                                                  #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING    #
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND       #
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,     #
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          #
####################################################################################################

# Load argument parser
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--output", default="output.csv", metavar="<path.csv>",
                    help="The CSV output path")
parser$add_argument("--ranking-count", type="integer", default=10000, metavar="<seconds>",
                    help="Number of rankings after which to stop calculating [default=10000]")
parser$add_argument("--length-range", nargs=2, type="integer", default=c(4, 8), metavar="<n>",
                    help="The range in which to uniformly pick the ranking length [default=4-8]")
parser$add_argument("--item-count", type="integer", default=12, metavar="<n>",
                    help="The number of ranked items [default=12]. Must exceed the maximum ranking length.")
parser$add_argument("--seed", type="integer", metavar="<seed>", required=TRUE,
                    help="Use a seed for the RNG [required]")
parser$add_argument("--rbo-p", nargs='+', type="double", metavar="<p>", default=0.9,
                    help="p-parameter for RBO [default=0.9]")
parser$add_argument("--quantiles", nargs='+', type="double", metavar="<p>", default=c(0, 0.005, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.995, 1),
                    help="quantiles to store [default=0, 0.005, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.995, 1]")
parser$add_argument("--permutation-limit", type="integer", metavar="<count>", default=100000,
                    help="The maximum number of permutations a ranking pair can have for which it is not skipped. If 0, no limit [default=100000]")
parser$add_argument("--workers", type="integer", metavar="<count>", default=4,
                    help="The number of R workers for parallelized tasks [default=4]")

# Load libraries
library(plyr)
library(rlang)
library(future.apply)


args <- list(ranking_count=10000, length_range=c(10, 10), rbo_p=0.9, item_count=12,
             quantiles=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
             permutation_limit=100000, output="budget.csv", workers=4)

plan(multisession, workers=args$workers)

# Loading the library
suppressPackageStartupMessages(source("src/lib.R"))

results <- future_lapply(seq(args$ranking_count),  function(ranking_index) {

  # simulate the rankings
  pair <- lapply(list(x=csv$ranking_x[ranking_index], y=csv$ranking_y[ranking_index]), from_string)
  
  pair_result <- lapply(args$rbo_p, function(p) {
    
    orig_quantiles <- c(csv$q0.025[ranking_index],
                        csv$q0.05[ranking_index], 
                        csv$q0.01[ranking_index],
                        csv$q0.05[ranking_index], 
                        csv$q0.09[ranking_index],
                        csv$q0.095[ranking_index],
                        csv$q0.0975[ranking_index])
    
    aggregate_data <- data.frame()
    
    
    repeats <- 10
    
    for (truncation_threshold in seq(csv$length[ranking_index])) {
      
      step_time <- 0
      step_mses <- rep(0, 7)

      for (r in seq(repeats)) {
        a <- Sys.time()
        pmf <- estimate_pmf(pair, truncation_index = truncation_threshold)
        b <- Sys.time()
        
        quantiles <- pmf_quantile(pmf, c(0.025, 0.05, 0.1, 0.5, 0.9, 0.095, 0.0975))
      
        step_mses <- step_mses + (quantiles - orig_quantiles)^2 / repeats
        step_time <- step_time + (b - a)
      }

      msedf <- as.data.frame(as.list(step_mses))
      
      aggregate_data <- cbind(
        data.frame(as.list(step))
      )
      
    }
    
    
    pmf_data <- data.frame()
    sampling_data <- data.frame()

    for (warmup in seq(5)) {
      pmf <- estimate_pmf(pair, truncation_index = 5)
    }
    
    for (trunc in seq(15)) {
      times <- c()
      mse <- c()
      for (tr in seq(tries)) {
        a <- Sys.time()
        pmf <- estimate_pmf(pair, truncation_index = trunc)
        b <- Sys.time()
        burntime <- b - a
        
        qs <- pmf_quantile(pmf, args$quantiles)
        qdf <- as.data.frame(as.list(qs))
        names(qdf) <- args$quantiles

        times <- c(times, burntime)
        mse <- c(mse, )
        df <- cbind(
          qdf, data.frame(time=burntime, trunc=trunc)
        )
        if (is.null(df1)) {
          df1 <- df
        } else {
          df1 <- rbind(df1, df)
        }
       
      }
    }
    
    for (count in c(1, 2, 4, 8, 16, 32)) {
      for (tr in seq(tries)) {
        a <- Sys.time()
        dist <- sapply(seq(count), function(i) {
          rbo(randomize_ties(pair$x), randomize_ties(pair$y), p=0.9, score='min')
        })
        b <- Sys.time()
        burntime <- b - a
        
        qs <- quantile(dist, args$quantiles)
        qdf <- as.data.frame(as.list(qs))
        names(qdf) <- args$quantiles
        
        df <- cbind(
          qdf, data.frame(time=burntime, samp=count)
        )
        if (is.null(df2)) {
          df2 <- df
        } else {
          df2 <- rbind(df2, df)
        }
      }
    }
    
    cbind(data.frame(
      ranking_x=to_string(pair$x),
      ranking_y=to_string(pair$y),
      item_count=n,
      permutations=permutation_count,
      p=p,
      rbo_low=rbo(x=low_rbo_result$x, y=low_rbo_result$y, ties = "w", score='min', p=p),
      rbo_high=rbo(x=high_rbo_result$x, y=high_rbo_result$y, ties = "w", score='min', p=p),
      rbo_a=rbo(x=pair$x, y=pair$y, ties = "a", score='min', p=p),
      rbo_b=rbo(x=pair$x, y=pair$y, ties = "b", score='min', p=p),
      rbo_w=rbo(x=pair$x, y=pair$y, ties = "w", score='min', p=p)#,
      #earth_movers=emdw(A=dist, B=est$rbo, wA=1/length(dist), wB=est$prob, distance="euclidean"),
      #dist_mean=mean(dist)
    ), #dist_quantiles,
    sample_mses,
    est_mses)
  })
  
  message(paste(ranking_index, "complete"))
  pair_result
  
})

out <- as.data.frame(do.call(rbind, do.call(rbind, results)))
rownames(out) <- NULL

write.csv(out, file=args$output)