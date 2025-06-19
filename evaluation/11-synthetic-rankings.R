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
parser$add_argument("--quantiles", nargs='+', type="double", metavar="<p>", default=c((0:100)/100, 0.025, 0.975),
                    help="quantiles to compute [default=0-1 in 0.01 increments plus 0.025, 0.975]")
parser$add_argument("--earth-movers-distance", action="store_true",
                    help="Compute the Earth Mover's distance")
parser$add_argument("--permutation-limit", type="integer", metavar="<count>", default=100000,
                    help="The maximum number of permutations a ranking pair can have for which it is not skipped. If 0, no limit [default=100000]")
parser$add_argument("--workers", type="integer", metavar="<count>", default=1,
                    help="The number of R workers for parallelized tasks [default=1 (sequential)]")

# Load libraries
library(rlang)
library(emdist)
library(readr)

args <- parser$parse_args()

setwd("~/RBO")

# Loading the library
suppressPackageStartupMessages(source("src/lib.R"))

if (args$workers == 1) {
  message("Using sequential runner.")
  plan(sequential)
} else {
  message("Using parallel runner.")
  plan(multisession, workers=args$workers)
}

message(paste("Using seed:", args$seed))
set.seed(args$seed)

if (args$length_range[1] < args$length_range[2]) {
  ls <- sample(seq(args$length_range[1], args$length_range[2]), args$ranking_count, replace=TRUE)
} else {
  ls <- rep(args$length_range[1], args$ranking_count) 
}


message("Generating rankings.")

pairs <- list()
while(length(pairs) < args$ranking_count) {
  curr_idx <- length(pairs) + 1
  pair <-  try(simulate_rankings(len_x=ls[curr_idx], len_y=ls[curr_idx], n=args$item_count))
  if(inherits(pair, "try-error")) {
    message("Check your item configuration..")
    stop()
  }
  if (tie_permutations(pair) < args$permutation_limit) {
    message(".", appendLF=FALSE)
    pairs[[curr_idx]] <-  lapply(pair, inum_ranking_to_letter)
  } else {
    message("x", appendLF=FALSE)
  }
}

message("Rankings generated.")

is_first <<- TRUE
output_df <- function(df) {
  if (is_first) {
    is_first <<- FALSE
    cat(format_csv(df, col_names=TRUE))
  } else {
    cat(format_csv(df, col_names=FALSE))
  }
}

message(paste("Generated", args$ranking_count, "rankings"))
message("Now computing permutations...")
  
results <- lapply(seq(args$ranking_count), function(ranking_index) {
  
  l <- ls[ranking_index]
  pair <- pairs[[ranking_index]]
  n <- args$item_count
  
  permutation_count <- tie_permutations(pair)
  
  message(paste(ranking_index, ":", permutation_count, "permutations"))
  
  pair_result <- lapply(args$rbo_p, function(p) {
  
    message(paste(ranking_index, "evaluating at p", p))
    
    low_rbo_result <- low_rankings(pair$x, pair$y)
    high_rbo_result <- hig_rankings(pair$x, pair$y)
    
  
    pmf <- complete_pmf(pair, p=p, score='min')
    est <- estimate_pmf(pair, p=p, truncation_tol = NULL, convolution_fn = cull_all_convolution)
    
    pmf_quantiles <- as.data.frame(as.list(pmf_quantile(pmf, probs=args$quantiles)))
    colnames(pmf_quantiles) <- paste("pmf.q", args$quantiles, sep="")
    
    est_quantiles <- as.data.frame(as.list(pmf_quantile(est, probs=args$quantiles)))
    colnames(est_quantiles) <- paste("est.q", args$quantiles, sep="")
    
    if (args$earth_movers_distance) {
      emd <- emdw(A=pmf$rbo, B=est$rbo, wA=pmf$prob, wB=est$prob, distance="euclidean")
    } else {
      emd <- NA
    }
    
    # try to compute KL divergence
    indices <- sapply(pmf$rbo, function(val) {
      which(est$rbo == val)[1]
    })
    if (anyNA(indices)) {
      # nonconjoint: can't compute KL
      kld <- NA
    } else {
      kld <- sum(
        pmf$prob * log( pmf$prob / est$prob[indices])
      ) 
    }
    
    cbind(data.frame(
      ranking_x=to_string(pair$x),
      ranking_y=to_string(pair$y),
      length=l,
      item_count=n,
      permutations=permutation_count,
      p=p,
      rbo_low=rbo(x=low_rbo_result$x, y=low_rbo_result$y, ties = "w", score='min', p=p),
      rbo_high=rbo(x=high_rbo_result$x, y=high_rbo_result$y, ties = "w", score='min', p=p),
      rbo_a=rbo(x=pair$x, y=pair$y, ties = "a", score='min', p=p),
      rbo_b=rbo(x=pair$x, y=pair$y, ties = "b", score='min', p=p),
      rbo_w=rbo(x=pair$x, y=pair$y, ties = "w", score='min', p=p),
      kl_divergence=kld,
      earth_movers_distance=emd,
      pmf.mean=pmf_mean(pmf),
      est.mean=pmf_mean(est),
      pmf.var=pmf_variance(pmf),
      est.var=pmf_variance(est),
      pmf.count=length(pmf$rbo),
      est.count=length(est$rbo)
    ), 
    pmf_quantiles,
    est_quantiles)
  })

  message(paste(ranking_index, "complete"))

  lapply(pair_result, output_df)
  
  pair_result
  
})

message("All done.")