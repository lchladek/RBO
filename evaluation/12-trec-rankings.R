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
#parser$add_argument("--time-limit", type="integer", default=1200, metavar="<seconds>",
#                    help="Time (in seconds) after which to stop calculating [default=1200]")
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
                    help="quantiles to store [default=0, 0.005, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.995, 1")
#parser$add_argument("--rbo-score", metavar="ext|res|min|max", default="min",
#                    help="Score to use for RBO (one of 'ext', 'res', 'min', and 'max') [default='min']")
#parser$add_argument("--samples", type="integer", metavar="<count>", default=0,
#                    help="Number of samples to take of the permutation space of each ranking pair. If 0, all permutations are evaluated [default=0]")
parser$add_argument("--permutation-limit", type="integer", metavar="<count>", default=100000,
                    help="The maximum number of permutations a ranking pair can have for which it is not skipped. If 0, no limit [default=100000]")
parser$add_argument("--cores", type="integer", metavar="<count>", default=0,
                    help="The number of cores for parallelized tasks [default=0/auto]")

####
#args <- parser$parse_args()
args <- list(output="output.csv", ranking_count=500, length_range=c(15, 15),
             item_count=20, seed=19, rbo_p=c(0.9), permutation_limit=10000,
             quantiles=c(0, 0.005, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.995, 1))

# Load libraries
library(rlang)
library(foreach)
library(future.apply)
plan()

# Loading the library
suppressPackageStartupMessages(source("src/lib.R"))

message(paste("Using seed:", args$seed))

library(progressr)
library(beepr)
handlers("progress", "beepr")

with_progress({
  prog <- progressor(along = seq(args$ranking_count))
  results <- future_lapply(seq(args$ranking_count), future.seed=args$seed, function(ranking_index) {
    
    if (args$length_range[1] < args$length_range[2])
      l <- sample(seq(args$length_range[1], args$length_range[2]), 1)
    else
      l <- args$length_range[1]
    n <- args$item_count
    
    # simulate the rankings
    pair <- simulate_rankings(len_x = l, len_y = l, n = n)
    
    # convert to letters
    pair <- lapply(pair, inum_ranking_to_letter)
    
    # print
    #message(paste("x: ", to_string(pair$x)))
    #message(paste("y: ", to_string(pair$y)))
    
    permutation_count <- tie_permutations(pair)
    
    pair_result <- lapply(args$rbo_p, function(p) {
    
      if ((args$permutation_limit != 0) && (args$permutation_limit < permutation_count)) {
          #message(paste("Skipping", permutation_count, "permutation pair"))
          return (NULL)
      } else {
          #message(paste("Evaluating the full distribution for", permutation_count, "tie permutations"))
          dist <- tie_distribution_all(pair, p=p, score='min')
      }
      
      low_rbo_result <- low_rankings(pair$x, pair$y)
      high_rbo_result <- hig_rankings(pair$x, pair$y)
      
      est <- estimate_pmf(pair, p=p, collapse_tol = NULL, truncation_tol = NULL)
      
      dist_quantiles <- as.data.frame(as.list(quantile(dist, probs=args$quantiles, type=1)))
      colnames(dist_quantiles) <- paste("q", args$quantiles, sep="")
      
      est_quantiles <- as.data.frame(as.list(pmf_quantile(est, probs=args$quantiles)))
      colnames(est_quantiles) <- paste("est.q", args$quantiles, sep="")
      
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
        #earth_movers=emdw(A=dist, B=est$rbo, wA=1/length(dist), wB=est$prob, distance="euclidean"),
        dist_mean=mean(dist)
      ), dist_quantiles, est_quantiles)
    })
    
    prog()
    pair_result
    
  })
  
})

out <- as.data.frame(do.call(rbind, do.call(rbind, results)))
rownames(out) <- NULL

write.csv(out, file=args$output)