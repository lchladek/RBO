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

message(paste("Loading library dependencies...", sep=""))
# libraries required here
library(combinat)
library(purrr)

# libraries required for Corsi2024b
library(dplyr)
library(extraDistr)
library(future.apply)
library(ggplot2)
library(glue)
library(latex2exp)
library(mvtnorm)
library(rio)

corsi_dir <- "sigir_ap24/"
message(paste("Loading Corsi and Urbano (2024) libraries from '", corsi_dir, "'...", sep=""))

if (!dir.exists(corsi_dir)) {
  stop("The library directory cannot be accessed.")
}

wd1 <- getwd()
setwd(corsi_dir)
  source("rbo/rbo.R")
  source("rbo/low.R")
  source("rbo/hig.R")
  source("rbo/avg.R")
  source("src/common.R")
  source("src/simulate.R")
setwd(wd1)

message("Loaded Corsi and Urbano (2024).")

#' Encode a number as letters
#' 
#' This function maps positive integers to letters as follows:
#' 1, 2, ..., 26, 27, 28, ... => a, b, ..., z, aa, ab, ...
num_to_letter <- function(num) {
  prefix = floor((num - 1) / 26)
  if (prefix == 0) {
    return (letters[num])
  } else {
    return (paste(num_to_letter(prefix), num_to_letter(((num - 1) %% 26) + 1), sep=""))
  }
}

#' Convert a ranking to letter form
#' 
#' i1 (i2 i26) i27 i28 => a (b z) aa ab
inum_ranking_to_letter <- function(ranking) {
  lapply(ranking, function(rank) {
    bad_items <- which(substr(rank, 1 ,1) != "i")
    if (length(bad_items) > 0) {
      stop(paste("Unexpected item format: ", rank[bad_items]))
    }
    nums <- as.numeric(substr(rank, 2, 1000))
    sapply(nums, num_to_letter, USE.NAMES = FALSE)
  })
}

#' Get the number of tie permutations per truncation depth
#' 
#' This function takes a ranking pair as a list and returns a list of vectors,
#' where each vector is the number of permutations of ties for that ranking up
#' to depth i.
tie_permutations_per_depth <- function(ranking_pair) {
  lapply(ranking_pair, function(ranking) {
    cumprod(unlist(lapply(ranking, function(rank) {seq(length(rank))})))
  })
}

#' Get the number of tie permutations
#' 
#' This function takes a list of rankings and counts the total number of
#' permutations produced by arranging all ties in the rankings.
tie_permutations <- function(ranking_pair) {
  prod(sapply(tie_permutations_per_depth(ranking_pair), last))
}

#' Break ties in a ranking at random
#' 
#' This function takes a ranking (use `lapply(pair, randomize_ties)` for a 
#' ranking pair) and replaces all ties with a random permutation of single ranks.
randomize_ties <- function(ranking) {
  list_flatten(sapply(ranking, function(rank) {
    as.list(sample(rank, length(rank), replace=FALSE))
  }))
}

#' Get all permutations of a tied ranking
#' 
#' This function takes a ranking (not a pair) and returns a matrix of all
#' tie permutations, such that row i is all possibilities for rank i and column
#' i is the ith possible ranking with broken ties.
get_ranking_tie_combinations <- function(ranking) {
  combos <- expand.grid(
    lapply(ranking, permn)
  )
  apply(combos, MARGIN=1, FUN=unlist)
}

#' Gets the RBO for all permutations of a tied ranking pair
#' 
#' This function takes a ranking pair and returns a vector of all RBO values
#' stemming from every possible permutation of ties in the two rankings.
tie_distribution_all <- function(ranking_pair, p=0.9, score='ext') {
  combinations <- lapply(ranking_pair, get_ranking_tie_combinations)
  as.vector(
    apply(combinations$x, MARGIN=2, FUN=function(xcomb) {
      apply(combinations$y, MARGIN=2, FUN=function(ycomb) {
        rbo(xcomb, ycomb, ties="w", p, score)
      })
    }))
}

#' Gets the RBO for all permutations of a tied ranking pair as a PMF
#' 
#' This function takes a ranking pair and returns a PMF data frame of all tie
#' permutations.
complete_pmf <- function(pair, p=0.9, score="min") {
  lx <- length(pair$x)
  ly <- length(pair$y)
  permutation <- rep(1, lx + ly)
  max_permutation <- factorial(c(sapply(pair$x, length), sapply(pair$y, length)))
  incrementable_indices <- seq_along(permutation)[which(max_permutation != 1)]
  distribution <- data.frame()
  while (TRUE) {
    # next permutation
    x <- unlist(lapply(seq(lx), function(i) {
      get_indexed_permutation(permutation[i], pair$x[[i]])
    }))
    y <- unlist(lapply(seq(ly), function(i) {
      get_indexed_permutation(permutation[lx + i], pair$y[[i]])
    }))
    
    rbo_value <- rbo(x, y, p, ties="w", score)
    names(rbo_value) <- NULL
    rbo_index <- which(distribution$rbo == rbo_value)
    if (length(rbo_index) != 0) {
      distribution$count[rbo_index] <- distribution$count[rbo_index] + 1
    } else {
      distribution <- rbind(distribution, data.frame(rbo=rbo_value, count=1))
    }
    
    incremented <- FALSE
    for (bit in incrementable_indices) {
      if (permutation[bit] < max_permutation[bit]) {
        permutation[bit] <- permutation[bit] + 1
        incremented <- TRUE
        break
      }
      permutation[bit] <- 1
    }
    if (!incremented) {
      sorted <- distribution[order(distribution$rbo), ]
      rownames(sorted) <- NULL
      return(data.frame(rbo=sorted$rbo, prob=sorted$count/prod(max_permutation)))
    }
  }
}

get_indexed_permutation <- function(index, options) {
  if (length(options) >= 14) {
    stop("Permutations of length 14 or higher are not supported.")
  }
  if (length(options) <= 1) {
    return (options)
  }
  selected_index <- ((index - 1)  %% length(options)) + 1
  rest_index <- as.integer(1 + (index - 1) / length(options))
  return (c(options[selected_index],
           get_indexed_permutation(rest_index, options[-selected_index])))
}

bits_to_str <- function(bits) {
  paste(bits, collapse = '')
}

str_to_bits <- function(str) {
  sapply(strsplit(str, "")[[1]], as.integer, USE.NAMES=FALSE)
}

basic_convolution <- function(probabilities, ms, probs) {
  new_probabilities <- list()
  for (i in seq_along(probabilities)) {
    orig_bits <- str_to_bits(names(probabilities)[i])
    orig_prob <- probabilities[[i]]
    num_bits <- length(orig_bits)
    for (j in seq_along(ms)) {
      m <- ms[j]
      p <- probs[j]
      
      new_bits <- orig_bits
      new_bits[m] <- orig_bits[m] + 1
      
      new_bits_str <- bits_to_str(new_bits)
      curr_p <- new_probabilities[[ new_bits_str ]]
      if (!is.null(curr_p))
        new_probabilities[[ new_bits_str ]] <- orig_prob * p + curr_p
      else
        new_probabilities[[ new_bits_str ]] <- orig_prob * p
    }
  }
  new_probabilities
}

cull_twos_convolution <- function(probabilities, ms, probs) {
  new_probabilities <- list()
  for (i in seq_along(probabilities)) {
    orig_bits <- str_to_bits(names(probabilities)[i])
    orig_prob <- probabilities[[i]]
    num_bits <- length(orig_bits)
    for (j in seq_along(ms)) {
      m <- ms[j]
      p <- probs[j]
      
      if (orig_bits[m] >= 2)
        next
      
      new_bits <- orig_bits
      new_bits[m] <- orig_bits[m] + 1
      
      new_bits_str <- bits_to_str(new_bits)
      curr_p <- new_probabilities[[ new_bits_str ]]
      if (!is.null(curr_p))
        new_probabilities[[ new_bits_str ]] <- orig_prob * p + curr_p
      else
        new_probabilities[[ new_bits_str ]] <- orig_prob * p
    }
  }
  
  # normalize
  psum <- sum(unlist(new_probabilities))
  lapply(new_probabilities, function(p) {
    p / psum
  })
}

cull_all_convolution <- function(probabilities, ms, probs) {
  new_probabilities <- list()
  for (i in seq_along(probabilities)) {
    orig_bits <- str_to_bits(names(probabilities)[i])
    orig_prob <- probabilities[[i]]
    num_bits <- length(orig_bits)
    for (j in seq_along(ms)) {
      m <- ms[j]
      p <- probs[j]
      
      new_bits <- orig_bits
      new_bits[m] <- orig_bits[m] + 1
      
      exceeds_sequence <- cumsum(new_bits) > seq(num_bits)
      exceeds_two <- new_bits > 2
      
      if (sum( exceeds_two | exceeds_sequence ) > 0)
        next
      
      new_bits_str <- bits_to_str(new_bits)
      curr_p <- new_probabilities[[ new_bits_str ]]
      if (!is.null(curr_p))
        new_probabilities[[ new_bits_str ]] <- orig_prob * p + curr_p
      else
        new_probabilities[[ new_bits_str ]] <- orig_prob * p
    }
  }
  
  if (sum(unlist(new_probabilities)) == 0)
    return (probabilities)
  
  # normalize
  psum <- sum(unlist(new_probabilities))
  lapply(new_probabilities, function(p) {
    p / psum
  })
}

contribution_constants <- function(max_possible_m, p=0.9, truncation_tol=NULL) {
  max_m <- 1
  contribution_curr <- (1 - p) / p * (- log1p(-p) )
  contribution_cumulative <- contribution_curr
  contribution_by_m <- c(contribution_curr)
  
  while (max_m < max_possible_m) {
    max_m <- max_m + 1
    
    d <- max_m - 1
    contribution_curr <- contribution_curr - (1 - p) / p * ((p^d) / d)
    contribution_by_m <- c(contribution_by_m, contribution_curr)
    contribution_cumulative <- contribution_cumulative + contribution_curr
    if (!is.null(truncation_tol)) {
      if (1 - contribution_cumulative < truncation_tol) break;
    }
  }
  
  contribution_by_m
}

estimate_pmf <- function(pair, p=0.9, convolution_fn=cull_all_convolution, truncation_tol=NULL, truncation_index=NULL) {
  
  flat <- lapply(pair, unlist)
  len <- lapply(flat, length)
  max_possible_m <- max(len$x, len$y)
  if (!is.null(truncation_index)) {
    max_possible_m <- min(max_possible_m, truncation_index)
  }
  tie_lengths <- lapply(pair, sapply, length)
  tie_last_indices <- lapply(tie_lengths, cumsum)
  tie_first_indices <- list(x = tie_last_indices$x - tie_lengths$x + 1L,
                            y = tie_last_indices$y - tie_lengths$y + 1L)
    
  max_m <- 1
  contribution_curr <- (1 - p) / p * (- log1p(-p) )
  contribution_cumulative <- contribution_curr
  contribution_by_m <- c(contribution_curr)
  
  while (max_m < max_possible_m) {
    max_m <- max_m + 1
    
    d <- max_m - 1
    contribution_curr <- contribution_curr - (1 - p) / p * ((p^d) / d)
    contribution_by_m <- c(contribution_by_m, contribution_curr)
    contribution_cumulative <- contribution_cumulative + contribution_curr
    if (!is.null(truncation_tol)) {
      if (1 - contribution_cumulative < truncation_tol) break;
    }
  }
  
  #l <- min(lx, ly)
  #items <- c(items_in_x, items_in_y)[which(c(lx, ly) == min(lx, ly))]
  
  # last element that can appear below m_max in this ranking
  top_x_index <- tie_last_indices$x[which(tie_last_indices$x >= max_m)[1]]
  if (is.na(top_x_index))
    top_x_index <- len$x
  
  probabilities <- list()
  probabilities[[bits_to_str( bitstring(max_m) )]] <- 1
  
  for (x_index in seq(top_x_index)) {
    
    # starting with item i of x
    i <- flat$x[x_index]
    
    # get i's position in flat y
    y_index <- which(flat$y == i)[1]
    if (is.na(y_index))
      # i is not in y
      next;
    
    
    
    y_rank <- which(tie_last_indices$y >= y_index)[1]
    y_first <- tie_first_indices$y[y_rank]
    y_last <- tie_last_indices$y[y_rank]
    
    x_rank <- which(tie_last_indices$x >= x_index)[1]
    x_first <- tie_first_indices$x[x_rank]
    x_last <- tie_last_indices$x[x_rank]
    
    # lower bound on M_i
    first <- max(x_first, y_first)
    if (first > max_m)
      next;
    
    # upper bound on M_i
    last <- max(x_last, y_last)
    # if (first == last) {
    #   # M_i can only take one value
    #   rbos <- rbos + contribution_by_m[first]
    #   next;
    # }
    
    possible_ms <- seq(first, min(max_m, last))
    possible_ms_probs <- sapply(possible_ms, function(m) {
      (
        min(x_last - x_first, max(0, m - x_first)) * (y_first <= m && y_last >= m)
        + min(y_last - y_first, max(0, m - y_first)) * (x_first <= m && x_last >= m)
        + 1
      ) / 
        ((x_last - x_first + 1) * (y_last - y_first + 1))
    })
    
    probabilities <- convolution_fn( probabilities, possible_ms, possible_ms_probs )
    
  }
  
  rbos <- sapply(names(probabilities), function(str) {
    bits <- str_to_bits(str)
    sum(bits * contribution_by_m)
  }, USE.NAMES=FALSE)
  probs <- unlist(probabilities)
  
  df <- data.frame(rbo=rbos, prob=probs)
  df <- df[order(rbos), ]
  rownames(df) <- seq_along(rbos)
  df
  
}

pmf_quantile <- function(pmf, probs) {
  cum_p <- cumsum(pmf$prob)
  
  sapply(probs, function(prob) {
    index <- head(n=1, which(cum_p > prob))
    if (length(index) == 0) {
      pmf$rbo[length(pmf$rbo)]
    } else {
      pmf$rbo[index]
    }
  })
  
}

pmf_mean <- function(pmf) {
  weighted.mean(pmf$rbo, pmf$prob)
}

pmf_variance <- function(pmf) {
  mean <-  weighted.mean(pmf$rbo, pmf$prob)
  squared_deviation <- (pmf$rbo - mean)^2
  weighted.mean(squared_deviation, pmf$prob)
}

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


message("")
message("lib.R by lchladek is ready!")
