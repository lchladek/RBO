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

#' Finds the index of an item in a ranking
index_of <- function(ranking, item) {

  # find all ranks...
  which(laply(ranking, function(i) {

    # where at least one tied element is equal to the item
    length(which(i == item)) != 0
  }))
}


#' Checks if a ranking has ties
has_ties <- function(ranking) {
  length(unlist(ranking)) != length(ranking)
}


#' Gets the item contribution
rbo_contribution <- function(item, pair, p=0.9, ties=NULL) {
  if (is.null(ties) && sum(sapply(pair, has_ties)) != 0) {
    stop("There are ties!")
  }
  indices <- lapply(pair, index_of, item)
  if (sum(sapply(indices, is_empty)) != 0) {
    return (0)
  } else {
    m <- max(unlist(indices))
    return (
      (1 - p) / p *
        (- log1p(-p) -
          sum(sapply(seq(m - 1), function(d) {(p^d)/d}))
        )
    )
  }
}

rbo_from_contribution <- function(pair, ...) {
  sum(sapply(unlist(pair$x), rbo_contribution, pair, ...))
}

rbo_a_min <- function(pair, p=0.9, score='min', ties='a',...) {
  rbo(pair$x, pair$y, ties=ties, p=p, score=score)
}

rbo_low_min <- function(pair, p=0.9, score='min', ...) {
  low_pair <- low_rankings(pair$x, pair$y)
  rbo(low_pair$x, low_pair$y, p=p, score=score, ...)
}