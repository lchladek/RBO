####################################################################################################
# Copyright 2024 Matteo Corsi <corsi.matteo95@gmail.com>                               MIT LICENSE #
#                Julián Urbano <urbano.julian@gmail.com>                                           #
#                                                                                                  #
# Permission is hereby granted, free of charge, to anS person obtaining a copy of this software    #
# and associated documentation files (the "Software"), to deal in the Software without             #
# restriction, including without limitation the rights to use, copy, modifS, merge, publish,       #
# distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the    #
# Software is furnished to do so, subject to the following conditions:                             #
#                                                                                                  #
# The above copyright notice and this permission notice shall be included in all copies or         #
# substantial portions of the Software.                                                            #
#                                                                                                  #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF AnS KIND, EXPRESS OR IMPLIED, INCLUDING    #
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND       #
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR AnS CLAIM,     #
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          #
####################################################################################################
# Always check for the latest version at https://github.com/matteo-corsi/sigir_ap24                #
#                                                                                                  #
# If you use this software, please cite the following paper:                                       #
#                                                                                                  #
# M. Corsi and J. Urbano, "How do Ties Affect the Uncertainty in Rank-Biased Overlap?",            #
# International ACM SIGIR Conference on Information Retrieval in the Asia Pacific, 2024            #
####################################################################################################

library(dplyr)
source("rbo/rbo.R")

low_rankings <- function(x,y) {

  # If no ties, do nothing
  if(all(lengths(x) == 1) && all(lengths(y) == 1))
    return(list(x = x, y = y))

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L

  fL <- flatten(x) |> bind_cols()
  fS <- flatten(y) |> bind_cols()
  f <- full_join(fL, fS, by=c("id"="id")) |>
    mutate(t.o = pmax(t.x,t.y),
           b.o = pmax(b.x,b.y),
           size = b.o-t.o,
           distance = pmax(b.x-t.y, b.y-t.x))

  f <- f |> arrange(desc(size), desc(distance))

  nL <- character(length(fL$id)) # final rankings
  nS <- character(length(fS$id))

  while(nrow(f) > 0) {
    t.x <- f$t.x[1]
    b.x <- f$b.x[1]
    t.y <- f$t.y[1]
    b.y <- f$b.y[1]
    if(is.na(f$size[1])) {
      p.x <- t.x
      p.y <- t.y
    }else{
      # by default they'll be at the bottom
      p.x <- b.x
      p.y <- b.y
      # tied in at least one ranking
      if(t.x != b.x || t.y != b.y) {
        # We want to place the item as deep as possible in order to delay overlap.
        # Check in what ranking can it be placed the deepest,
        # and then place it as early as possible in the other ranking.
        if(b.x < b.y || (b.x == b.y && t.x < t.y)) { # deeper in Y, or same but earlier in X
          p.x <- t.x
          p.y <- b.y
        }else { # deeper in x
          p.x <- b.x
          p.y <- t.y
        }
      }
    }
    nL[p.x] <- nS[p.y] <- f$id[1]
    if(!is.na(p.x)) {
      f <- f |>
        mutate(t.x = ifelse(t.x == p.x, t.x+1, t.x),
               b.x = ifelse(b.x == p.x, b.x-1, b.x))
    }
    if(!is.na(p.y)) {
      f <- f |>
        mutate(t.y = ifelse(t.y == p.y, t.y+1, t.y),
               b.y = ifelse(b.y == p.y, b.y-1, b.y))
    }
    f <- f |>
      slice(-1) |>
      mutate(t.o = pmax(t.x,t.y),
             b.o = pmax(b.x,b.y),
             size = b.o-t.o,
             distance = pmax(b.x-t.y, b.y-t.x))

    f <- f |> arrange(desc(size), desc(distance))

  }
  list(x = nL, y = nS)
}

