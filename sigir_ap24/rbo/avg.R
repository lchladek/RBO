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

rbo_avg_ext <- function(x, y, p){

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L
  s <- length(fS$id)
  l <- length(fL$id)


  cc <- cc.a(fS,fL)

  X <- colSums(cc$cS*cc$cL)

  d <- seq_len(l);      r <- sum(X/d * p^d)
  d <- safeseq(s+1,l);  r <- r + sum(X[s]*(d-s)/d/s * p^d)
  r <- (1-p)/p*r
  r + ((X[l]-X[s])/l + X[s]/s) * p^l
}


rbo_avg_min <- function(x, y, p) {

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L
  s <- length(fS$id)
  l <- length(fL$id)

  cc <- cc.a(fS,fL)

  X <- colSums(cc$cS*cc$cL)

  d <- seq_len(l); r <- sum((X-X[l])/d * p^d)
  r <- r - X[l]*log(1-p)
  (1-p)/p*r
}


rbo_avg_max <- function(x, y, p) {

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L
  s <- length(fS$id)
  l <- length(fL$id)

  cc <- cc.a(fS,fL)

  X <- colSums(cc$cS*cc$cL)

  f <- l+s-X[l]
  d <- seq_len(s);        r <- sum( X[d]/d * p^d )
  d <- safeseq(s+1,l);    r <- r + sum( (X[d]/d + (d-s)/d) * p^d )
  d <- safeseq(l+1,f);    r <- r + sum( (2*d-l-s+X[l])/d * p^d )
  r <- r + (p^(f+1))/(1-p)

  (1-p)/p * r
}
