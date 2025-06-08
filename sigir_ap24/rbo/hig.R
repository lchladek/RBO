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

hig_rankings <- function(x,y){

  # If no ties, do nothing
  if(all(lengths(x) == 1) && all(lengths(y) == 1))
    return(list(x = x, y = y))

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L

  fS <- fS |> bind_cols()
  fL <- fL |> bind_cols()

  s <- length(fS$id)
  l <- length(fL$id)


  # still unmatched items at depth d
  uL <- uS <- character(0)

  for(d in safeseq(1,l)) {
    # what items will be fixed at depth d
    fixL <- fixS <- NA
    # candidates at depth d
    iL <- fL |> filter(t <=d & b >= d) |> pull(id) #those that can be active
    iS <- fS |> filter(t <=d & b >= d) |> pull(id) # items at d or crossed by d

    if(length(iL) > 1) { #we are in a tied group
      # if more than 1 candidate, check whether any one matches the unmatched in S
      cL <- iL[iL %in% uS] #if any of the tied el matches one in the unmatched ones
      if(length(cL) > 0)
        fixL <- cL[1] #arbitrarly choose one of them
    }else{
      # if just one candidate (ie. untied), just fix it
      fixL <- iL
    }
    if(length(iS) > 1) {
      # if more than 1 candidate, check whether any one matches the unmatched in L
      cS <- iS[iS %in% uL]
      if(length(cS) > 0)
        fixS <- cS[1]
    }else{
      # if just one candidate (ie. untied), just fix it
      fixS <- iS
    }

    # if still unknown in both L and S, check whether we could fix the same item in both
    if(d <= s && is.na(fixL) && is.na(fixS)) {
      cLS <- intersect(iL, iS)
      if(length(cLS) > 0) {
        fixS <- fixL <- cLS[1]
      }
    }
    # if still unknown in L, check whether we can match what will be fixed in S
    if(is.na(fixL)){
      if(d <= s && fixS %in% iL)
        fixL <- fixS
      else
        fixL <- iL[1]
    }
    # if still unknown in S, check whether we can match what will be fixed in L
    if(is.na(fixS) && d <= s) {
      if(fixL %in% iS)
        fixS <- fixL
      else
        fixS <- iS[1]
    }

    # Update unmatched
    if(!fixL %in% c(uS, fixS))
      uL <- c(uL, fixL)
    if(d <= s && !fixS %in% c(uL, fixL))
      uS <- c(uS, fixS)

    uL <- setdiff(uL, fixS)
    uS <- setdiff(uS, fixL)

    # Do fix
    i <- which(fL$id == fixL)
    fL$t[fL$t == d] <- d+1
    fL$t[i] <- fL$b[i] <- d

    if(d <= s) {
      i <- which(fS$id == fixS)
      fS$t[fS$t == d] <- d+1
      fS$t[i] <- fS$b[i] <- d
    }
  }
  list(x = fL |> arrange(t) |> pull(id),
       y = fS |> arrange(t) |> pull(id))
}
