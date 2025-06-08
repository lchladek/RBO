####################################################################################################
# Modified code from: M. Corsi and J. Urbano, "The Treatment of Ties in Rank-Biased Overlap",      #
# International ACM SIGIR, Conference on Research and Development in Information Retrieval, 2024.  #
#                                                                                                  #
# For efficiency in the experiments, we vectorized function rbo over the persistence parameter p.  #
#                                                                                                  #
# Original code is as in rbo/rbo.R                                                                 #
####################################################################################################

source("rbo/rbo.R")

rbo <- function(x, y, pp, ties = c("a", "b", "w"), score = c("ext", "min", "max", "res")) {
  #check_rbo(x, y, pp)

  score <- match.arg(score, several.ok = TRUE)
  ties <- match.arg(ties)

  # Flat representations
  sl <- flatten(x, y)
  S <- sl$S
  L <- sl$L
  s <- length(S$id)
  l <- length(L$id)

  # If there are no ties, use w-variant for efficiency
  if(all(S$t == S$b) && all(L$t == L$b))
    ties <- "w"

  # Calculate individual item contributions
  if(ties == "w")
    CC <- cc.w(S, L)
  else
    CC <- cc.a(S, L)

  # First section: 1 to s
  d <- safeseq(1, s)
  X1 <- colSums(CC$cS[,d,drop=FALSE] * CC$cL[,d,drop=FALSE])
  if(ties == "w")
    A1 <- 2*X1 / ( colSums(CC$cS[,d,drop=FALSE]) + colSums(CC$cL[,d,drop=FALSE]) )
  else if(ties == "a")
    A1 <- X1 / d
  else
    A1 <- X1 / sqrt(colSums(CC$cS[,d,drop=FALSE]^2)) / sqrt(colSums(CC$cL[,d,drop=FALSE]^2))

  # Second section: s+1 to l
  d <- safeseq(s+1, l)
  X2_seen <- colSums(CC$cS[,d,drop=FALSE] * CC$cL[,d,drop=FALSE])

  if("min" %in% score || "res" %in% score) {
    X2.min <- X2_seen
    if(ties == "w") {
      A2.min <- 2*X2.min / ( d + colSums(CC$cL[,d,drop=FALSE]) )
    }else if(ties == "a")
      A2.min <- X2.min / d
    else
      A2.min <- X2.min / sqrt(d) / sqrt(colSums(CC$cL[,d,drop=FALSE]^2))
  }
  if("max" %in% score || "res" %in% score) {
    LdS <- setdiff(L$id, S$id) # uniques in L
    cL_unseen.max <- CC$cL[LdS[d-s],d,drop=FALSE]
    cL_unseen.max[lower.tri(cL_unseen.max)] <- NA
    X2_unseen.max <- colSums(cL_unseen.max, na.rm = TRUE)
    X2.max <- X2_seen + X2_unseen.max
    if(ties == "w") {
      A2.max <- 2*X2.max / ( d + colSums(CC$cL[,d,drop=FALSE]) )
    }else if(ties == "a")
      A2.max <- X2.max / d
    else
      A2.max <- X2.max / sqrt(d) / sqrt(colSums(CC$cL[,d,drop=FALSE]^2))
  }
  if("ext" %in% score) {
    LdS <- setdiff(L$id, S$id) # uniques in L
    cL_unseen <- CC$cL[LdS,d,drop=FALSE]
    cL_unseen[cL_unseen==0] <- NA
    cL_unseen <- colMeans(cL_unseen, na.rm = TRUE)
    X2_unseen.ext <- (d-s)*A1[s]*cL_unseen

    X2.ext <- X2_seen + X2_unseen.ext
    if(ties == "w") {
      A2.ext <- 2*X2.ext / ( d + colSums(CC$cL[,d,drop=FALSE]) )
    }else if(ties == "a")
      A2.ext <- X2.ext / d
    else
      A2.ext <- X2.ext / sqrt(d) / sqrt(colSums(CC$cL[,d,drop=FALSE]^2))
  }
  lapply(pp, function (p){
    # Third section: l+1 to inf
    X_seen <- c(X1, X2_seen)
    Xl <- X_seen[l] # X_l
    if("min" %in% score || "res" %in% score) {
      d <- safeseq(1, l)
      sec3.min <- Xl*( log(1/(1-p)) - sum(p^d / d) )
    }
    if("max" %in% score || "res" %in% score) {
      f <- l+s-Xl
      d <- safeseq(l+1,f)
      sec3.max <- sum((2*d-l-s+Xl)/d * p^d) + p^(f+1) / (1-p)
    }
    if("ext" %in% score) {
      sec3.ext <- (Xl + (l-s)*A1[s])/l * p^(l+1) / (1-p)
    }

    # All sections
    d <- safeseq(1,l)
    if("min" %in% score || "res" %in% score)
      rbo.min <- (1-p)/p*( sum(c(A1, A2.min)*p^d) + sec3.min )
    if("max" %in% score || "res" %in% score)
      rbo.max <- (1-p)/p*( sum(c(A1, A2.max)*p^d) + sec3.max )
    if("res" %in% score)
      rbo.res <- rbo.max-rbo.min
    if("ext" %in% score)
      rbo.ext <- (1-p)/p*( sum(c(A1, A2.ext)*p^d) + sec3.ext )

    r <- sapply(score, function(s) get(paste0("rbo.",s)))
    r <- pmin(1,pmax(0,r)) # for numerical stability
    names(r) <- score
    r
  }) |> bind_rows() |> mutate(p = pp, .before=1)
}


####################################################################################################
# Modified code from rbo/avg.R.                                                                    #
#                                                                                                  #
# For efficiency in the experiments, we vectorized function rbo over the persistence parameter p.  #
####################################################################################################

rbo_avg_ext <- function(x,y,pp){

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L
  s <- length(fS$id)
  l <- length(fL$id)


  cc <- cc.a(fS,fL)

  X <- colSums(cc$cS*cc$cL)

  sapply(pp, function(p) {
    d <- seq_len(l);      r <- sum(X/d * p^d)
    d <- safeseq(s+1,l);  r <- r + sum(X[s]*(d-s)/d/s * p^d)
    r <- (1-p)/p*r
    r + ((X[l]-X[s])/l + X[s]/s) * p^l
  })
}


rbo_avg_min <- function(x, y, pp) {

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L
  s <- length(fS$id)
  l <- length(fL$id)

  cc <- cc.a(fS,fL)

  X <- colSums(cc$cS*cc$cL)

  sapply(pp, function(p) {
    d <- seq_len(l); r <- sum((X-X[l])/d * p^d)
    r <- r - X[l]*log(1-p)
    (1-p)/p*r
  })
}


rbo_avg_max <- function(x, y, pp) {

  f <- flatten(x,y)
  fS <- f$S
  fL <- f$L
  s <- length(fS$id)
  l <- length(fL$id)

  cc <- cc.a(fS,fL)

  X <- colSums(cc$cS*cc$cL)

  sapply(pp, function(p){
    f <- l+s-X[l]
    d <- seq_len(s);        r <- sum( X[d]/d * p^d )
    d <- safeseq(s+1,l);    r <- r + sum( (X[d]/d + (d-s)/d) * p^d )
    d <- safeseq(l+1,f);    r <- r + sum( (2*d-l-s+X[l])/d * p^d )
    r <- r + (p^(f+1))/(1-p)

    (1-p)/p * r
  })
}
