# lambda3 changed from regmex in that nop is served and lambda is lifted to the seq length
lambda3 <- function (pattern, seq, nop, freq = FALSE, nt.null = 1, mode = FALSE,overlap = FALSE) 
{
  tm <- Regmex:::transition.matrix(pattern$matrix, seq$freq.mono)
  #nop <- Regmex:::n.obs.pat(pattern$pattern, seq$seq, overlap = overlap)
  #nop <- ifelse(nop %in% c(0, 1), 2, nop)
  finl.st <- pattern$endState
  n.states <- dim(tm)[1]
  lam <- matrix(0, n.states * nop, n.states * nop)
  for (i in 1:nop) {
    bgn.indx <- (i - 1) * n.states + 1
    end.indx <- (i * n.states)
    blck.indx <- bgn.indx:end.indx
    lam[blck.indx, blck.indx] <- tm
  }
  for (k in finl.st) {
    for (i in 1:(nop - 1)) {
      bgn.finl <- (i - 1) * n.states + k
      bgn.indx <- (i - 1) * n.states + 1
      end.indx <- (i * n.states)
      blck.indx <- bgn.indx:end.indx
      lam[bgn.finl, blck.indx + n.states] <- lam[bgn.finl,blck.indx]
      lam[bgn.finl, blck.indx] <- 0
    }
  }
  end.st <- n.states * (nop - 1) + finl.st
  for (i in end.st) {
    lam[i, ] <- 0
    lam[i, i] <- 1
  }
  #    return(lam)
  return(lam %^% seq$length)
}
#pd
pd <- function(pattern, seq, maxnop=10){
  nop <- Regmex:::n.obs.pat(pattern$pattern, seq$seq, overlap = FALSE)
  nop <- min(nop,maxnop)
  if (nop == 0) return(list(prob.n.or.more = 1, n.obs.patterns = nop, prob.dist = NA, prob.1.or.more = NA))
  nop2 <- ifelse(nop %in% c(0,1),2,nop)
  prb.init.st <- lambda3(pattern, seq,nop=nop2)[pattern$startState, ]
  st.st <- pattern$startState
  end.st <- pattern$endState
  nSt <- dim(pattern$matrix)[1]
  prb.dst <- rep(0, nop2 + 1)
  for (i in 1:nop2) {
    if (i == 1) {prb.dst[i] <- sum(prb.init.st[1:nSt][-end.st])
    if (nop == 1){return(list(prob.n.or.more = 1 - prb.dst[1], n.obs.patterns = nop, prob.dist = prb.dst[1], prob.1.or.more = 1 - prb.dst[1]))}    
    next}
    prb.dst[i] <- sum(c(prb.init.st[(nSt * (i - 1) + 1):(nSt*i)][-end.st], prb.init.st[nSt * (i - 2) + end.st]))
  }
  prb.dst[nop2 + 1] <- sum(prb.init.st[(nSt * (nop2 - 1) + end.st)]) # the final end state
  pnom <- prb.dst[nop + 1]
  return(list(prob.n.or.more = pnom, n.obs.patterns = nop, prob.dist = prb.dst, prob.1.or.more = 1 - prb.dst[1]))
}
pd.mrs <- function(pattern, seq){
  prb.init.st <- lambda3(pattern, seq,nop=1)[pattern$startState, ]
  st.st <- pattern$startState
  end.st <- pattern$endState
  nSt <- dim(pattern$matrix)[1]
  prb.dst <- sum(prb.init.st[1:nSt][-end.st])
  return(list(prob.n.or.more = 1 - prb.dst, n.obs.patterns = 1, prob.dist = prb.dst, prob.1.or.more = 1 - prb.dst))
}

pd.mrs2 <- function(pattern, seq){
  tm <- Regmex:::transition.matrix(pattern$matrix, seq$freq.mono)
  finl.st <- pattern$endState
  tm[finl.st,] <- 0
  tm[finl.st,finl.st] <- 1
  return(1-sum((tm %^% seq$length)[pattern$startState,-finl.st]))
}