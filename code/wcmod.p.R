wcmod.p <- function(p,n,alpha=1e-10){
    l<-(-log(1-ifelse(p<alpha,alpha,p)))
    c.l<-cumsum(l)
    r<-(c(0,c.l[-length(c.l)])+c.l)/2	
    n.t<-sum(n)
    l.t=sum(l)
    t<-sqrt(sum(n))*((sum(n*r))/n.t-l.t/2)/l.t
    return(-log10(2*pnorm(-abs(t),mean=0,sd=sqrt(1/12)))*sign(-t))
}

wcmod.p2 <- function(l,n){
  c.l<-cumsum(l)
  r<-(c(0,c.l[-length(c.l)])+c.l)/2	
  n.t<-sum(n) # it gives nothing to feed this
  l.t=sum(l) # it gives nothing to feed this
  t<-sqrt(n.t)*(sum(n*r)/n.t-l.t/2)/l.t
  return(-log10(2*pnorm(-abs(t),mean=0,sd=sqrt(1/12)))*sign(-t))
}
