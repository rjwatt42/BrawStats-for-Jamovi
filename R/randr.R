rr<-function(nr,r,design=braw.def$design) {
  n<-round(
    braw.env$minN+rgamma(nr,shape=design$sNRandK,scale=(design$sN-braw.env$minN)/design$sNRandK)
  )
  r<-tanh(rnorm(nr,atanh(r),1/sqrt(n-3)))
  return(data.frame(rIV=r,nval=n))
}
