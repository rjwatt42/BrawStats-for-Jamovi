
# META-ANALYSIS
# calculations
# graphs (sample, describe, infer)
# report (sample, describe, infer)
#    

#' @return metaResult object 
#' @examples
#' doMetaAnalysis<-function(nsims=50,metaAnalysis=makeMetaAnalysis(),
#'                          hypothesis=braw.def$hypothesis,design=braw.def$design,evidence=braw.def$evidence,
#'                          metaResult=NULL)
#' @export
doMetaAnalysis<-function(nsims=100,metaResult=braw.res$metaMultiple,metaAnalysis=braw.def$metaAnalysis,
                         hypothesis=braw.def$hypothesis,design=braw.def$design,evidence=braw.def$evidence
                         ) {
  if (is.null(metaAnalysis)) metaAnalysis<-makeMetaAnalysis()
  if (nsims==1) metaResult<-NULL
  evidence$sigOnly<-metaAnalysis$sigOnlySource
  
  for (i in 1:nsims) {
    studies<-multipleAnalysis(metaAnalysis$nstudies,hypothesis,design,evidence)
    metaResult<-runMetaAnalysis(metaAnalysis,studies,metaResult)
  }
  metaResult$hypothesis<-hypothesis
  metaResult$design<-design
  setBrawRes("metaResult",metaResult)
  if (nsims>1) setBrawRes("metaMultiple",metaResult)
  else         setBrawRes("metaSingle",metaResult)
  metaResult
}


getMaxLikelihood<-function(zs,ns,df1,dist,metaAnalysis) {
  # param1 is kvals
  # param2 is normally nullvals
  
  np1points<-13
  np2points<-13
  
  niterations<-1
  reInc<-(np1points-1)/2/3
  
  if (metaAnalysis$includeNulls) {
    param2<-seq(0,1,length.out=np2points)
  } else {
    param2<-0
  }
  if (dist=="Single") {
    param1<-seq(-1,1,length.out=np1points)*0.95
  } else {
    param1<-seq(0.01,1,length.out=np1points)
  }
  
  if (dist=="fixed") {
    param1<-seq(-1,1,length.out=np1points)*0.95
    param0<-0
  }
  if (dist=="random") {
    param1<-seq(-1,1,length.out=np1points)*0.95
    param2<-seq(0,0.5,length.out=np2points)*0.95
  }
  
  remove_nonsig<-metaAnalysis$includeBias
  for (re in 1:niterations) {
    S<-getLogLikelihood(zs,ns,df1,dist,param1,param2,remove_nonsig)
    Smax<-max(S,na.rm=TRUE)
    
    use<-which(S==Smax, arr.ind = TRUE)
    param2Max<-param2[use[1,2]]
    lb2<-param2[max(1,use[1,2]-reInc)]
    ub2<-param2[min(length(param2),use[1,2]+reInc)]
    param1Max<-param1[use[1,1]]
    lb1<-param1[max(1,use[1,1]-reInc)]
    ub1<-param1[min(length(param1),use[1,1]+reInc)]
    
    param1<-seq(lb1,ub1,length.out=np1points)
    if (metaAnalysis$includeNulls) {
      param2<-seq(lb2,ub2,length.out=np2points)
    }
  }
  
  if (niterations==1) {      
    fun<-function(x) { -getLogLikelihood(zs,ns,df1,dist,x[1],x[2],remove_nonsig)}
    result<-fmincon(c(param1Max,param2Max),fun,ub=c(ub1,ub2),lb=c(lb1,lb2))
    param1Max<-result$par[1]
    param2Max<-result$par[2]
    Smax<- -result$value
  }
  return(list(param1Max=param1Max,param2Max=param2Max,Smax=Smax))
}


runMetaAnalysis<-function(metaAnalysis,studies,metaResult){
  rs<-studies$rIV
  zs<-atanh(rs)
  ns<-studies$nval
  df1<-studies$df1
  
  fixed<-list(param1Max=NA,param2Max=NA,Smax=NA)
  random<-list(param1Max=NA,param2Max=NA,Smax=NA)
  single<-list(param1Max=NA,param2Max=NA,Smax=NA)
  gauss<-list(param1Max=NA,param2Max=NA,Smax=NA)
  exp<-list(param1Max=NA,param2Max=NA,Smax=NA)
  switch(metaAnalysis$analysisType,
         "fixed"={
           # a fixed analysis finds a single effect size
           metaAnalysis$includeNulls<-FALSE
           fixed<-getMaxLikelihood(zs,ns,df1,"fixed",metaAnalysis)
         },
         "random"={
           metaAnalysis$includeNulls<-FALSE
           random<-getMaxLikelihood(zs,ns,df1,"random",metaAnalysis)
         },
         "world"={
           # doing world effects analysis
           
           # find best Single 
           if (metaAnalysis$modelPDF=="Single" || (metaAnalysis$modelPDF=="All")) 
             single<-getMaxLikelihood(zs,ns,df1,"Single",metaAnalysis)

           # find best Gauss
           if (metaAnalysis$modelPDF=="Gauss" || metaAnalysis$modelPDF=="All") 
             gauss<-getMaxLikelihood(zs,ns,df1,"Gauss",metaAnalysis)

           # find best Exp
           if (metaAnalysis$modelPDF=="Exp" || metaAnalysis$modelPDF=="All") 
             exp<-getMaxLikelihood(zs,ns,df1,"Exp",metaAnalysis)

         })
if (is.element(metaAnalysis$analysisType,c("fixed","random"))) {
  if (metaAnalysis$analysisType=="fixed") {
    fixed$param1Max<-c(metaResult$fixed$param1Max,fixed$param1Max)
    fixed$Smax<-c(metaResult$fixed$Smax,fixed$Smax)
    fixed$param2Max<-c(metaResult$fixed$param2Max,fixed$param2Max)
    count<-length(fixed$Smax)
  } else {
    random$param1Max<-c(metaResult$random$param1Max,random$param1Max)
    random$Smax<-c(metaResult$random$Smax,random$Smax)
    random$param2Max<-c(metaResult$random$param2Max,random$param2Max)
    count<-length(random$Smax)
  }  
  bestDist<-NA
  bestK<-NA
  bestNull<-NA
  bestS<-NA
  
} else {
  use<-which.max(c(single$Smax,gauss$Smax,exp$Smax))
  bestDist<-c("Single","Gauss","Exp")[use]
  bestK<-c(single$param1Max,gauss$param1Max,exp$param1Max)[use]
  bestNull<-c(single$param2Max,gauss$param2Max,exp$param2Max)[use]
  bestS<-c(single$Smax,gauss$Smax,exp$Smax)[use]
  
  if (!is.null(metaResult)) {
    bestDist<-c(metaResult$bestDist,bestDist)
    bestK<-c(metaResult$bestK,bestK)
    bestNull<-c(metaResult$bestNull,bestNull)
    bestS<-c(metaResult$bestS,bestS)
    count<-length(bestS)
    single$param1Max<-c(metaResult$single$param1Max,single$param1Max)
    single$Smax<-c(metaResult$single$Smax,single$Smax)
    single$param2Max<-c(metaResult$single$param2Max,single$param2Max)
    gauss$param1Max<-c(metaResult$gauss$param1Max,gauss$param1Max)
    gauss$Smax<-c(metaResult$gauss$Smax,gauss$Smax)
    gauss$param2Max<-c(metaResult$gauss$param2Max,gauss$param2Max)
    exp$param1Max<-c(metaResult$exp$param1Max,exp$param1Max)
    exp$Smax<-c(metaResult$exp$Smax,exp$Smax)
    exp$param2Max<-c(metaResult$exp$param2Max,exp$param2Max)
  }
  
}
  
  metaResult<-list(fixed=fixed,
                   random=random,
                   single=single,
                   gauss=gauss,
                   exp=exp,
                   bestDist=bestDist,
                   bestK=bestK,
                   bestNull=bestNull,
                   bestS=bestS,
                   count=count,
                   metaAnalysis=metaAnalysis,
                   result=studies
  )
  return(metaResult)
}

