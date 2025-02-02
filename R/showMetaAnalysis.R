
makeMetaHist<-function(vals,use,xlim) {
  nbins<-10
  bins<-seq(xlim[1],xlim[2],length.out=nbins+1)
  dens<-hist(vals[use],bins,plot=FALSE)$counts/length(vals)
  h<-list(bins=bins,dens=dens)
}

worldLabel<-function(metaResult,whichMeta=NULL) {
  if (is.null(whichMeta)) whichMeta<-metaResult$bestDist
    Dist<-tolower(whichMeta)
    p1<-metaResult[[Dist]]$param1Max
    p2<-metaResult[[Dist]]$param2Max

    if (is.element(Dist,c("random","fixed"))) label1<-"r[est]" else label1<-Dist
    lb<-paste0(label1,"=",brawFormat(mean(p1,na.rm=TRUE),digits=2))
    # if (length(p1)>1)
    #   lb<-paste0(lb,"\u00B1",brawFormat(std(p1),digits=2))
    if (!is.null(p2)) {
      if (is.element(Dist,c("random","fixed"))) label2<-"sd(r)[est]" else label2<-"p(null)"
      lb<-paste0(lb,"\n",label2,"=",brawFormat(mean(p2,na.rm=TRUE),digits=2))
      # if (length(p2)>1)
      #   lb<-paste0(lb,"\u00B1",brawFormat(std(p2),digits=2))
    }
    return(lb)
}

#' show a single meta-analysis 
#' 
#' @return ggplot2 object - and printed
#' @examples
#' showSingleMeta(metaResult=doMetaAnalysis(),showTheory=FALSE)
#' @export
showMetaSingle<-function(metaResult=braw.res$metaSingle,showTheory=FALSE) {
  if (is.null(metaResult)) metaResult<-doMetaAnalysis()
  
  metaAnalysis<-metaResult$metaAnalysis
  hypothesis<-metaResult$hypothesis
  design<-metaResult$design
  
  d1<-metaResult$result$rIV
  d1n<-(metaResult$result$rpIV==0)
  x<-plotAxis("rs",hypothesis)
  xlim<-x$lim
  disp1<-x$label
  # if (all(d1>=0)) xlim[1]<-0
  
  d2<-metaResult$result$nval
  y<-plotAxis("n",hypothesis)
  disp2<-y$label
  ylim<-y$lim
  
  if (y$logScale) d2<-log10(d2)
  useAll<-(d2>ylim[1]) & (d2<ylim[2])
  ptsAll<-data.frame(x=d1[useAll],y=d2[useAll])
  useNull<-(d2>ylim[1]) & (d2<ylim[2] & d1n)
  ptsNull<-data.frame(x=d1[useNull],y=d2[useNull])
  
  assign("plotArea",c(0,0,1,1),braw.env)
  g<-startPlot(xlim,ylim,
               xticks=makeTicks(x$ticks),xlabel=makeLabel(disp1),
               yticks=makeTicks(y$ticks,10^y$ticks),ylabel=makeLabel(disp2),
               top=1,g=NULL)

  if (length(d1)>=1200) {
    nbins<-diff(ylim)/(2*IQR(d2[use])*length(d2[use])^(-0.33))*2
    nbins<-min(nbins,101)
    g<-addG(g,stat_bin2d(data=pts,aes(x=x,y=y),bins=nbins)+scale_fill_gradientn(colours=c(braw.env$plotColours$graphBack,braw.env$plotColours$descriptionC)))
  }
  
  g<-drawWorld(hypothesis,design,metaResult,g,
               darken(braw.env$plotColours$metaAnalysis,off=0.1),showTheory=showTheory)
  if (metaAnalysis$includeBias) {
    rv<-seq(min(xlim),max(xlim),length.out=101)
    nv<-rw2n(rv,0.5)
    if (braw.env$nPlotScale=="log10") {nv<-log10(nv)}
    nv[nv>ylim[2] | nv<ylim[1]]<-NA
    g<-addG(g,dataLine(data.frame(x=rv,y=nv),
                       colour=darken(braw.env$plotColours$infer_sigC,off=0.1),linewidth=0.5))
  }
    
  dotSize<-braw.env$dotSize/2
  
  # show individual studies
  if (length(d1)<1200) {
    colgain<-1-min(1,sqrt(max(0,(length(d1)-50))/200))
    dotSize<-dotSize/(ceil(length(d1)/400))
    cl<-"black"
    col1<-braw.env$plotColours$metaAnalysis
    col2<-braw.env$plotColours$infer_nsigC
    g<-addG(g,dataPoint(data=ptsAll, shape=braw.env$plotShapes$study, colour = darken(col1,off=-colgain), fill = col1, size = dotSize))
    g<-addG(g,dataPoint(data=ptsNull,shape=braw.env$plotShapes$study, colour = darken(col2,off=-colgain), fill = col2,  size = dotSize))
  }
  
  lb<-worldLabel(metaResult,metaAnalysis$analysisType)
  names=strsplit(lb,"\n")[[1]]
  if (length(names)==1) colours=braw.env$plotColours$metaAnalysis else colours=c(braw.env$plotColours$metaAnalysis,NA)
  g<-addG(g,dataLegend(data.frame(names=names,colours=colours),title="",shape=22))
  # g<-addG(g,plotTitle(lb,"left",size=1))
  
  if (braw.env$graphHTML && braw.env$autoShow) {
    showHTML(g)
    return(invisible(g))
  }
  else return(g)  

}

#' show a multiple meta-analyses
#' 
#' @return ggplot2 object - and printed
#' @examples
#' showMultipleMeta<-function(metaResult=doMetaAnalysis(),showType="n-k")
#' @export
showMetaMultiple<-function(metaResult=braw.res$metaMultiple,showType="n-k") {
  if (is.null(metaResult)) metaResult<-doMetaMultiple()

  if (metaResult$metaAnalysis$analysisType=="fixed") {
    if (metaResult$hypothesis$effect$world$worldOn) showType<-"k-rp"
    else showType<-"S-k"
  }
  if (metaResult$metaAnalysis$analysisType=="random") showType<-"sd-k"
  if (showType=="S-S") {
    braw.env$plotArea<-c(0,0,1,1)
    g<-drawMeta(metaResult=metaResult,showType=showType,g=NULL)
  } else {
    g<-nullPlot()
    if (is.element(metaResult$metaAnalysis$analysisType,c("fixed","random"))) {
      braw.env$plotArea<-c(0,0,1,1)
      g<-drawMeta(metaResult=metaResult,showType=showType,whichMeta=metaResult$metaAnalysis$analysisType,g=g)
    }  else {
      braw.env$plotArea<-c(0,0,0.4,1)
      if (!all(is.na(metaResult$single$Smax)))
        g<-drawMeta(metaResult=metaResult,whichMeta="Single",showType=showType,g)
      braw.env$plotArea<-c(0.4,0,0.3,1)
      if (!all(is.na(metaResult$gauss$Smax)))
        g<-drawMeta(metaResult=metaResult,whichMeta="Gauss",showType=showType,g)
      braw.env$plotArea<-c(0.7,0,0.3,1)
      if (!all(is.na(metaResult$exp$Smax)))
        g<-drawMeta(metaResult=metaResult,whichMeta="Exp",showType=showType,g)
    }
  }
  if (braw.env$graphHTML && braw.env$autoShow) {
    showHTML(g)
    return(invisible(g))
  }
  else return(g)  
}

drawMeta<-function(metaResult=doMetaMultiple(),whichMeta="Single",showType="n-k",g=NULL) {
  
  metaAnalysis<-metaResult$metaAnalysis

  xlim<-c(-1,1)
  xticks<-seq(-1,1,0.5)
  
  if (is.element(whichMeta,c("Single","Gauss","Exp"))) {
  n1<-sum(metaResult$bestDist=="Single")
  n2<-sum(metaResult$bestDist=="Gauss")
  n3<-sum(metaResult$bestDist=="Exp")
  sAll<-c(metaResult$single$Smax,metaResult$gauss$Smax,metaResult$exp$Smax)
  
  if (showType=="S-S") {
    use<-order(c(n1,n2,n3))
    use1<-c("Single","Gauss","Exp")[use[2]]
    use2<-c("Single","Gauss","Exp")[use[3]]
    metaX<-metaResult[[tolower(use1)]]
    metaY<-metaResult[[tolower(use2)]]
    x<-metaX$Smax
    yS<-metaY$Smax
    y1<-yS
    xticks<-c()
  } else {
    switch (whichMeta,
            "Single"={
              x<-metaResult$single$param1Max
              yS<-metaResult$single$Smax
              y1<-metaResult$single$param2Max
            },
            "Gauss"={
              x<-metaResult$gauss$param1Max
              yS<-metaResult$gauss$Smax
              y1<-metaResult$gauss$param2Max
            },
            "Exp"={
              x<-metaResult$exp$param1Max
              yS<-metaResult$exp$Smax
              y1<-metaResult$exp$param2Max
            }
    )
  }
    keep<- !is.na(x) & !is.na(yS)
    best<-metaResult$bestS[keep]
    yS<-yS[keep]
    y1<-y1[keep]
    x<-x[keep]
    useBest<-yS==best
    
    if (isempty(x)) {return(nullPlot())}
  }
  
    yticks<-c()
    switch (showType,
            "S-k"={
              x<-metaResult$fixed$param1Max
              y<-metaResult$fixed$Smax
              y1<-0
              sAll<-metaResult$fixed$Smax
              ylim<-c(min(sAll,na.rm=TRUE),max(sAll,na.rm=TRUE))+c(-1,1)*(max(sAll,na.rm=TRUE)-min(sAll,na.rm=TRUE))/4
              ylabel<-"log(lk)"
              xlabel<-"r[est]"
              useBest<-1:length(x)
            },
            "n-k"={
              y<-y1
              ylim<-c(-0.02,1.1)
              ylabel<-"p[null]"
              xlabel<-braw.env$Llabel
            },
            "k-rp"={
              x<-metaResult$fixed$rpIV
              y<-metaResult$fixed$param1Max
              y1<-0
              ylim<-c(-1,1)
              ylabel<-"r[est]"
              xlabel<-"r[p]"
              useBest<-1:length(x)
            },
            "sd-k"={
              x<-metaResult$random$param1Max
              y<-metaResult$random$param2Max
              y1<-0
              ylim<-c(-0.02,0.5)
              ylabel<-"sd(r)[est]"
              xlabel<-"r[est]"
              useBest<-1:length(x)
            },
            "S-S"={
              y<-yS
              xlim<-c(min(sAll,na.rm=TRUE),max(sAll,na.rm=TRUE))
              if (length(x)==1) xlim<-xlim+c(-1,1)
              else xlim=xlim+c(-1,1)*(max(sAll,na.rm=TRUE)-min(sAll,na.rm=TRUE))/4
              xlabel<-paste0("log(lk ",use1,")")
              ylim<-xlim
              ylabel<-paste0("log(lk ",use2,")")
              useBest<- (y>x & metaResult$hypothesis$effect$world$populationPDF==use2) | (y<x & metaResult$hypothesis$effect$world$populationPDF==use1)
            }
    )
    pts<-data.frame(x=x,y=y)
      
      if (braw.env$plotArea[1]==0)  
        g<-startPlot(xlim,ylim,
                     xticks=makeTicks(xticks),xlabel=makeLabel(xlabel),
                     yticks=makeTicks(yticks),ylabel=makeLabel(ylabel),
                     top=TRUE,g=g)
    else  
        g<-startPlot(xlim,ylim,
                     xticks=makeTicks(xticks),xlabel=makeLabel(xlabel),
                     top=TRUE,g=g)

      dotSize=16*min(0.25,2.5/sqrt(length(x)))
      
      g<-addG(g,dataPoint(data=pts,shape=braw.env$plotShapes$meta, 
                          colour="black", fill="grey", alpha=min(1,2.5/sqrt(length(x))), 
                          size = dotSize))
      pts<-data.frame(x=x[useBest],y=y[useBest])
      g<-addG(g,dataPoint(data=pts,shape=braw.env$plotShapes$meta,
                          colour="black", fill=braw.env$plotColours$metaAnalysis, alpha=min(1,2.5/sqrt(length(x))), 
                          size = dotSize))
      
      if (showType=="S-S") {
        g<-addG(g,dataPath(data=data.frame(x=xlim,y=ylim),colour="red"))
      }

    if (mean(y1)>0.5) {
      yp<-ylim[1]+diff(ylim)/10
      vj<-0
    } else {
        yp<-ylim[2]-diff(ylim)/10
        vj<-1
        }
    if (showType=="S-S") {
      fullText<-paste0(use2,"(",format(mean(metaY$param1Max),digits=3))
      if (length(metaY$param1Max)>1) fullText<-paste0(fullText,"\u00B1",format(std(metaY$param1Max),digits=2),")")
      else fullText<-paste0(fullText,")")
      if (metaAnalysis$includeNulls) {
        fullText<-paste0(fullText,"\nnull=",format(mean(metaY$param2Max),digits=3))
        if (length(metaY$param2Max)>1) fullText<-paste0(fullText,"\u00B1",format(std(metaY$param2Max),digits=2),")")
      }
      fullText<-paste0(fullText,"\nS= ",format(mean(metaY$Smax),digits=2))
      if (length(metaY$Smax)>1) fullText<-paste0(fullText,"\u00B1",format(std(metaY$Smax),digits=2),")")
      fullText<-paste0(fullText," (",format(sum(y>x)),"/",length(metaResult$bestDist),")")
      
      if (mean(y>x)) colM=braw.env$plotColours$metaAnalysis  else colM="grey"
      names<-strsplit(fullText,"\n")[[1]]
      g<-addG(g,dataLegend(data.frame(names=names,colours=c(colM,rep(NA,length(names)-1))),title="",shape=braw.env$plotShapes$meta))
      
      fullText<-paste0(use1,"(",format(mean(metaX$param1Max),digits=3))
      if (length(metaX$param1Max)>1) fullText<-paste0(fullText,"\u00B1",format(std(metaX$param1Max),digits=2),")")
      else fullText<-paste0(fullText,")")
      if (metaAnalysis$includeNulls) {
        fullText<-paste0(fullText,"\nnull=",format(mean(metaX$param2Max),digits=3))
        if (length(metaX$param2Max)>1) fullText<-paste0(fullText,"\u00B1",format(std(metaX$param2Max),digits=2),")")
      }
      fullText<-paste0(fullText,"\nS= ",format(mean(metaX$Smax),digits=2))
      if (length(metaX$Smax)>1) fullText<-paste0(fullText,"\u00B1",format(std(metaX$Smax),digits=2),")")
      fullText<-paste0(fullText," (",format(sum(x>y)),"/",length(metaResult$bestDist),")")
      
      if (mean(y>x)) colM="grey"  else colM=braw.env$plotColours$metaAnalysis
      names<-strsplit(fullText,"\n")[[1]]
      g<-addG(g,dataLegend(data.frame(names=names,colours=c(colM,rep(NA,length(names)-1))),title="",shape=braw.env$plotShapes$meta))
    } else {

      if (showType!="k-rp") {
      if (is.element(showType,c("S-k","sd-k"))) {
        colM=braw.env$plotColours$metaAnalysis
        lb<-worldLabel(metaResult,whichMeta)
        names<-strsplit(lb,"\n")[[1]]

        if (whichMeta=="fixed") {names<-names[1]; colours<-colM;}
        else colours<-c(colM,NA)
        g<-addG(g,dataLegend(data.frame(names=names,colours=colours),title="",
                             shape=braw.env$plotShapes$meta))
      } else {
        use<-which.max(c(n1,n2,n3))
        bestD<-c("Single","Gauss","Exp")[use]
        if (whichMeta==bestD)  colM=braw.env$plotColours$metaAnalysis else colM="grey"
        lb<-worldLabel(metaResult,whichMeta)
        g<-addG(g,dataLegend(data.frame(names=strsplit(lb,"\n")[[1]],colours=c(colM,NA)),title="",shape=braw.env$plotShapes$meta))
      }
      }
    }
    return(g)

}

makeWorldDist<-function(metaResult,design,world,x,y,sigOnly=FALSE,doTheory=FALSE) {
  if (doTheory) {
    lambda<-world$populationPDFk
    nullP<-world$populationNullp
    offset<-0
    shape<-0
    nullP<-world$populationNullp
    if (metaResult$metaAnalysis$analysisType=="random") {
      lambda<-metaResult$hypothesis$effect$rSD
      offset<-metaResult$hypothesis$effect$rIV
      nullP<-0
      world$populationPDF<-"Gauss"
    }
    if (metaResult$metaAnalysis$analysisType=="fixed") {
      lambda<-metaResult$hypothesis$effect$rIV
      nullP<-0
      world$populationPDF<-"Single"
    }
  } else {
    lambda<-metaResult$bestK
    nullP<-metaResult$bestNull
    offset<-0
    shape<-0
    if (metaResult$metaAnalysis$analysisType=="random") {
      lambda<-metaResult$random$param1Max
      shape<-metaResult$random$param2Max
      nullP<-0
      world$populationPDF<-"Single"
    }
    if (metaResult$metaAnalysis$analysisType=="fixed") {
      lambda<-metaResult$fixed$param1Max
      nullP<-0
      world$populationPDF<-"Single"
    }
  }
  sigma<-1/sqrt(y-3)
  gain<-dgamma(y-braw.env$minN,shape=design$sNRandK,scale=(design$sN-braw.env$minN)/design$sNRandK)
  
  z<-c()
  switch (world$populationPDF,
          "Single"={
            for (i in 1:length(y)) {
              z1<-SingleSamplingPDF(atanh(x),lambda,sigma[i],shape)$pdf*(1-nullP)+
                SingleSamplingPDF(atanh(x),0,sigma[i])$pdf*nullP
              if (sigOnly) {
                zcrit<-qnorm(1-braw.env$alphaSig/2,0,sigma[i])
                z1[atanh(abs(x))<zcrit]<-0
              }
              z<-rbind(z,zdens2rdens(z1,x)*gain[i])
            }
          },
          "Gauss"={
            for (i in 1:length(y)) {
              z1<-GaussSamplingPDF(atanh(x),lambda,sigma[i],offset)$pdf*(1-nullP)+
                SingleSamplingPDF(atanh(x),0,sigma[i])$pdf*nullP
              if (sigOnly) {
                zcrit<-qnorm(1-braw.env$alphaSig/2,0,sigma[i])
                z1[atanh(abs(x))<zcrit]<-0
              }
              z<-rbind(z,zdens2rdens(z1,x)*gain[i])
            }
          },
          "Exp"={
            for (i in 1:length(y)) {
              z1<-ExpSamplingPDF(atanh(x),lambda,sigma[i])$pdf*(1-nullP)+
                SingleSamplingPDF(atanh(x),0,sigma[i])$pdf*nullP
              if (sigOnly) {
                zcrit<-qnorm(1-braw.env$alphaSig/2,0,sigma[i])
                z1[atanh(abs(x))<zcrit]<-0
              }
              z<-rbind(z,zdens2rdens(z1,x)*gain[i])
            }
          }
  )
  return(z)
}

drawWorld<-function(hypothesis,design,metaResult,g,colour="white",showTheory=FALSE) {
  world<-hypothesis$effect$world
  if (!world$worldOn) {
    world<-makeWorld(worldOn=TRUE,populationPDF="Single",populationRZ="r",
                     populationPDFk=hypothesis$effect$rIV,populationNullp=0)
  }
  x<-seq(-1,1,length.out=101)*braw.env$r_range
  y<-seq(5,braw.env$maxN,length.out=101)

  za<-makeWorldDist(metaResult,design,world,x,y,doTheory=TRUE)

  if (!is.element(metaResult$metaAnalysis$analysisType,c("fixed","random"))) {
    world$populationPDF<-metaResult$bestDist
    world$populationPDFk<-metaResult$bestK
    world$populationNullp<-metaResult$bestNull
  }
  zb<-makeWorldDist(metaResult,design,world,x,y,doTheory=FALSE)
  
  if (braw.env$nPlotScale=="log10") {y<-log10(y)}
  
  x1<-as.vector(matrix(rep(x,101),101,byrow=TRUE))
  y1<-as.vector(matrix(rep(y,101),101,byrow=FALSE))
  za<-as.vector(za)
  zb<-as.vector(zb)

  use<-is.finite(za) & y1<braw.env$maxRandN*design$sN
  za<-za/max(za,na.rm=TRUE)
  ptsa<-data.frame(x=x1[use],y=y1[use],z=za[use])
  
  use<-is.finite(zb) & y1<braw.env$maxRandN*design$sN
  zb<-zb/max(zb,na.rm=TRUE)
  ptsb<-data.frame(x=x1[use],y=y1[use],z=zb[use])

  # white is the actual world
  # black is the best fit world
  if (showTheory) {
    g<-addG(g,dataContour(data=ptsa,colour="white"))
  }
  g<-addG(g,dataContour(data=ptsb,colour=colour,linewidth=0.5))
  return(g)
}
