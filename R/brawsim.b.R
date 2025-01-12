
# This file is a generated template, your changes will not be overwritten

BrawSimClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "BrawSimClass",
  inherit = BrawSimBase,
  private = list(
    # .htmlwidget = NULL,  # Add instance for HTMLWidget
    
    .init = function() {
      # private$.htmlwidget <- HTMLWidget$new() # Initialize the HTMLWidget instance 
      
      # initialization code 
      if (!exists("braw.env")) {
        BrawOpts(graphC="transparent",reducedOutput=TRUE,reportHTML=TRUE,autoShow=FALSE,
                 height=400,fullGraphSize=0.5)
        statusStore<-list(lastOutput="System",
                          showSampleType="Variables",
                          showInferParam="Basic",
                          showInferDimension="1D",
                          showMultipleParam="Basic",
                          showMultipleDimension="1D",
                          whichShowMultipleOut="all",
                          showExploreParam="Basic",
                          showExploreDimension="1D",
                          whichShowExploreOut="all",
                          exploreMode="Design",
                          basicHelpWhich=1,
                          demoHelpWhich=c(0,0),
                          simHelpWhich=0,
                          openJamovi=0,
                          planMode=NULL,
                          nrowTableLM=1,
                          nrowTableSEM=1
        )
        braw.env$statusStore<<-statusStore
        # braw.env$counter<<-0
      }
    },
    
    .run = function() {
      # debug information

      # changing the mode selector triggers a call here which we can ignore
      if (!is.null(braw.env$statusStore$planMode) && braw.env$statusStore$planMode!=self$options$planOptions) {
        braw.env$statusStore$planMode<<-self$options$planOptions
        return()
      }
      
      doDemos<-4 # set to 5 to include path demos
      systemAsHTML<-TRUE
      nestedHelp<-TRUE
      joinSystem<-TRUE
      
      if (joinSystem) self$results$SystemHTML$setVisible(FALSE)
      else self$results$SystemHTML$setVisible(TRUE)
      
      statusStore<-braw.env$statusStore
      if (self$options$showHTML) {
        if (self$results$simGraph$visible) self$results$simGraph$setVisible(FALSE)
        if (self$results$simReport$visible) self$results$simReport$setVisible(FALSE)
        if (!self$results$simGraphHTML$visible) self$results$simGraphHTML$setVisible(TRUE)
      } else {
        if (self$results$simGraphHTML$visible) self$results$simGraphHTML$setVisible(FALSE)
        if (!self$results$simGraph$visible) self$results$simGraph$setVisible(TRUE)
        if (!self$results$simReport$visible) self$results$simReport$setVisible(TRUE)
      }

      old_simHelpWhich<-statusStore$simHelpWhich
      if (self$options$simHelp) statusStore$simHelpWhich<-1
      if (self$options$simPlanHelp) statusStore$simHelpWhich<-2
      if (self$options$simSingleHelp) statusStore$simHelpWhich<-3
      if (self$options$simMultipleHelp) statusStore$simHelpWhich<-4
      if (self$options$simExploreHelp) statusStore$simHelpWhich<-5
      if (statusStore$simHelpWhich==old_simHelpWhich) statusStore$simHelpWhich<-0

      old_demoHelpWhich<-statusStore$demoHelpWhich
      if (self$options$demosHelp) statusStore$demoHelpWhich<-c(1,0)
      if (self$options$demo1Help) statusStore$demoHelpWhich<-c(2,1)
      if (self$options$demo2Help) statusStore$demoHelpWhich<-c(3,1)
      if (self$options$demo3Help) statusStore$demoHelpWhich<-c(4,1)
      if (self$options$doProject1AhBtn) statusStore$demoHelpWhich<-c(2,2)
      if (self$options$doProject1BhBtn) statusStore$demoHelpWhich<-c(2,3)
      if (self$options$doProject1ChBtn) statusStore$demoHelpWhich<-c(2,4)
      if (self$options$doProject2AhBtn) statusStore$demoHelpWhich<-c(3,2)
      if (self$options$doProject2BhBtn) statusStore$demoHelpWhich<-c(3,3)
      if (self$options$doProject2ChBtn) statusStore$demoHelpWhich<-c(3,4)
      if (self$options$doProject3AhBtn) statusStore$demoHelpWhich<-c(4,2)
      if (self$options$doProject3BhBtn) statusStore$demoHelpWhich<-c(4,3)
      if (self$options$doProject3ChBtn) statusStore$demoHelpWhich<-c(4,4)
      if (doDemos>4 && self$options$doProject4AhBtn) statusStore$demoHelpWhich<-c(5,2)
      if (all(statusStore$demoHelpWhich==old_demoHelpWhich))statusStore$demoHelpWhich<-c(0,0)
      
      statusStore$openJamovi<-0
      if (self$options$doProject1A2Btn || self$options$doProject1A3Btn) statusStore$openJamovi<-1
      
      # get some display parameters for later
      # single sample
      makeSampleNow<-self$options$makeSampleBtn
      showSampleType<-self$options$showSampleType
      showInferParam<-paste0(self$options$inferVar1,";",self$options$inferVar2)
      showInferDimension<-self$options$showInferDimension
      # self$results$debug$setContent(showInferParam)
      # self$results$debug$setVisible(TRUE)
      
      # multiple samples
      makeMultipleNow<-self$options$makeMultipleBtn
      showMultipleParam<-self$options$showMultipleParam
      if (is.element(showMultipleParam,c("Basic","Custom"))) {
        showMultipleParam<-paste0(self$options$inferVar1,";",self$options$inferVar2)
      } 
      showMultipleDimension<-self$options$showMultipleDimension
      whichShowMultipleOut<-self$options$whichShowMultiple
      
      # explore
      makeExploreNow<-self$options$makeExploreBtn
      showExploreParam<-self$options$showExploreParam
      if (is.element(showExploreParam,c("Single"))) {
        showExploreParam<-self$options$inferVar1
      } 
      if (is.element(showExploreParam,c("Basic","Custom"))) {
        showExploreParam<-paste0(self$options$inferVar1,";",self$options$inferVar2)
      } 
      showExploreDimension<-self$options$showExploreDimension
      whichShowExploreOut<-self$options$whichShowMultiple
      
      # make all the standard things we need
      # the hypothesis has IV, DV and effect
      oldH<-braw.def$hypothesis
      #
      DV<-makeVariable(self$options$DVname,self$options$DVtype,
                       mu=self$options$DVmu,sd=self$options$DVsd,skew=self$options$DVskew,kurtosis=self$options$DVkurt,
                       ncats=self$options$DVncats,cases=self$options$DVcases,proportions=self$options$DVprops,
                       nlevs=self$options$DVnlevs,iqr=self$options$DViqr)
      IV<-makeVariable(self$options$IVname,self$options$IVtype,
                       mu=self$options$IVmu,sd=self$options$IVsd,skew=self$options$IVskew,kurtosis=self$options$IVkurt,
                       ncats=self$options$IVncats,cases=self$options$IVcases,proportions=self$options$IVprops,
                       nlevs=self$options$IVnlevs,iqr=self$options$IViqr)
      if (self$options$presetIV2!="none") {
        IV2<-makeVariable(self$options$IV2name,self$options$IV2type,
                          mu=self$options$IV2mu,sd=self$options$IV2sd,skew=self$options$IV2skew,kurtosis=self$options$IV2kurt,
                          ncats=self$options$IV2ncats,cases=self$options$IV2cases,proportions=self$options$IV2props,
                          nlevs=self$options$IV2nlevs,iqr=self$options$IV2iqr)
      } else {
        IV2<-NULL
      }
      
      effect<-makeEffect(rIV=self$options$EffectSize1,
                         rIV2=self$options$EffectSize2,
                         rIVIV2=self$options$EffectSize3,
                         rIVIV2DV=self$options$EffectSize12,
                         rSD=self$options$rSD,
                         Heteroscedasticity=self$options$Heteroscedasticity,
                         ResidDistr=self$options$Residuals,
                         world=makeWorld(worldOn=self$options$WorldOn,
                                         populationPDF=self$options$WorldPDF,
                                         populationRZ = self$options$WorldRZ,
                                         populationPDFk = self$options$WorldLambda,
                                         populationNullp = self$options$WorldNullP
                         )
      )
      
      hypothesis<-makeHypothesis(IV,IV2,DV,effect,layout=self$options$EffectConfig)
      changedH<- !identical(oldH,hypothesis)

      # set up the design variable
      oldD<-braw.def$design
      #
      design<-makeDesign(sN=self$options$SampleSize,
                         sNRand=self$options$SampleSpreadOn,sNRandK=self$options$SampleGamma,
                         sMethod=makeSampling(self$options$SampleMethod),
                         sMethodSeverity=self$options$PoorSamplingAmount,
                                sIV1Use=self$options$SampleUsage1,
                                sIV2Use=self$options$SampleUsage2,
                                sDependence=self$options$Dependence,
                                sOutliers=self$options$Outliers,
                                sNonResponse=self$options$NonResponse,
                                sRangeOn=self$options$LimitRange=="yes", 
                                sIVRange=(c(self$options$RangeMin,self$options$RangeMax)-hypothesis$IV$mu)/hypothesis$IV$sd, 
                                sCheating=self$options$Cheating,sCheatingLimit="Budget",sCheatingBudget=self$options$CheatingBudget,
                         Replication=makeReplication(On=self$options$ReplicationOn,
                                                     Power=self$options$ReplicationPower,
                                                     Repeats=self$options$ReplicationAttempts,
                                                     Keep=self$options$ReplicationDecision,
                                                     forceSigOriginal=self$options$ReplicationSigOriginal=="yes",forceSign=(self$options$ReplicationSign=="yes"),
                                                     RepAlpha=self$options$ReplicationAlpha,
                                                     PowerPrior=self$options$ReplicationPrior
                                                     )
                         )
      if (design$sNRand) design$sN<-self$options$SampleSizeM
      changedD<- !identical(oldD,design)
      
      # set up the evidence variable
      oldE<-braw.def$evidence
      #
      evidence<-makeEvidence(rInteractionOn=self$options$interaction=="yes",
                             ssqType=self$options$ssq,sigOnly=FALSE,
                             Welch=self$options$equalVar=="no",
                             Transform=self$options$Transform,
                             doSEM=self$options$doSEM,useAIC=self$options$useAIC
      )
      changedE<- !identical(oldE,evidence)
      braw.env$alphaSig<<-self$options$alphaSig
      
      # set up the explore variable
      oldX<-braw.def$explore
      #
      switch(self$options$exploreMode,
             "hypothesisExplore"={
               typeExplore<-self$options$hypothesisExploreList
               minV<-as.numeric(self$options$exploreMinValH)
               maxV<-as.numeric(self$options$exploreMaxValH)
               xlog<-self$options$exploreXLogH
               exploreNPoints<-as.numeric(self$options$exploreNPointsH)
             },
             "designExplore"={
               typeExplore<-self$options$designExploreList
               minV<-as.numeric(self$options$exploreMinValD)
               maxV<-as.numeric(self$options$exploreMaxValD)
               xlog<-self$options$exploreXLogD
               exploreNPoints<-as.numeric(self$options$exploreNPointsD)
             },
             "analysisExplore"={
               typeExplore<-self$options$analysisExploreList
               minV<-as.numeric(self$options$exploreMinValA)
               maxV<-as.numeric(self$options$exploreMaxValA)
               xlog<-self$options$exploreXLogA
               exploreNPoints<-as.numeric(self$options$exploreNPointsA)
             },
             "moreExplore"={
               typeExplore<-self$options$moreExploreList
               minV<-as.numeric(self$options$exploreMinValM)
               maxV<-as.numeric(self$options$exploreMaxValM)
               xlog<-self$options$exploreXLogM
               exploreNPoints<-as.numeric(self$options$exploreNPointsM)
             }
      )
      explore<-makeExplore(exploreType=typeExplore,
                           exploreNPoints=exploreNPoints,
                           minVal=minV,maxVal=maxV,
                           xlog=xlog)
      changedX<- !identical(oldX,explore)

      # set up the metaAnalysis variable
      oldM<-braw.def$metaAnalysis
      #
      metaAnalysis<-makeMetaAnalysis(nstudies=self$options$MetaAnalysisNStudies,
                                 analysisType=self$options$MetaAnalysisType,
                                 modelPDF="All",
                                 sigOnlySource=self$options$MetaAnalysisStudiesSig=="sigOnly",
                                 includeNulls=self$options$MetaAnalysisNulls=="yes",
                                 includeBias=self$options$MetaAnalysisBias=="yes")
      changedM<- !identical(oldM,metaAnalysis)
      
      # store the option variables inside the braw package
      braw.def$hypothesis<<-hypothesis
      braw.def$design<<-design
      braw.def$evidence<<-evidence
      braw.def$explore<<-explore
      braw.def$metaAnalysis<<-metaAnalysis

      # we pressed the "show" hypothesis button
      openSystem<-0
      if (self$options$showHypothesisBtn) openSystem<-1
      # now deal with a request for help/instructions
      # after we have set up the hypothesis
      indent<-0
      titleWidth<-135
      if (nestedHelp) {
        indent<-50
        titleWidth<-0
      }
      if (self$options$brawHelp) {
        basicHelp<-brawbasicHelp(open=statusStore$basicHelpWhich,indent,titleWidth)
        simHelp<-brawsimHelp(open=statusStore$simHelpWhich,indent,titleWidth)
        jamoviHelp<-brawJamoviHelp(open=statusStore$openJamovi,indent,titleWidth,
                                   hypothesis,design)

        demoHelp<-brawDemosHelp(statusStore$demoHelpWhich,indent,titleWidth,doDemos)

        if (nestedHelp) {
          open0<-max(0,statusStore$basicHelpWhich,any(statusStore$demoHelpWhich>0)*2,(statusStore$simHelpWhich>0)*3,statusStore$openJamovi*4)
          help<-generate_tab(
            title="Help:",
            plainTabs=TRUE,
            titleWidth=50,
            tabs=c("Basic","Demos","Simulations","Jamovi","Key"),
            tabContents=c(basicHelp,demoHelp,simHelp,jamoviHelp,BrawInstructions("Key")),
            open=open0
          )
        } else help<-paste0(basicHelp,demoHelp,simHelp,jamoviHelp) 
        statusStore$basicHelpWhich<-0 # so closed next time round
      } else help<-''
      
      if (systemAsHTML) {
        if (changedD) openSystem<-2
        if (changedH) openSystem<-1
        if (changedE) openSystem<-1
        assign("graphHTML",TRUE,braw.env)
        svgBox(200)
        h<-showHypothesis()
        e<-showPrediction()
        svgBox(180)
        d<-joinHTML(showDesign(),reportDesign())
        systemHTML<-generate_tab(
          titleWidth=50,
          title="Plan:",
          plain=(nchar(help)>0),
          # topMargin=15,
          # titleWidth=30,
          tabs=c("Hypothesis","Design","Expected"),
          tabContents = c(
            h,
            d,
            e
          ),
          open=openSystem
        )
        assign("graphHTML",self$options$showHTML,braw.env)
        svgBox(400)
        if (joinSystem) {
          help<-paste0(help,'<div><br></div>',systemHTML)
        } else {
          self$results$SystemHTML$setContent(systemHTML)
        }
      }

      if (nchar(help)>0) {
        self$results$BrawStatsInstructions$setContent(help)
        if (!self$results$BrawStatsInstructions$visible) self$results$BrawStatsInstructions$setVisible(TRUE)
      } else {
        if (self$results$BrawStatsInstructions$visible)  self$results$BrawStatsInstructions$setVisible(FALSE)
      }
      
      # unless we have done something, we will make the same output as last time
      outputNow<-NULL
      
      # are we asking for a different display of the current explore?
      if (!is.null(braw.res$explore) && statusStore$lastOutput=="Explore") {
        if (showExploreParam != statusStore$showExploreParam ||
            showExploreDimension != statusStore$showInferDimension ||
            whichShowExploreOut != statusStore$whichShowExploreOut)
          outputNow<-"Explore"
      }
      # or multiple?
      if (!is.null(braw.res$multiple) && statusStore$lastOutput=="Multiple") {
        if (showMultipleParam != statusStore$showMultipleParam ||
            showMultipleDimension != statusStore$showInferDimension || 
            whichShowMultipleOut != statusStore$whichShowMultipleOut)
          outputNow<-"Multiple"
      }
      # or result?
      if (!is.null(braw.res$result) && is.element(statusStore$lastOutput,c("Compact","Sample","Describe","Infer","Variables","Likelihood"))) {
        if (showInferParam != statusStore$showInferParam ||
            showInferDimension != statusStore$showInferDimension)
          outputNow<-"Infer"
        if (showSampleType != statusStore$showSampleType)
          outputNow<-showSampleType
      }
      
      # are any of the existing stored results now invalid?
      if (changedH || changedD) {
        braw.res$result<<-NULL
        braw.res$multiple<<-NULL
        braw.res$explore<<-NULL
        braw.res$metaAnalysis<<-NULL
        outputNow<-"System"
      }
      if (changedE) {
        braw.res$multiple<<-NULL
        braw.res$explore<<-NULL
        braw.res$metaAnalysis<<-NULL
        if (!is.null(braw.res$result) && is.element(statusStore$lastOutput,c("Compact","Sample","Describe","Infer","Variables","Likelihood"))) {
          braw.res$result<<-doAnalysis(sample=braw.res$result)
          outputNow<-showSampleType
        }
      }
      if (changedX) {
        braw.res$explore<<-NULL
      }
      if (changedM) {
        braw.res$metaAnalysis<<-NULL
      }

      if (systemAsHTML && !is.null(outputNow) && outputNow=="System") outputNow<-NULL

      # now we start doing things
      # did we ask for a new sample?
      if (makeSampleNow) {
        # make a sample - either striaght sample or single metaAnalysis
        if (self$options$MetaAnalysisOn) {
          # do we need to do this, or are we just returning to the existing one?
          if (is.null(braw.res$metaSingle) || statusStore$lastOutput=="MetaSingle") 
            doMetaAnalysis(1,NULL)
          outputNow<-"MetaSingle"
        } else {
          # do we need to do this, or are we just returning to the existing one?
          if (is.null(braw.res$result) || statusStore$lastOutput==showSampleType) 
            doResult()
          outputNow<-showSampleType
        }
      }

      # did we ask for new multiples?
      if (makeMultipleNow) {
        numberSamples<-self$options$numberSamples
        if (self$options$MetaAnalysisOn) {
          if (statusStore$lastOutput=="MetaMultiple") 
            metaMultipleResult<-doMetaAnalysis(numberSamples,braw.res$metaAnalysis)
          outputNow<-"MetaMultiple"
        } else {
          # do we need to do this, or are we just returning to the existing one?
          if (is.null(braw.res$multiple) || statusStore$lastOutput=="Multiple") 
          doMultiple(nsims=numberSamples,multipleResult=braw.res$multiple)
          outputNow<-"Multiple"
        }
      }
      
      # did we ask for new explore?
      if (makeExploreNow) {
        numberExplores<-self$options$numberExplores
        # do we need to do this, or are we just returning to the existing one?
        if (is.null(braw.res$explore) || statusStore$lastOutput=="Explore") 
          exploreResult<-doExplore(nsims=numberExplores,exploreResult=braw.res$explore)
        outputNow<-"Explore"
      }

      # self$results$debug$setVisible(TRUE)
      # self$results$debug$setContent(outRep)
      
      # what are we showing?
      # main results graphs/reports
      if (!is.null(outputNow))  {

        if (outputNow=="Likelihood") {
          possible<-makePossible(UsePrior=self$options$likelihoodUsePrior,
                                 prior=makeWorld(worldOn=TRUE,
                                                 populationPDF=self$options$priorPDF,
                                                 populationRZ=self$options$priorRZ,
                                                 populationPDFk=self$options$priorLambda,
                                                 populationNullp=self$options$priorNullP)
          )
          possibleResult<-doPossible(possible)
          likelihoodCutaway<-(self$options$likelihoodCutaway=="cutaway")
        }
        
          if (self$options$showHTML) {
            svgBox(height=350,aspect=1.5)
            # if (outputNow=="System") svgBox(aspect=1.5)
            # if (is.element(outputNow,c("Compact","Describe"))) svgBox(height=300)
            # if (outputNow=="Likelihood") svgBox(aspect=1.5)
            # if (is.element(outputNow,c("Infer","Multiple")) && is.element(showMultipleParam,c("Basic","Custom")) && showMultipleDimension=="2D") svgBox(height=250)
            # if (is.element(outputNow,c("Infer","Multiple")) && showMultipleDimension=="1D") svgBox(aspect=1.6)
            # if (outputNow=="Explore" && is.element(showExploreParam,c("Basic","Custom")) && showExploreDimension=="2D") svgBox(height=250)
            # if (outputNow=="Explore" && is.element(showExploreParam,c("Basic","Custom")) && showExploreDimension=="1D") svgBox(aspect=1.6)
            # if (outputNow=="Explore" && showExploreParam=="NHST") svgBox(aspect=1.6)
            assign("graphHTML",TRUE,braw.env)
            graphSystem<-paste0(showSystem(),reportSystem())
            if (!is.null(braw.res$result))
              switch(self$options$showSampleType,
                     "Compact"= graphSingle<-paste0(showDescription(),reportInference()),
                     "Variables" = graphSingle<-paste0(showMarginals(style="all"),reportSample()),
                     "Sample"= graphSingle<-paste0(showSample(),reportSample()),
                     "Describe"= graphSingle<-paste0(showDescription(),reportDescription()),
                     "Infer"= graphSingle<-paste0(showInference(showType=showInferParam,dimension=showInferDimension),reportInference()),
                     "Likelihood"=graphSingle<-paste0(showPossible(showType=self$options$likelihoodType,cutaway=likelihoodCutaway),reportLikelihood())
                     )
            else graphSingle<-nullPlot()
            if (!is.null(braw.res$multiple))
              graphMultiple<-paste0(showMultiple(showType=showMultipleParam,dimension=showMultipleDimension,effectType=whichShowMultipleOut),
                                    reportMultiple(showType=showMultipleParam,reportStats=self$options$reportInferStats)
              )
            else graphMultiple<-nullPlot()
            if (!is.null(braw.res$explore))
              graphExplore<-paste0(showExplore(showType=showExploreParam,dimension=showExploreDimension,effectType=whichShowExploreOut),
                                   reportExplore(showType=showExploreParam,effectType=whichShowExploreOut,reportStats=self$options$reportInferStats)
                                   )
            else graphExplore<-nullPlot()
            switch(outputNow,
                   "System"= open<-0,
                   "Compact"= open<-1,
                   "Variables"= open<-1,
                   "Sample"= open<-1,
                   "Describe"= open<-1,
                   "Infer"= open<-1,
                   "Likelihood"=open<-1,
                   "Multiple"= open<-2,
                   "Explore"= open<-3,
                   {open<-0}
            )
      if (systemAsHTML) {
        brawResults<-generate_tab(
          title="Results:",
          titleWidth=50,
          tabs=c("Single Sample","Multiple Samples","Explore"),
          tabContents = c(
            graphSingle,
            graphMultiple,
            graphExplore
          ),
          open=open
        )
      } else {
        brawResults<-generate_tab(
          title="BrawStats Results:",
          tabs=c("Plan","Single Sample","Multiple Samples","Explore"),
          tabContents = c(
            graphSystem,
            graphSingle,
            graphMultiple,
            graphExplore
          ),
          open=open
        )
      }
      self$results$simGraphHTML$setContent(brawResults)
            # switch(outputNow,
            #        "System"= self$results$graphHTML$setContent(showSystem()),
            #        "Compact"= self$results$graphHTML$setContent(showDescription()),
            #        "Sample"= self$results$graphHTML$setContent(showMarginals(style=self$options$showSample)),
            #        "Describe"= self$results$graphHTML$setContent(showDescription()),
            #        "Infer"= self$results$graphHTML$setContent(showInference(showType=showInferParam,dimension=showInferDimension)),
            #        "Likelihood"=self$results$graphHTML$setContent(showPossible(showType=self$options$likelihoodType,cutaway=likelihoodCutaway)),
            #        "Multiple"= self$results$graphHTML$setContent(showMultiple(showType=showMultipleParam,dimension=showMultipleDimension,effectType=whichShowMultipleOut)),
            #        "Explore"= self$results$graphHTML$setContent(showExplore(showType=showExploreParam,dimension=showExploreDimension,effectType=whichShowExploreOut)),
            #        self$results$graphHTML$setContent(NULL)
            # )
            # svgBox(300)
          } else {

          switch(outputNow,
                 "System"= self$results$simGraph$setState(outputNow),
                 "Compact"= self$results$simGraph$setState("Describe"),
                 "Variables"= self$results$simGraph$setState(c(outputNow,"all")),
                 "Sample"= self$results$simGraph$setState(outputNow),
                 "Describe"= self$results$simGraph$setState(outputNow),
                 "Infer"= self$results$simGraph$setState(c(outputNow,showInferParam,showInferDimension)),
                 "Likelihood"= self$results$simGraph$setState(c(outputNow,self$options$likelihoodType,likelihoodCutaway)),
                 "Multiple"= self$results$simGraph$setState(c(outputNow,showMultipleParam,showMultipleDimension,whichShowMultipleOut)),
                 "Explore"= self$results$simGraph$setState(c(outputNow,showExploreParam,showExploreDimension,whichShowExploreOut)),
                 "MetaSingle"  =self$results$simGraph$setState(outputNow),
                 "MetaMultiple"  =self$results$simGraph$setState(outputNow),
                 self$results$simGraph$setState(outputNow)
          )
            
            switch(outputNow,
                   "System"= self$results$simReport$setContent(reportPlot(NULL)),
                   "Compact"= self$results$simReport$setContent(reportInference()),
                   "Variables"= self$results$simReport$setContent(reportSample()),
                   "Sample"= self$results$simReport$setContent(reportSample()),
                   "Describe"= self$results$simReport$setContent(reportDescription()),
                   "Infer"= self$results$simReport$setContent(reportInference()),
                   "Likelihood"=self$results$simReport$setContent(reportLikelihood()),
                   "Multiple"= self$results$simReport$setContent(reportMultiple(showType=showMultipleParam,reportStats=self$options$reportInferStats)),
                   "Explore"= self$results$simReport$setContent(reportExplore(showType=showExploreParam,reportStats=self$options$reportInferStats)),
                   "MetaSingle"  =self$results$simReport$setContent(reportMetaAnalysis()),
                   "MetaMultiple"  =self$results$simReport$setContent(reportMetaAnalysis()),
                   self$results$simReport$setContent(reportPlot(NULL))
            )
          }
      } else {
        if (!self$options$showHTML) {
          # self$results$simReport$setContent(reportPlot(NULL))
          # self$results$simGraph$setVisible(FALSE)
        }
      }
      
      # save everything for the next round      
      statusStore$showSampleType<-showSampleType
      statusStore$showInferParam<-showInferParam
      statusStore$showInferDimension<-showInferDimension
      statusStore$showMultipleParam<-showMultipleParam
      statusStore$showMultipleDimension<-showMultipleDimension
      statusStore$whichShowMultipleOut<-whichShowMultipleOut
      statusStore$showExploreParam<-showExploreParam
      statusStore$showExploreDimension<-showExploreDimension
      statusStore$whichShowExploreOut<-whichShowExploreOut
      if (!is.null(outputNow))  statusStore$lastOutput<-outputNow
      braw.env$statusStore<<-statusStore
      
      # now we save any results to the Jamovi spreadsheet
      # single result first
      if (!is.null(braw.res$result)) {
        if (is.null(IV2)) {
          newVariables<-data.frame(braw.res$result$participant,braw.res$result$dv,braw.res$result$iv,braw.res$result$dv+NA)
          names(newVariables)<-c("ID",DV$name,IV$name,"-")
        } else {
          newVariables<-data.frame(braw.res$result$participant,braw.res$result$dv,braw.res$result$iv,braw.res$result$iv2)
          names(newVariables)<-c("ID",DV$name,IV$name,IV2$name)
        }

        keys<-1:length(newVariables)
        measureTypes<-sapply(newVariables,function(x) { if (is.character(x)) "Nominal" else "Continuous"})
        
        self$results$sendSample$set(keys=keys,titles=names(newVariables),
                                    descriptions=rep("simulated",length(newVariables)),
                                    measureTypes=measureTypes
        )
        self$results$sendSample$setValues(newVariables)
      }
      # then multiple result
      q<-NULL
      if (!is.null(outputNow) && outputNow=="MetaSingle") {
          q<-braw.res$metaSingle$result
      } 
      if (!is.null(outputNow) && outputNow=="Multiple") {
        q<-mergeMultiple(braw.res$multiple$result,braw.res$multiple$nullresult)
      }
      if (!is.null(q)) {
        newMultiple<-data.frame(q$rIV,q$nval+0.0,q$pIV)
        newMultiple<-newMultiple[!is.na(q$rIV),]
        names(newMultiple)<-c("rs","n","p")
        nvars<-ncol(newMultiple)
        
        keys<-1:nvars
        self$results$sendMultiple$set(keys=keys,titles=names(newMultiple),
                                      descriptions=rep("simulated",nvars),
                                      measureTypes=rep("Continuous",nvars)
        )
        self$results$sendMultiple$setValues(newMultiple)
      }
      # end of .run()
    },
    
    .plotSimGraph=function(image, ...) {
      outputGraph <- image$state[1]
      if (!is.null(outputGraph)) {
        assign("graphHTML",FALSE,braw.env)
        switch(outputGraph,
               "System"    =outputGraph<-showSystem(),
               "Hypothesis"=outputGraph<-showHypothesis(),
               "Design"    =outputGraph<-showDesign(),
               "Population"=outputGraph<-showPopulation(),
               "Prediction"=outputGraph<-showPrediction(),
               "Variables"    =outputGraph<-showMarginals(style=image$state[2]),
               "Sample"    =outputGraph<-showSample(),
               "Describe"  =outputGraph<-showDescription(),
               "Infer"     =outputGraph<-showInference(showType=image$state[2],dimension=image$state[3]),
               "Likelihood"=outputGraph<-showPossible(showType=image$state[2],cutaway=as.logical(image$state[3])),
               "Multiple"  =outputGraph<-showMultiple(showType=image$state[2],dimension=image$state[3],effectType=image$state[4]),
               "MetaSingle"  =outputGraph<-showMetaSingle(),
               "MetaMultiple"  =outputGraph<-showMetaMultiple(),
               "Explore"   =outputGraph<-showExplore(showType=image$state[2],dimension=image$state[3],effectType=image$state[4])
        )
        print(outputGraph)
        return(TRUE)
      } else {
        return(FALSE)
      }
    },
    
    .plotReport=function(image, ...) {
      outputGraph <- image$state[1]
      if (!is.null(outputGraph)) {
        switch(outputGraph,
               "System"    =outputGraph<-reportTerms(),
               "Prediction"    =outputGraph<-reportTerms(),
               "Variables"    =outputGraph<-reportSample(),
               "Sample"    =outputGraph<-reportSample(),
               "Describe"  =outputGraph<-reportDescription(),
               "Infer"     =outputGraph<-reportInference(),
               "MetaSingle"  =outputGraph<-reportMetaAnalysis(),
               "MetaMultiple"  =outputGraph<-reportMetaAnalysis(),
               "Multiple"  =outputGraph<-reportMultiple(showType=image$state[2]),
               "Explore"   =outputGraph<-reportExplore(showType=image$state[2])
        )
        print(outputGraph)
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  )
)
