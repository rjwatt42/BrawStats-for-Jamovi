const events =  {

    // in here is where the event functions go
    onChange_exploreMode: function(ui) {
      var newRange = {min:0.3,max:0.7,xlog:false,np:8};
      let mode = ui.exploreMode.value();
      let value="n";
      switch(mode) {
        case "hypothesisExplore":
          value = ui.hypothesisExploreList.value();
          var newRange = updateRange(value);
              ui.exploreMinValH.setValue(newRange.min);
              ui.exploreMaxValH.setValue(newRange.max);
              ui.exploreXLogH.setValue(newRange.xlog);
              ui.exploreNPointsH.setValue(newRange.np);
          break;
        case "designExplore":
          value = ui.designExploreList.value();
          var newRange = updateRange(value);
              ui.exploreMinValD.setValue(newRange.min);
              ui.exploreMaxValD.setValue(newRange.max);
              ui.exploreXLogD.setValue(newRange.xlog);
              ui.exploreNPointsD.setValue(newRange.np);
          break;
        case "analysisExplore":
          value = ui.analysisExploreList.value();
          var newRange = updateRange(value);
              ui.exploreMinValA.setValue(newRange.min);
              ui.exploreMaxValA.setValue(newRange.max);
              ui.exploreXLogA.setValue(newRange.xlog);
              ui.exploreNPointsA.setValue(newRange.np);
          break;
        case "moreExplore":
          value = ui.moreExploreList.value();
          var newRange = updateRange(value);
              ui.exploreMinValM.setValue(newRange.min);
              ui.exploreMaxValM.setValue(newRange.max);
              ui.exploreXLogM.setValue(newRange.xlog);
              ui.exploreNPointsM.setValue(newRange.np);
          break;
      };
      return;
  },

    onChange_hypothesisExploreList: function(ui) {
          let value = ui.hypothesisExploreList.value();
          let newRange = updateRange(value)
          ui.exploreMinValH.setValue(newRange.min);
          ui.exploreMaxValH.setValue(newRange.max);
          ui.exploreXLogH.setValue(newRange.xlog);
          ui.exploreNPointsH.setValue(newRange.np);
      return;
  },
    onChange_designExploreList: function(ui) {
          let value = ui.designExploreList.value();
          let newRange = updateRange(value)
          ui.exploreMinValD.setValue(newRange.min);
          ui.exploreMaxValD.setValue(newRange.max);
          ui.exploreXLogD.setValue(newRange.xlog);
          ui.exploreNPointsD.setValue(newRange.np);
      return;
  },
    onChange_analysisExploreList: function(ui) {
          let value = ui.analysisExploreList.value();
          let newRange = updateRange(value)
          ui.exploreMinValA.setValue(newRange.min);
          ui.exploreMaxValA.setValue(newRange.max);
          ui.exploreXLogA.setValue(newRange.xlog);
          ui.exploreNPointsA.setValue(newRange.np);
  },
    onChange_moreExploreList: function(ui) {
          let value = ui.moreExploreList.value();
          let newRange = updateRange(value)
          ui.exploreMinValM.setValue(newRange.min);
          ui.exploreMaxValM.setValue(newRange.max);
          ui.exploreXLogM.setValue(newRange.xlog);
          ui.exploreNPointsM.setValue(newRange.np);
  },

    onChange_presetDV: function(ui) {
      let presetDV = ui.presetDV.value();
      let variable = makeVar(presetDV);
          ui.DVname.setValue(variable.name);
          ui.DVtype.setValue(variable.type);
          ui.DVmu.setValue(variable.mu);
          ui.DVsd.setValue(variable.sd);
          ui.DVskew.setValue(variable.skew);
          ui.DVkurt.setValue(variable.kurt);
          ui.DVnlevs.setValue(variable.nlevs);
          ui.DViqr.setValue(variable.iqr)
          ui.DVncats.setValue(variable.ncats);
          ui.DVcases.setValue(variable.cases)
          ui.DVprops.setValue(variable.props);
    },

    onChange_presetIV: function(ui) {
      let presetIV = ui.presetIV.value();
      let oldmu = ui.IVmu.value();
      let oldsd = ui.IVsd.value();
      let variable = makeVar(presetIV);
          ui.IVname.setValue(variable.name);
          ui.IVtype.setValue(variable.type);
          ui.IVmu.setValue(variable.mu);
          ui.IVsd.setValue(variable.sd);
          ui.IVskew.setValue(variable.skew);
          ui.IVkurt.setValue(variable.kurt);
          ui.IVnlevs.setValue(variable.nlevs);
          ui.IViqr.setValue(variable.iqr);
          ui.IVncats.setValue(variable.ncats);
          ui.IVcases.setValue(variable.cases)
          ui.IVprops.setValue(variable.props);
          ui.RangeMin.setValue(variable.mu+variable.sd*(ui.RangeMin.value()-oldmu)/oldsd);
          ui.RangeMax.setValue(variable.mu+variable.sd*(ui.RangeMax.value()-oldmu)/oldsd);
    },

    onChange_presetIV2: function(ui) {
      let presetIV2 = ui.presetIV2.value();
      if (presetIV2!="none") {
      let variable = makeVar(presetIV2);
          ui.IV2name.setValue(variable.name);
          ui.IV2type.setValue(variable.type);
          ui.IV2mu.setValue(variable.mu);
          ui.IV2sd.setValue(variable.sd);
          ui.IV2skew.setValue(variable.skew);
          ui.IV2kurt.setValue(variable.kurt);
          ui.IV2nlevs.setValue(variable.nlevs);
          ui.IV2iqr.setValue(variable.iqr)
          ui.IV2ncats.setValue(variable.ncats);
          ui.IV2cases.setValue(variable.cases)
          ui.IV2props.setValue(variable.props);
      }
    },

    onChange_presetWorld: function(ui) {
      let presetH = ui.presetWorld.value();
      switch(presetH) {
        case "psych":
          ui.WorldOn.setValue(true);
          ui.WorldPDF.setValue("Exp");
          ui.WorldRZ.setValue("z");
          ui.WorldLambda.setValue(0.3);
          ui.WorldNullP.setValue(0.75);
          ui.SampleSpreadOn.setValue(true);
          ui.SampleSizeM.setValue(52);
          ui.SampleGamma.setValue(1.56);
          break;
        case "simple":
          ui.WorldOn.setValue(true);
          ui.WorldPDF.setValue("Single");
          ui.WorldRZ.setValue("r");
          ui.WorldLambda.setValue(0.3);
          ui.WorldNullP.setValue(0.5);
          ui.SampleSpreadOn.setValue(false);
          break;
      }
    },
    
    onChange_Single: function(ui) {
      ui.whichGraph.setValue("Single")
    },
    
    onChange_Project1sH: function(ui) {
        demo1SetUp(ui)
        ui.showHypothesisBtn.setValue(true)
    },
    
    onChange_Project1sD: function(ui) {
        demo1SetUp(ui)
        ui.showHypothesisBtn.setValue(true)
    },

    onChange_Project2sH: function(ui) {
        demo2SetUp(ui)
        ui.showHypothesisBtn.setValue(true)
    },
    
    onChange_Project2sD: function(ui) {
        demo2SetUp(ui)
        ui.showHypothesisBtn.setValue(true)
    },
    
    onChange_Project3sH: function(ui) {
        demo3SetUp(ui)
        ui.showHypothesisBtn.setValue(true)
    },
    
    onChange_Project3sD: function(ui) {
        demo3SetUp(ui)
        ui.showHypothesisBtn.setValue(true)
    },
    
    onChange_Project1As: function(ui) {
      let BtnOn1 = ui.doProject1AsBtn.value();
      let BtnOn2A = ui.doProject2AsBtn.value();
      let BtnOn2B = ui.doProject2BsBtn.value();
      let BtnOn3A = ui.doProject3AsBtn.value();
      let BtnOn3B = ui.doProject3BsBtn.value();
      let BtnOn3C = ui.doProject3CsBtn.value();
      if (BtnOn1==true) {
        demo1SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
      if (BtnOn2A==true) {
        demo2SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
      if (BtnOn2B==true) {
        demo2SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
      if (BtnOn3A==true) {
        demo3SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
      if (BtnOn3B==true) {
        demo3SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
      if (BtnOn3C==true) {
        demo3SetUp(ui)
        ui.makeExploreBtn.setValue(true)
      }
    },
    
    onChange_project1A1: function(ui) {
      let BtnOn = ui.doProject1A1Btn.value();
      if (BtnOn==true) {
        ui.doProject1sLstC.setValue("Compact")
        demo1SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
    },

    onChange_project1B1: function(ui) {
      let BtnOn = ui.doProject1B1Btn.value();
      if (BtnOn==true) {
        ui.doProject1sLstC.setValue("Variables")
        demo1SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project1B2: function(ui) {
      let BtnOn = ui.doProject1B2Btn.value();
      if (BtnOn==true) {
        ui.doProject1sLstC.setValue("Sample")
        demo1SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project1B3: function(ui) {
      let BtnOn = ui.doProject1B3Btn.value();
      if (BtnOn==true) {
        ui.doProject1sLstD.setValue("Convenience")
      }
    },

    onChange_project1C1: function(ui) {
      let BtnOn = ui.doProject1C1Btn.value();
      if (BtnOn==true) {
        ui.doProject1sLstC.setValue("Describe")
        demo1SetUp(ui)
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project1C2: function(ui) {
      let BtnOn = ui.doProject1C2Btn.value();
      if (BtnOn==true) {
        ui.doProject1sLstE.setValue(500)
      }
    },
    
    onChange_project1C3: function(ui) {
      let BtnOn = ui.doProject1C3Btn.value();
      if (BtnOn==true) {
        ui.doProject1sLstE.setValue(42)
        ui.doProject1sLstD.setValue("Convenience")
      }
    },
    
    onChange_project2As: function(ui) {
      let BtnOn = ui.doProject2AsBtn.value();
      if (BtnOn==true) {
        demo2SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project2A1: function(ui) {
      let BtnOn = ui.doProject2A1Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstC.setValue("Infer")
        demo2SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
    },
        
    onChange_project2A2: function(ui) {
      let BtnOn = ui.doProject2A2Btn.value();
      if (BtnOn==true) {
        demo2SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
    },
    
    onChange_project2A3: function(ui) {
      let BtnOn = ui.doProject2A3Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstD.setValue("Convenience")
      }
    },
    
    onChange_Project2B1: function(ui) {
      let BtnOn = ui.doProject2B1Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstG.setValue("eq")
        ui.doProject2sLstC.setValue("Infer")
        demo2SetUp(ui)
        ui.makeSampleBtn.setValue(true)
      }
    },
        
    onChange_Project2B2: function(ui) {
      let BtnOn = ui.doProject2B2Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstG.setValue("eq")
        ui.doProject2sLstF.setValue("Basic")
        demo2SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
    },
    
    onChange_Project2B3: function(ui) {
      let BtnOn = ui.doProject2B3Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstG.setValue("eq")
        ui.doProject2sLstE.setValue("n320")
      }
    },

    onChange_project2C1: function(ui) {
      let BtnOn = ui.doProject2C1Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstG.setValue("eq")
        ui.doProject2sLstF.setValue("NHST")
        demo2SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
    },
    
    onChange_project2C2: function(ui) {
      let BtnOn = ui.doProject2C2Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstB.setValue("r05")
      }
    },
    
    onChange_project2C3: function(ui) {
      let BtnOn = ui.doProject2C3Btn.value();
      if (BtnOn==true) {
        ui.doProject2sLstB.setValue("r03")
        ui.doProject2sLstD.setValue("Convenience")
        demo2SetUp(ui)
      }
    },
    
    onChange_project3A1: function(ui) {
      let BtnOn = ui.doProject3A1Btn.value();
      if (BtnOn==true) {
        ui.doProject3sLstF.setValue("p(sig)")
        demo3SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
    },

    onChange_project3A3: function(ui) {
      let BtnOn = ui.doProject3A3Btn.value();
      if (BtnOn==true) {
        ui.doProject3sLstH.setValue("n")
        ui.doProject3sLstI.setValue("p(sig)")
        demo3SetUp(ui)
        ui.makeExploreBtn.setValue(true)
      }
    },

    onChange_project3B1: function(ui) {
      let BtnOn = ui.doProject3B1Btn.value();
      if (BtnOn==true) {
        ui.doProject3sLstH.setValue("rIV")
        ui.doProject3sLstI.setValue("p(sig)")
        demo3SetUp(ui)
        ui.makeExploreBtn.setValue(true)
      }
    },

    onChange_project3C1: function(ui) {
      let BtnOn = ui.doProject3C1Btn.value();
      if (BtnOn==true) {
        ui.doProject3sLst.setValue("Categorical")
        ui.doProject3sLstF.setValue("p(sig)")
        demo3SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
    },

    onChange_project3C2: function(ui) {
      let BtnOn = ui.doProject3C1Btn.value();
      if (BtnOn==true) {
        ui.doProject3sLst.setValue("Categorical")
        ui.doProject3sLstJ.setValue("Within")
        ui.doProject3sLstF.setValue("p(sig)")
        demo3SetUp(ui)
        ui.makeMultipleBtn.setValue(true)
      }
    },

    onChange_project3C3: function(ui) {
      let BtnOn = ui.doProject3C1Btn.value();
      if (BtnOn==true) {
        ui.doProject3sLstH.setValue("rIV")
        ui.doProject3sLstI.setValue("p(sig)")
        demo3SetUp(ui)
        ui.makeExploreBtn.setValue(true)
      }
    },
    
    onChange_project4As: function(ui) {
      let BtnOn = ui.doProject4AsBtn.value();
      if (BtnOn==true) {
      let hyp = ui.doProject4AsLst.value();
      defaultSetUp(ui)
      ui.presetIV2.setValue("IV2");
      ui.interaction.setValue("yes");
      switch(hyp) {
        case "iii": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Interval"); 
          ui.IV2type.setValue("Interval"); 
        break;
        case "ici": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Categorical"); 
          ui.IV2type.setValue("Interval"); 
        break;
        case "iic": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Interval"); 
          ui.IV2type.setValue("Categorical"); 
        break;
        case "icc": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Categorical"); 
          ui.IV2type.setValue("Categorical"); 
        break;
      }
      ui.EffectSize1.setValue(0.0)
      ui.EffectSize2.setValue(0.0)
      ui.EffectSize3.setValue(0.0)
      ui.SampleSize.setValue(1000)
      ui.showHypothesisBtn.setValue(true)
      }
    },
    
    onChange_project4A1: function(ui) {
      let BtnOn = ui.doProject4A1Btn.value();
      if (BtnOn==true) {
      let rVal = ui.doProject4A1Lst.value();
      switch(rVal) {
        case "r00": ui.EffectSize1.setValue(0.0); break;
        case "r01": ui.EffectSize1.setValue(0.1); break;
        case "r02": ui.EffectSize1.setValue(0.2); break;
        case "r03": ui.EffectSize1.setValue(0.3); break;
        case "r05": ui.EffectSize1.setValue(0.5); break;
        case "r06": ui.EffectSize1.setValue(0.6); break;
      }
      ui.EffectSize2.setValue(0.0)
      ui.EffectSize12.setValue(0.0)
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project4A2: function(ui) {
      let BtnOn = ui.doProject4A2Btn.value();
      if (BtnOn==true) {
      let rVal = ui.doProject4A2Lst.value();
      switch(rVal) {
        case "r00": ui.EffectSize2.setValue(0.0); break;
        case "r01": ui.EffectSize2.setValue(0.1); break;
        case "r02": ui.EffectSize2.setValue(0.2); break;
        case "r03": ui.EffectSize2.setValue(0.3); break;
        case "r05": ui.EffectSize2.setValue(0.5); break;
        case "r06": ui.EffectSize2.setValue(0.6); break;
      }
      ui.EffectSize1.setValue(0.0)
      ui.EffectSize12.setValue(0.0)
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project4A3: function(ui) {
      let BtnOn = ui.doProject4A3Btn.value();
      if (BtnOn==true) {
      let rVal = ui.doProject4A3Lst.value();
      switch(rVal) {
        case "r00": ui.EffectSize12.setValue(0.0); break;
        case "r01": ui.EffectSize12.setValue(0.1); break;
        case "r02": ui.EffectSize12.setValue(0.2); break;
        case "r03": ui.EffectSize12.setValue(0.3); break;
        case "r05": ui.EffectSize12.setValue(0.5); break;
        case "r075": ui.EffectSize12.setValue(0.75); break;
      }
      ui.EffectSize1.setValue(0.0)
      ui.EffectSize2.setValue(0.0)
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project4A4: function(ui) {
      let BtnOn = ui.doProject4A4Btn.value();
      if (BtnOn==true) {
      let rVal1 = ui.doProject4A1Lst.value();
      switch(rVal1) {
        case "r00": ui.EffectSize1.setValue(0.0); break;
        case "r01": ui.EffectSize1.setValue(0.1); break;
        case "r02": ui.EffectSize1.setValue(0.2); break;
        case "r03": ui.EffectSize1.setValue(0.3); break;
        case "r05": ui.EffectSize1.setValue(0.5); break;
        case "r075": ui.EffectSize1.setValue(0.75); break;
      }
      let rVal2 = ui.doProject4A2Lst.value();
      switch(rVal2) {
        case "r00": ui.EffectSize2.setValue(0.0); break;
        case "r01": ui.EffectSize2.setValue(0.1); break;
        case "r02": ui.EffectSize2.setValue(0.2); break;
        case "r03": ui.EffectSize2.setValue(0.3); break;
        case "r05": ui.EffectSize2.setValue(0.5); break;
        case "r075": ui.EffectSize2.setValue(0.75); break;
      }
      let rVal12 = ui.doProject4A3Lst.value();
      switch(rVal12) {
        case "r00": ui.EffectSize12.setValue(0.0); break;
        case "r01": ui.EffectSize12.setValue(0.1); break;
        case "r02": ui.EffectSize12.setValue(0.2); break;
        case "r03": ui.EffectSize12.setValue(0.3); break;
        case "r05": ui.EffectSize12.setValue(0.5); break;
        case "r075": ui.EffectSize12.setValue(0.75); break;
      }
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },

    onChange_project4Bs: function(ui) {
      let BtnOn = ui.doProject4BsBtn.value();
      if (BtnOn==true) {
      let hyp = ui.doProject4BsLst.value();
      defaultSetUp(ui)
      ui.presetIV2.setValue("IV2");
      ui.EffectConfig.setValue("normal")
      ui.interaction.setValue("yes");
      switch(hyp) {
        case "iii": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Interval"); 
          ui.IV2type.setValue("Interval"); 
        break;
        case "ici": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Categorical"); 
          ui.IV2type.setValue("Interval"); 
        break;
        case "iic": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Interval"); 
          ui.IV2type.setValue("Categorical"); 
        break;
        case "icc": 
          ui.DVtype.setValue("Interval"); 
          ui.IVtype.setValue("Categorical"); 
          ui.IV2type.setValue("Categorical"); 
        break;
      }
      ui.EffectSize1.setValue(0.0)
      ui.EffectSize2.setValue(0.0)
      ui.EffectSize12.setValue(0.0)
      ui.SampleSize.setValue(1000)
      ui.showHypothesisBtn.setValue(true)
      }
    },
    
    onChange_project4B1: function(ui) {
      let BtnOn = ui.doProject4B1Btn.value();
      if (BtnOn==true) {
      let rVal = ui.doProject4B1Lst.value();
      switch(rVal) {
        case "r00": 
          ui.EffectSize1.setValue(0.0); 
          ui.EffectSize12.setValue(0.0); 
          break;
        case "r01": 
          ui.EffectSize1.setValue(0.1); 
          ui.EffectSize12.setValue(0.1); 
          break;
        case "r02": 
          ui.EffectSize1.setValue(0.2); 
          ui.EffectSize12.setValue(0.2); 
          break;
        case "r03": 
          ui.EffectSize1.setValue(0.3); 
          ui.EffectSize12.setValue(0.3); 
          break;
        case "r05": 
          ui.EffectSize1.setValue(0.5); 
          ui.EffectSize12.setValue(0.5); 
          break;
        case "r075": 
          ui.EffectSize1.setValue(0.75); 
          ui.EffectSize12.setValue(0.75); 
          break;
      }
      ui.EffectSize2.setValue(0.0)
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project4B2: function(ui) {
      let BtnOn = ui.doProject4B2Btn.value();
      if (BtnOn==true) {
      let rVal = ui.doProject4B2Lst.value();
      switch(rVal) {
        case "r00": 
          ui.EffectSize12.setValue(0.0); 
          break;
        case "r01": 
          ui.EffectSize12.setValue(0.1); 
          break;
        case "r02": 
          ui.EffectSize12.setValue(0.2); 
          break;
        case "r03": 
          ui.EffectSize12.setValue(0.3); 
          break;
        case "r05": 
          ui.EffectSize12.setValue(0.5); 
          break;
        case "r075": 
          ui.EffectSize12.setValue(0.75); 
          break;
      }
          ui.EffectSize1.setValue(0.0); 
      ui.EffectSize2.setValue(0.0)
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },

    onChange_Project4Cs: function(ui) {
      let BtnOn = ui.doProject4CsBtn.value();
      if (BtnOn==true) {
      defaultSetUp(ui)
      ui.presetDV.setValue("ExamGrade");
      ui.presetIV.setValue("Perfectionism");
      ui.presetIV2.setValue("Anxiety");
      ui.interaction.setValue("no");
      ui.EffectSize1.setValue(0.3)
      ui.EffectSize2.setValue(-0.5)
      ui.EffectSize3.setValue(0.6)
      ui.SampleSize.setValue(1000)
      ui.showHypothesisBtn.setValue(true)
      }
    },
    
    onChange_project4C1: function(ui) {
      let BtnOn = ui.doProject4C1Btn.value();
      if (BtnOn==true) {
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project4C2: function(ui) {
      let BtnOn = ui.doProject4C2Btn.value();
      if (BtnOn==true) {
      let rVal = ui.doProject4C2Lst.value();
      switch(rVal) {
        case "r00": 
          ui.EffectSize3.setValue(0.0); 
          break;
        case "r01": 
          ui.EffectSize3.setValue(0.1); 
          break;
        case "r02": 
          ui.EffectSize3.setValue(0.2); 
          break;
        case "r03": 
          ui.EffectSize3.setValue(0.3); 
          break;
        case "r05": 
          ui.EffectSize3.setValue(0.5); 
          break;
        case "r075": 
          ui.EffectSize3.setValue(0.75); 
          break;
      }
        ui.showSampleType.setValue("Describe")
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project4C3: function(ui) {
      let BtnOn = ui.doProject4C3Btn.value();
      if (BtnOn==true) {
        ui.SampleSize.setValue(42)
        ui.whichShowMultiple.setValue("all")
        ui.exploreMode.setValue("hypothesisExplore")
        ui.hypothesisExploreList.setValue("rIVIV2")
        ui.showExploreParam.setValue("Single")
      ui.makeExploreBtn.setValue(true)
      }
    },
    
    
    onChange_Project5As: function(ui) {
      let BtnOn = ui.doProject5AsBtn.value();
      if (BtnOn==true) {
      defaultSetUp(ui)
      ui.SampleSize.setValue(1000)
          ui.presetDV.setValue("ExamGrade"); 
          ui.presetIV.setValue("HoursSleep"); 
      ui.EffectSize1.setValue(0.0)
          ui.presetIV2.setValue("Anxiety"); 
      ui.EffectSize2.setValue(0.7)
      ui.EffectSize3.setValue(0.7)
      ui.interaction.setValue("no");
      ui.EffectConfig.setValue("path")
      ui.showHypothesisBtn.setValue(true)
      }
    },
    
    onChange_project5A1: function(ui) {
      let BtnOn = ui.doProject5A1Btn.value();
      if (BtnOn==true) {
          ui.EffectSize3.setValue(0.7); 
        ui.showSampleType.setValue("Infer")
      ui.makeSampleBtn.setValue(true)
      }
    },
    
    onChange_project5A2: function(ui) {
      let BtnOn = ui.doProject5A2Btn.value();
      if (BtnOn==true) {
      let rVal = ui.doProject5A2Lst.value();
      switch(rVal) {
        case "rn07": 
          ui.EffectSize3.setValue(-0.7); 
          break;
        case "rn03": 
          ui.EffectSize3.setValue(-0.3); 
          break;
        case "rn01": 
          ui.EffectSize3.setValue(-0.1); 
          break;
        case "r00": 
          ui.EffectSize3.setValue(0.0); 
          break;
        case "r01": 
          ui.EffectSize3.setValue(0.1); 
          break;
        case "r03": 
          ui.EffectSize3.setValue(0.3); 
          break;
        case "r07": 
          ui.EffectSize3.setValue(0.7); 
          break;
      }
        ui.showSampleType.setValue("Infer")
      ui.makeSampleBtn.setValue(true)
      }
    },

    onChange_project5A3: function(ui) {
      let BtnOn = ui.doProject5A3Btn.value();
      if (BtnOn==true) {
        ui.SampleSize.setValue(42)
        ui.whichShowMultiple.setValue("all")
        ui.exploreMode.setValue("hypothesisExplore")
        ui.hypothesisExploreList.setValue("rIVIV2")
        ui.showExploreParam.setValue("Single")
      ui.makeExploreBtn.setValue(true)
      }
    }

}

let defaultSetUp = function(ui) {
      ui.DVtype.setValue("Interval")
        ui.presetIV2.setValue("none")
        ui.IVtype.setValue("Interval")
      ui.EffectSize1.setValue(0.3)
        ui.EffectSize2.setValue(0)
        ui.EffectSize3.setValue(0)
        ui.EffectSize12.setValue(0)
        ui.WorldOn.setValue(false)
      ui.EffectConfig.setValue("normal")
      ui.SampleSpreadOn.setValue(false)
      ui.SampleSize.setValue(42)
        ui.SampleMethod.setValue("Random")
      ui.doSEM.setValue(false)
      ui.inferVar1.setValue("rs")
      ui.inferVar2.setValue("p")
      ui.showInferDimension.setValue("1D")
        ui.showMultipleParam.setValue("Basic")
        ui.showMultipleDimension.setValue("1D")
}

    
let demo1SetUp = function(ui) {
//        defaultSetUp(ui)    
        let variable1 = ui.doProject1sLst.value();
        variable1 = variable1.replace("?","")
        ui.presetIV.setValue(variable1)
        let variable2 = ui.doProject1sLstA.value();
        variable2 = variable2.replace("?","")
        ui.presetDV.setValue(variable2)
        let variable3 = ui.doProject1sLstB.value();
        if (variable3>0.95) {
          variable3 = 0.95
        }
        if (variable3< -0.95) {
          variable3 = -0.95
        }
        ui.EffectSize1.setValue(variable3);
        let variable4 = ui.doProject1sLstC.value();
        ui.showSampleType.setValue(variable4)
        let variable5 = ui.doProject1sLstD.value()
        ui.SampleMethod.setValue(variable5)
        let variable6 = ui.doProject1sLstE.value()
        ui.SampleSize.setValue(variable6);
    }
    
let demo2SetUp = function(ui) {
//        defaultSetUp(ui)    
        let variable1 = ui.doProject2sLst.value();
        variable1 = variable1.replace("?","")
        ui.presetIV.setValue(variable1)
        let variable2 = ui.doProject2sLstA.value();
        variable2 = variable2.replace("?","")
        ui.presetDV.setValue(variable2)
        let variable3 = ui.doProject2sLstB.value();
        if (variable3>0.95) {
          variable3 = 0.95
        }
        if (variable3< -0.95) {
          variable3 = -0.95
        }
        ui.EffectSize1.setValue(variable3);
        let variable4 = ui.doProject2sLstC.value();
        ui.showSampleType.setValue(variable4)
        let variable5 = ui.doProject2sLstD.value()
        ui.SampleMethod.setValue(variable5)
        let variable6 = ui.doProject2sLstE.value()
        ui.SampleSize.setValue(variable6);
       let variable7 = ui.doProject2sLstF.value()
        ui.showMultipleParam.setValue(variable7)
       let variable8 = ui.doProject2sLstG.value()
      switch(variable8) {
        case "off": 
          ui.WorldOn.setValue(false);
          ui.SampleSpreadOn.setValue(false);
          break;
        case "eq": 
          ui.WorldOn.setValue(true);
          ui.WorldPDF.setValue("Single");
          ui.WorldRZ.setValue("r");
          ui.WorldLambda.setValue(0.3);
          ui.WorldNullP.setValue(0.5);
          ui.SampleSpreadOn.setValue(false);
          break;
      }

    }

let demo3SetUp = function(ui) {
//        defaultSetUp(ui)    
        let variable1 = ui.doProject3sLst.value();
        switch (variable1) {
          case "Interval":
            ui.presetIV.setValue("IV")
            ui.IVtype.setValue("Interval")
            break;
          case "Categorical":
            ui.presetIV.setValue("IV")
            ui.IVtype.setValue("Categorical")
            break;
        }
        let variable2 = ui.doProject3sLstA.value();
        switch (variable2) {
          case "Interval":
            ui.presetDV.setValue("DV")
            ui.DVtype.setValue("Interval")
            break;
          case "Categorical":
            ui.presetDV.setValue("DV")
            ui.DVtype.setValue("Categorical")
            break;
        }
        let variable3 = ui.doProject3sLstB.value();
        if (variable3>0.95) {
          variable3 = 0.95
        }
        if (variable3< -0.95) {
          variable3 = -0.95
        }
        ui.EffectSize1.setValue(variable3);
        let variable4 = ui.doProject3sLstC.value();
        ui.showSampleType.setValue(variable4)
        let variable5 = ui.doProject3sLstD.value()
        ui.SampleMethod.setValue(variable5)
        let variable6 = ui.doProject3sLstE.value()
        ui.SampleSize.setValue(variable6);
       let variable7 = ui.doProject3sLstF.value()
        ui.showMultipleParam.setValue(variable7)
       let variable8 = ui.doProject3sLstG.value()
      switch(variable8) {
        case "off": 
          ui.WorldOn.setValue(false);
          ui.SampleSpreadOn.setValue(false);
          break;
        case "eq": 
          ui.WorldOn.setValue(true);
          ui.WorldPDF.setValue("Single");
          ui.WorldRZ.setValue("r");
          ui.WorldLambda.setValue(0.3);
          ui.WorldNullP.setValue(0.5);
          ui.SampleSpreadOn.setValue(false);
          break;
      }
       let variable9 = ui.doProject3sLstH.value()
       switch (variable9) {
         case "n":
          ui.exploreMode.setValue("designExplore")
          ui.designExploreList.setValue("n")
          break;
         case "rIV":
          ui.exploreMode.setValue("hypothesisExplore")
          ui.hypothesisExploreList.setValue("rIV")
          break;
       }
       let variable10 = ui.doProject3sLstI.value()
      ui.showExploreParam.setValue(variable10)
      ui.numberExplores.setValue(100)
      let variable11 = ui.doProject3sLstJ.value()
      ui.SampleUsage1.setValue(variable11)
    }

let makeVar = function(name) {
  switch (name) {
    case "DV":
      var variable={name:"DV",type:"Interval",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "IV":
      var variable={name:"IV",type:"Interval",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "IV2":
      var variable={name:"IV2",type:"Interval",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "IQ":
      var variable={name:"IQ",type:"Interval",mu:100,sd:15,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "Diligence":
      var variable={name:"Diligence",type:"Interval",mu:0,sd:2,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "Perfectionism":
      var variable={name:"Perfectionism",type:"Interval",mu:0,sd:2,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "Anxiety":
      var variable={name:"Anxiety",type:"Interval",mu:5,sd:2,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "Happiness":
      var variable={name:"Happiness",type:"Interval",mu:50,sd:12,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "SelfConfidence":
      var variable={name:"SelfConfidence",type:"Interval",mu:50,sd:12,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "HoursSleep":
      var variable={name:"HoursSleep",type:"Interval",mu:7,sd:1,skew:-0.7,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "ExamGrade":
      var variable={name:"ExamGrade",type:"Interval",mu:55,sd:10,skew:-0.6,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "ExamPass":
      var variable={name:"ExamPass?",type:"Categorical",mu:55,sd:10,skew:-0.6,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"no,yes",props:"1,3"};
      break;
    case "RiskTaking":
      var variable={name:"RiskTaking",type:"Interval",mu:30,sd:6,skew:0.5,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "Interesting":
      var variable={name:"Interesting",type:"Interval",mu:10,sd:2,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"C1,C2",props:"1,1"};
      break;
    case "Musician":
      var variable={name:"Musician?",type:"Categorical",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"no,yes",props:"1,1"};
      break;
    case "RiskTaker":
      var variable={name:"RiskTaker?",type:"Categorical",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"no,yes",props:"1,1"};
      break;
    case "Smoker":
      var variable={name:"Smoker?",type:"Categorical",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"no,yes",props:"2,1"};
      break;
    case "StudySubject":
      var variable={name:"StudySubject",type:"Categorical",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:3,cases:"psych,phil,sports",props:"1.5,1,2"};
      break;
    case "BirthOrder":
      var variable={name:"BirthOrder",type:"Categorical",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:4,cases:"first,middle,last,only",props:"1,0.4,0.6,0.2"};
      break;
    case "Treatment":
      var variable={name:"Treatment?",type:"Categorical",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"no,yes",props:"1,1"};
      break;
    case "Phase":
      var variable={name:"Phase",type:"Categorical",mu:0,sd:1,skew:0,kurt:0,
      nlevs:7,iqr:3,
      ncats:2,cases:"before,after",props:"1,1"};
      break;
  }
  return variable;
}

let makeRange = function(min,max,xlog,np) {
  this.min=min;
  this.max=max;
  this.xlog=xlog;
  this.np=np;
}

  let  updateRange = function(value) {
      switch (value) {
        case "n":
           var range={min:10,max:250,xlog:false,np:13};
          break;
        case "rIV":
           var range={min:0,max:0.9,xlog:false,np:13};
          break;
        case "rIV2":
        case "rIVIV2":
        case "rIVIV2DV":
           var range={min:-0.9,max:0.9,xlog:false,np:13};
          break;
        case "IVskew":
        case "DVskew":
        case "Heteroscedasticity":
        case "Dependence":
        case "Outliers":
        case "IVRange":
        case "DVRange":
             var range={min:0,max:1,xlog:false,np:13};
          break;
        case "IVkurtosis":
        case "DVkurtosis":
           var range={min:1.5,max:100000,xlog:true,np:13};
          break;
        case "IVprop":
        case "DVprop":
             var range={min:0.2,max:0.8,xlog:false,np:13};
          break;
        case "IVcats":
        case "DVcats":
             var range={min:2,max:6,xlog:false,np:5};
          break;
        case "IVlevels":
        case "DVlevels":
             var range={min:3,max:10,xlog:false,np:8};
          break;
        case "WithinCorr":
             var range={min:0,max:1,xlog:false,np:13};
          break;
        case "&alpha;":
             var range={min:0.001,max:0.5,xlog:true,np:13};
          break;
        case "Power":
             var range={min:0.1,max:0.9,xlog:false,np:13};
          break;
        case "Repeats":
             var range={min:0,max:8,xlog:false,np:9};
          break;
        case "pNull":
             var range={min:0,max:1,xlog:false,np:13};
          break;
        case "lambda":
             var range={min:0.1,max:1,xlog:false,np:13};
          break;
        case "PoorSamplingAmount":
             var range={min:0,max:0.8,xlog:false,np:13};
          break;
        case "CheatingAmount":
             var range={min:0,max:0.8,xlog:false,np:13};
          break;
        case "ClusterRad":
             var range={min:0,max:1,xlog:false,np:13};
          break;
        case "SampleGamma":
             var range={min:1,max:10,xlog:false,np:13};
          break;
        case "IVType":
             var range={min:"",max:"",xlog:false,np:5};
          break;
        case "DVType":
             var range={min:"",max:"",xlog:false,np:5};
          break;
        case "PDF":
             var range={min:"",max:"",xlog:false,np:7};
          break;
        case "Method":
             var range={min:"",max:"",xlog:false,np:5};
          break;
        case "Usage":
             var range={min:"",max:"",xlog:false,np:2};
          break;
        case "IVRangeC":
             var range={min:0.1,max:3,xlog:false,np:13}
          break;
        case "IVRangeE":
             var range={min:-3,max:3,xlog:false,np:13}
          break;
        case "Cheating":
             var range={min:"",max:"",xlog:false,np:6};
          break;
        case "Transform":
             var range={min:"",max:"",xlog:false,np:3};
          break;
        case "Welch":
             var range={min:"",max:"",xlog:false,np:2};
          break;
        case "Keep":
             var range={min:"",max:"",xlog:false,np:5};
          break;
        case "NoStudies":
             var range={min:2,max:100,xlog:true,np:13};
          break;
        case "MetaType":
             var range={min:"",max:"",xlog:false,np:4};
          break;
        default: 
             var range={min:0,max:1,xlog:false,np:13};
          break;
      }
      return range;
    };

module.exports = events;
