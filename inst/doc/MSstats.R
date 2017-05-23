## ---- eval=FALSE---------------------------------------------------------
#  # 'MSstatsInput.csv' is the MSstats report from Skyline.
#  input <- read.csv(file="MSstatsInput.csv")
#  
#  raw <- SkylinetoMSstatsFormat(input)

## ---- eval=FALSE---------------------------------------------------------
#  # Read in MaxQuant files
#  proteinGroups <- read.table("proteinGroups.txt", sep="\t", header=TRUE)
#  
#  infile <- read.table("evidence.txt", sep="\t", header=TRUE)
#  
#  # Read in annotation including condition and biological replicates per run.
#  # Users should make this annotation file. It is not the output from MaxQuant.
#  annot <- read.csv("annotation.csv", header=TRUE)
#  
#  raw <- MaxQtoMSstatsFormat(evidence=infile,
#                             annotation=annot,
#                             proteinGroups=proteinGroups)

## ---- eval=FALSE---------------------------------------------------------
#  input <- read.csv("output_progenesis.csv", stringsAsFactors=F)
#  
#  # Read in annotation including condition and biological replicates per run.
#  # Users should make this annotation file. It is not the output from Progenesis.
#  annot <- read.csv('annotation.csv')
#  
#  raw <- ProgenesistoMSstatsFormat(input, annotation=annot)

## ---- eval=FALSE---------------------------------------------------------
#  input <- read.csv("output_spectronaut.csv", stringsAsFactors=F)
#  
#  quant <- SpectronauttoMSstatsFormat(input)

## ---- eval=FALSE---------------------------------------------------------
#  QuantData <- dataProcess(SRMRawData)

## ---- eval=FALSE---------------------------------------------------------
#  QuantData <- dataProcess(SRMRawData)
#  
#  # Profile plot
#  dataProcessPlots(data=QuantData, type="ProfilePlot")
#  
#  # Quality control plot
#  dataProcessPlots(data=QuantData, type="QCPlot")	
#  
#  # Quantification plot for conditions
#  dataProcessPlots(data=QuantData, type="ConditionPlot")

## ---- eval=FALSE---------------------------------------------------------
#  QuantData <- dataProcess(SRMRawData)
#  
#  levels(QuantData$ProcessedData$GROUP_ORIGINAL)
#  comparison <- matrix(c(-1,0,0,0,0,0,1,0,0,0), nrow=1)
#  row.names(comparison) <- "T7-T1"
#  
#  # Tests for differentially abundant proteins with models:
#  testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)

## ---- eval=FALSE---------------------------------------------------------
#  QuantData <- dataProcess(SRMRawData)
#  
#  # based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
#  comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
#  comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
#  comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
#  comparison<-rbind(comparison1,comparison2, comparison3)
#  row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")
#  
#  testResultMultiComparisons <- groupComparison(contrast.matrix=comparison, data=QuantData)
#  
#  # Volcano plot
#  groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="VolcanoPlot")
#  
#  # Heatmap
#  groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="Heatmap")
#  
#  # Comparison Plot
#  groupComparisonPlots(data=testResultMultiComparisons$ComparisonResult, type="ComparisonPlot")

## ---- eval=FALSE---------------------------------------------------------
#  testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)
#  
#  # normal quantile-quantile plots
#  modelBasedQCPlots(data=testResultOneComparison, type="QQPlots")
#  
#  # residual plots
#  modelBasedQCPlots(data=testResultOneComparison, type="ResidualPlots")

## ---- eval=FALSE---------------------------------------------------------
#  QuantData <- dataProcess(SRMRawData)
#  head(QuantData$ProcessedData)
#  
#  ## based on multiple comparisons  (T1 vs T3; T1 vs T7; T1 vs T9)
#  comparison1 <- matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
#  comparison2 <- matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
#  comparison3 <- matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
#  comparison <- rbind(comparison1,comparison2, comparison3)
#  row.names(comparison) <- c("T3-T1","T7-T1","T9-T1")
#  
#  testResultMultiComparisons <- groupComparison(contrast.matrix=comparison,data=QuantData)
#  
#  #(1) Minimal number of biological replicates per condition
#  designSampleSize(data=testResultMultiComparisons$fittedmodel, numSample=TRUE,
#    desiredFC=c(1.25,1.75), FDR=0.05, power=0.8)
#  
#  #(2) Power calculation
#  designSampleSize(data=testResultMultiComparisons$fittedmodel, numSample=2,
#    desiredFC=c(1.25,1.75), FDR=0.05, power=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  # (1) Minimal number of biological replicates per condition
#  result.sample <- designSampleSize(data=testResultMultiComparisons$fittedmodel, numSample=TRUE,
#                                  desiredFC=c(1.25,1.75), FDR=0.05, power=0.8)
#  designSampleSizePlots(data=result.sample)
#  
#  # (2) Power
#  result.power <- designSampleSize(data=testResultMultiComparisons$fittedmodel, numSample=2,
#                                 desiredFC=c(1.25,1.75), FDR=0.05, power=TRUE)
#  designSampleSizePlots(data=result.power)

## ---- eval=FALSE---------------------------------------------------------
#  QuantData <- dataProcess(SRMRawData)
#  
#  # Sample quantification
#  sampleQuant <- quantification(QuantData)
#  
#  # Group quantification
#  groupQuant <- quantification(QuantData, type="Group")

## ---- eval=FALSE---------------------------------------------------------
#  # Consider data from a spiked-in contained in an example dataset
#  head(SpikeInDataNonLinear)
#  
#  nonlinear_quantlim_out <- nonlinear_quantlim(SpikeInDataNonLinear)
#  
#  # Get values of LOB/LOD
#  nonlinear_quantlim_out$LOB[1]
#  nonlinear_quantlim_out$LOD[1]

## ---- eval=FALSE---------------------------------------------------------
#  # Consider data from a spiked-in contained in an example dataset
#  head(SpikeInDataLinear)
#  
#  linear_quantlim_out <- linear_quantlim(SpikeInDataLinear)
#  
#  # Get values of LOB/LOD
#  linear_quantlim_out$LOB[1]
#  linear_quantlim_out$LOD[1]

## ---- eval=FALSE---------------------------------------------------------
#  # Consider data from a spiked-in contained in an example dataset
#  head(SpikeInDataNonLinear)
#  
#  nonlinear_quantlim_out <- nonlinear_quantlim(SpikeInDataNonLinear, alpha = 0.05)
#  
#  #Get values of LOB/LOD
#  nonlinear_quantlim_out$LOB[1]
#  nonlinear_quantlim_out$LOD[1]
#  
#  plot_quantlim(spikeindata = SpikeInDataLinear, quantlim_out  = nonlinear_quantlim_out,
#  dir_output =  getwd(), alpha = 0.05)

