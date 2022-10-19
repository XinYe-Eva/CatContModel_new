#This example uses the Zhang & Luck (2008) model variant ("ZL") 
#Thanks to Dr. Hardman for providing intial code

# 1. Convert study:response pairs to response error (response - study):0 pairs.
# This can be done since the ZL model has no categories and can be fitted with just response error.
# 2. Treat response error in both tasks as linear, with bounded endpoints. Since the endpoints are in
# the tails of the response error distribution, treating data as linear rather than circular should 
# not bias the results.

setwd("~/Desktop/emotion_R") #or wherever you are working

#import library
library(CatContModel)
library(ggplot2)
library(CMBBHT)

#import data
raw_data = read.delim("sub42.txt")

#clean up data, exclude participants
hcData=subset(raw_data,response>0)
hcData=subset(hcData,pnum!=8)
hcData=subset(hcData,pnum!=10)
hcData=subset(hcData,pnum!=13)
hcData=subset(hcData,pnum!=14)
hcData=subset(hcData,pnum!=20)
hcData=subset(hcData,pnum!=21)
hcData=subset(hcData,pnum!=23)
hcData=subset(hcData,pnum!=24)
hcData=subset(hcData,pnum!=26)
hcData=subset(hcData,pnum!=29)
hcData=subset(hcData,pnum!=31)
hcData=subset(hcData,pnum!=32)
hcData=subset(hcData,pnum!=35)
hcData=subset(hcData,pnum!=39)
hcData=subset(hcData,pnum!=40)
hcData=subset(hcData,pnum!=41)

#----------------------------------------

plot(hcData$study, hcData$response)

# Multiply half circle by 2 just so the circular distance function works correctly.
# The data are divided by 2 before fitting.
hcData$study = hcData$study * 2
hcData$response = hcData$response * 2

# Calculate signed response error
clampAngle = function(xs, pm180=FALSE) {
  
  xs = xs %% 360
  
  if (pm180) {
    xs[ xs > 180 ] = xs[ xs > 180 ] - 360
  }
  
  xs
}
circDist = function(xs, ys) {
  clampAngle(ys - xs, pm180=TRUE)
}

hcData$error = circDist(hcData$study, hcData$response)

# Divide the half circle data, including error, by 2 to return to original space
hcData$study = hcData$study / 2
hcData$response = hcData$response / 2
hcData$error = hcData$error / 2

# Put the centered data into new data frames for fitting.
# Set study to 0 to reflect that response error has had study subtracted from it.
hcDataCentered = data.frame(pnum=hcData$pnum, cond=hcData$cond, study=0, response=hcData$error)

hist(hcDataCentered$response)

range(hcData$response)


# Make half circle config treating the data as linear.
hcConfig = list(modelVariant = "ZL", 
                iterations = 11000,
                dataType = "linear", 
                responseRange = c(-90,90))

mhTuning = list()
mhTuning$contSD = 1.5
mhTuning$contSD_cond = 0.41
mhTuning$pMem = 0.4
mhTuning$pMem_cond = 0.1

Z_results = runParameterEstimation(hcConfig, hcDataCentered,mhTuningOverrides = mhTuning)
saveRDS(Z_results, file="ZL_results2_bar180.RDS")
examineMHAcceptance(Z_results)


#--ZL----------------------------------------------------
Z_results = readRDS("ZL_results2_bar180.RDS")
Z_results = removeBurnIn(Z_results, burnIn=1000)


ssz=data.frame() 
for (i in  1:11){
  ssz[i,1] = mean(getParameterPosterior(Z_results, param="pMem", pnum=i, cond='neg'))
  ssz[i,2] = mean(getParameterPosterior(Z_results, param="pMem", pnum=i, cond='neu'))
  ssz[i,3] = mean(getParameterPosterior(Z_results, param="pMem", pnum=i, cond='pos'))
  ssz[i,4] = mean(getParameterPosterior(Z_results, param="contSD", pnum=i, cond='neg'))
  ssz[i,5] = mean(getParameterPosterior(Z_results, param="contSD", pnum=i, cond='neu'))
  ssz[i,6] = mean(getParameterPosterior(Z_results, param="contSD", pnum=i, cond='pos'))
}

pdf("ZL_results.pdf")  # change the file name accordingly
par(cex=1)
plotParameterSummary(Z_results)
dev.off

posteriorMeansAndCredibleIntervals(Z_results)

# Testing Main Effects and Interactions
mei = testMainEffectsAndInteractions(Z_results, doPairwise = TRUE)
# Examine bayes factors in favor of the effect.
mei[ mei$bfType == "10", ]
# Select tests with reasonably strong BFs (either for or against an effect).
mei[ mei$bf > 3, ]

waic = calculateWAIC(Z_results)
waic



