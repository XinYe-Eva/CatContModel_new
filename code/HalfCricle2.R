# This is the script for analyzing [0,180] bar dataset.
# I used the dataset sub42 gathered by XinYe 2021.
# Thanks to Dr. Hardman for providing the initial code.
# Method2

setwd("~/Desktop/emotion_R") #or wherever you are working

#import library
library(CatContModel)
library(ggplot2)
library(CMBBHT)

#import data
raw_data = read.delim("sub42.txt")

#clean up data, exclude participants
data=subset(raw_data,response>0)
data=subset(data,pnum!=8)
data=subset(raw_data,pnum!=10)
data=subset(raw_data,pnum!=13)
data=subset(raw_data,pnum!=14)
data=subset(raw_data,pnum!=20)
data=subset(raw_data,pnum!=21)
data=subset(raw_data,pnum!=23)
data=subset(raw_data,pnum!=24)
data=subset(raw_data,pnum!=26)
data=subset(raw_data,pnum!=29)
data=subset(raw_data,pnum!=31)
data=subset(raw_data,pnum!=32)
data=subset(raw_data,pnum!=35)
data=subset(raw_data,pnum!=39)
data=subset(raw_data,pnum!=40)
data=subset(raw_data,pnum!=41)

#----------------------------------------
# Data Visualization
head(data)
# Data scatterplot
plot(data$study, data$response, pch=20, col=rgb(0,0,0,0.3))

# Response histogram
hist(data$response, breaks=100)


########################
# Solution 2: Modify the data by multiplying by 2, changing the range from [0,180] to [0,360]
#
# This solution is good in that it treats the data as circular, just on a half-size circle.
#
# This solution has downsides when it comes to describing the results.
# In particular, because the range of the data has been increased by a factor of 2,
# the standard deviation parameters (contSD, catSD, and catSelectivity) and catMu
# will be double their true values.
#
# For plots and other descriptive statistics, these parameters must be scaled by dividing by 2.
# For hypothesis tests, however, the original, undivided parameter values must be used.
# The code below gives more detail.

data2x = data # copy the data before modification

# Modify the data by multiplying by 2
data2x$study = data2x$study * 2
data2x$response = data2x$response * 2

# Check the ranges and plot the data
range(data2x$study)
range(data2x$response)
plot(data2x$study, data2x$response)

# Set up a basic configuration and run the model
config2x = list(modelVariant="withinItem", 
                iterations=5, 
                maxCategories=10,
                responseRange = c(1,360))

#within-----------------------
#configLin$modelVariant="withinItem"
# 0.4-0.6
mhTuning = list()
mhTuning$catMu = 15
mhTuning$catSD = 8
mhTuning$catSelectivity = 6.5
mhTuning$contSD = 3.9
mhTuning$contSD_cond = 1.1
mhTuning$pCatGuess  = 1.9
mhTuning$pContWithin = 0.4
mhTuning$pContWithin_cond = 0.09
mhTuning$pMem  = 0.4
mhTuning$pMem_cond = 0.15

W_results = runParameterEstimation(config2x, data2x, mhTuningOverrides = mhTuning)
#saveRDS(W_results, file="WithinItem_results_bar180_method2.RDS")
examineMHAcceptance(W_results)

#=================================================================================

#--WI-------------------------------------------------
W_results = readRDS("WithinItem_results_bar180_method2.RDS")

# Brief convergence analysis and remove burn in
plot(res2x$posteriors[[ "pMem_cond[3]" ]], type='l')
plot(res2x$posteriors[[ "pContBetween_cond[2]" ]], type='l')
W_results = removeBurnIn(W_results, burnIn=1000)

ssw=data.frame() 
for (i in  1:11){
  ssw[i,1] = mean(getParameterPosterior(W_results, param="pMem", pnum=i, cond='pos'))
  ssw[i,2] = mean(getParameterPosterior(W_results, param="pMem", pnum=i, cond='neu'))
  ssw[i,3] = mean(getParameterPosterior(W_results, param="pMem", pnum=i, cond='neg'))
  ssw[i,4] = mean(getParameterPosterior(W_results, param="contSD", pnum=i, cond='pos'))
  ssw[i,5] = mean(getParameterPosterior(W_results, param="contSD", pnum=i, cond='neu'))
  ssw[i,6] = mean(getParameterPosterior(W_results, param="contSD", pnum=i, cond='neg'))
  ssw[i,7] = mean(getParameterPosterior(W_results, param="pContWithin", pnum=i, cond='pos'))
  ssw[i,8] = mean(getParameterPosterior(W_results, param="pContWithin", pnum=i, cond='neu'))
  ssw[i,9] = mean(getParameterPosterior(W_results, param="pContWithin", pnum=i, cond='neg'))
}


####
# [0,360] degrees
clampAngle_R = function(xs, pm180=FALSE) {
  
  xs = xs %% 360
  
  if (pm180) {
    xs[ xs > 180 ] = xs[ xs > 180 ] - 360
  }
  
  xs
}

# signed distance in degrees
circDist_R = function(xs, ys) {
  clampAngle_R(ys - xs, pm180=TRUE)
}
#####
# Descriptive Analysis
# To describe the standard deviation and catMu parameter posteriors with respect to the 
# original half-circle data space, those parameter posteriors must be divided by 2.

# This function to convert posterior parameter chains to/from half circle space
# The standard deviation parameters (contSD, catSD, and catSelectivity) and catMu 
# posteriors will be multiplied by multiple.

halfCircle_convertParameters = function(res, multiple) {
  
  fullParamNames = names(res$posteriors)
  
  # Select standard deviation parameters from all parameters
  sdParamNames = getSdParams(res)
  fullSdParam = NULL
  for (sdp in sdParamNames) {
    matchSd = grepl(sdp, fullParamNames, fixed=TRUE)
    fullSdParam = c(fullSdParam, fullParamNames[ matchSd ])
  }
  
  # Modify standard deviation parameters
  for (sdp in fullSdParam) {
    
    post = res$posteriors[[ sdp ]]
    
    if (grepl(".var", sdp, fixed=TRUE)) {
      # Variance parameters are in squared units, so sqrt, divide, then square
      post = sqrt(post)
      post = post * multiple
      post = post^2
      # Adjusting variance parameters like this may not be totally correct, but approximate.
    } else {
      # Participant and mean (.mu) parameters divide by 2
      post = post * multiple
    }
    
    res$posteriors[[ sdp ]] = post
  }
  
  # Modify catMu
  fullCatMu = fullParamNames[ grepl("catMu", fullParamNames, fixed=TRUE) ]
  for (cm in fullCatMu) {
    
    post = res$posteriors[[ cm ]]
    
    # The catMu parameters are free to wander outside of [0,360].
    # Bring them back before multiplying
    #post = CatContModel::clampAngle(post, pm180 = FALSE, degrees = TRUE)
    post = clampAngle_R(post, pm180 = FALSE)
    
    post = post * multiple
    
    res$posteriors[[ cm ]] = post
  }
  
  res
}

# res2x was estimated with data that were multiplied by 2 to get to the full circle space.
# Get back to the half circle space by dividing by 2 (i.e. multiply by 1/2)
res2x_conv = halfCircle_convertParameters(W_results, multiple=1/2)

# Means and credible intervals are descriptive and use converted parameters
posteriorMeansAndCredibleIntervals(res2x_conv)

# Plots are also descriptive and use converted parameters.
# The catMu plot has some problems, but is still readable
plotParameterSummary(res2x_conv)

# You can plot individual panels of the overall plot with plotParameter
plotParameter(res2x_conv, param="contSD")
plotParameter(res2x_conv, param="catSelectivity")
plotParameter(res2x_conv, param="catSD")
plotParameter(res2x_conv, param="pMem")
plotParameter(res2x_conv, param="pContWithin")
plotParameter(res2x_conv, param="pCatGuess")

# This function is a possible way to plot catMu for half circle data
plotCatMu_halfCircle = function(res, precision = 1, pnums = res$pnums) {
  
  post = convertPosteriorsToMatrices(res, param = c("catMu", "catActive"))
  
  catMu = post$catMu[ pnums, , ]
  catActive = post$catActive[ pnums, , ]
  
  plotData = CatContModel:::catMu_getPlotData(catMu, catActive, precision = precision, 
                                              dataType = res$config$dataType, 
                                              responseRange = res$config$responseRange)
  
  # Ignore all points past 180
  plotData = plotData[ plotData$x <= 180, ]
  
  # The last point is the first point, so copy the value from the first point.
  plotData$y[ plotData$x == 180 ] = plotData$y[1]
  
  CatContModel:::catMu_plotSetup(plotData, res$config$dataType)
  CatContModel:::catMu_plotData(plotData$x, plotData$y, col = plotData$color, type = "polygon")
  
  invisible(plotData)
}

plotCatMu_halfCircle(res2x_conv)


###
# Although converted parameters are used for describing the results, other kinds of 
# analysis will not work correctly if you use the converted parameter values.
# Here are some major examples of where you should use the original parameter values.

# Hypothesis tests depend on priors, so use unconverted parameters
# Testing Main Effects and Interactions
mei = testMainEffectsAndInteractions(W_results,doPairwise = TRUE)
# Examine bayes factors in favor of the effect.
mei[ mei$bfType == "10", ]
# Select tests with reasonably strong BFs (either for or against an effect).
mei[ mei$bf > 3, ]


ce = testConditionEffects(W_results)
ce[ ce$bfType == "10", ]


# Posterior predictive (sampling data based on parameters) uses unconverted parameters
posteriorPredictivePlot(W_results)

# WAIC also uses unconverted parameters
calculateWAIC(W_results)

# There are many functions in this package that summarize results.
# When using each, think carefully about whether the function wants
# the original real posteriors or converted posteriors.
# Anything that uses priors must use original posteriors.
