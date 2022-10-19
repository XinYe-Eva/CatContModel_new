library(CatContModel)
library(ggplot2)
library(CMBBHT)
raw_data = read.delim("sub2.txt")
data=subset(raw_data,response>0)
#data=subset(raw_data,pnum!=8)
#data=subset(raw_data,pnum!=10)
#data=subset(raw_data,pnum!=13)
#data=subset(raw_data,pnum!=14)
#data=subset(raw_data,pnum!=20)
#data=subset(raw_data,pnum!=21)
#data=subset(raw_data,pnum!=23)
#data=subset(raw_data,pnum!=24)
#data=subset(raw_data,pnum!=26)
#data=subset(raw_data,pnum!=29)
#data=subset(raw_data,pnum!=31)
#data=subset(raw_data,pnum!=32)
#data=subset(raw_data,pnum!=35)
#data=subset(raw_data,pnum!=39)
#data=subset(raw_data,pnum!=40)
#data=subset(raw_data,pnum!=41)
data=subset(raw_data,pnum!=2)
data=subset(raw_data,pnum!=4)
data=subset(raw_data,pnum!=10)
data=subset(raw_data,pnum!=11)
data=subset(raw_data,pnum!=12)
data=subset(raw_data,pnum!=20)
data=subset(raw_data,pnum!=23)
data=subset(raw_data,pnum!=27)
#----------------------------------------
config = list(iterations=300, 
              maxCategories = 10,
              responseRange = c(1,180))

#within-----------------------
config$modelVariant="withinItem"
mhTuning = list()
mhTuning$catMu = 6
mhTuning$catSD = 1.8
mhTuning$catSelectivity = 4
mhTuning$contSD = 2.5
mhTuning$contSD_cond = 0.7
mhTuning$pCatGuess  = 1.5
mhTuning$pContWithin = 0.3
mhTuning$pContWithin_cond = 0.08
mhTuning$pMem  = 0.3
mhTuning$pMem_cond = 0.1

W_results1 = runParameterEstimation(config, data, mhTuningOverrides = mhTuning)
saveRDS(W_results1, file="WithinItem_results_newmodel3_180.RDS")
examineMHAcceptance(W_results1)

#between------------------------------

config$modelVariant="betweenItem"
mhTuning = list()
mhTuning$catMu = 6
mhTuning$catSD = 1.5
mhTuning$catSelectivity = 2.5#?ĳ?2????
mhTuning$contSD = 1.6
mhTuning$contSD_cond = 0.47
mhTuning$pCatGuess  = 1.7
mhTuning$pContBetween = 0.72
mhTuning$pContBetween_cond = 0.21
mhTuning$pMem  = 0.36
mhTuning$pMem_cond = 0.1

B_results = runParameterEstimation(config, data,mhTuningOverrides = mhTuning)

saveRDS(B_results, file="BetweenItem_results.RDS")
examineMHAcceptance(B_results)

#zl----------------------------
config$modelVariant="ZL"
mhTuning = list()
mhTuning$contSD = 1.5
mhTuning$contSD_cond = 0.41
mhTuning$pMem = 0.3
mhTuning$pMem_cond = 0.1

Z_results = runParameterEstimation(config, data,mhTuningOverrides = mhTuning)
saveRDS(Z_results, file="ZL_results.RDS")
examineMHAcceptance(Z_results)


#=================================================================================

#--WI-------------------------------------------------
W_results = readRDS("WithinItem_results_newmodel2_180.RDS")
W_results = removeBurnIn(W_results, burnIn=250)

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

pdf("WithinItem_results_new2.pdf", 10, 10)  # change the file name accordingly
par(cex=1)
plotParameterSummary(W_results)
dev.off
posteriorMeansAndCredibleIntervals(W_results)
testConditionEffects(W_results)

# Testing Main Effects and Interactions
mei = testMainEffectsAndInteractions(W_results,doPairwise = TRUE)
# Examine bayes factors in favor of the effect.
mei[ mei$bfType == "10", ]
# Select tests with reasonably strong BFs (either for or against an effect).
mei[ mei$bf > 3, ]

waic = calculateWAIC(W_results)
waic


#--BI------------------------------------------------------
B_results = readRDS("BetweenItem_results.RDS")
B_results = removeBurnIn(B_results, burnIn=1000)


ssb=data.frame() 
for (i in  1:11){
  ssb[i,1] = mean(getParameterPosterior(B_results, param="pMem", pnum=i, cond='neu'))
  ssb[i,2] = mean(getParameterPosterior(B_results, param="pMem", pnum=i, cond='neg'))
  ssb[i,3] = mean(getParameterPosterior(B_results, param="contSD", pnum=i, cond='neu'))
  ssb[i,4] = mean(getParameterPosterior(B_results, param="contSD", pnum=i, cond='neg'))
  ssb[i,5] = mean(getParameterPosterior(B_results, param="pContBetween", pnum=i, cond='neu'))
  ssb[i,6] = mean(getParameterPosterior(B_results, param="pContBetween", pnum=i, cond='neg'))
}

pdf("BetweenItem_results.pdf", 10, 10)  # change the file name accordingly
par(cex=1)
plotParameterSummary(B_results)
dev.off

posteriorMeansAndCredibleIntervals(B_results)

# Testing Main Effects and Interactions
mei = testMainEffectsAndInteractions(B_results, doPairwise = TRUE)
# Examine bayes factors in favor of the effect.
mei[ mei$bfType == "10", ]
# Select tests with reasonably strong BFs (either for or against an effect).
mei[ mei$bf > 3, ]

waic = calculateWAIC(B_results)
waic



#--ZL----------------------------------------------------
Z_results = readRDS("ZL_results.RDS")
Z_results = removeBurnIn(Z_results, burnIn=1000)


ssz=data.frame() 
for (i in  1:11){
  ssz[i,1] = mean(getParameterPosterior(Z_results, param="pMem", pnum=i, cond='neg'))
  ssz[i,2] = mean(getParameterPosterior(Z_results, param="pMem", pnum=i, cond='neu'))
  ssz[i,3] = mean(getParameterPosterior(Z_results, param="contSD", pnum=i, cond='neg'))
  ssz[i,4] = mean(getParameterPosterior(Z_results, param="contSD", pnum=i, cond='neu'))
  }

pdf("zl_results.pdf")  # change the file name accordingly
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







#===============================================================
raw_data = read.delim("data_sub10.txt")
subdata=subset(raw_data,response>0)
data=subdata

data$error<-data$response-data$study
for (i in 1:1195){
  if (data$error[i]>90){
    data$error2[i]=(180-data$error[i])}
  
  else if (data$error[i]<(-90)){
    data$error2[i]=180+data$error[i]}
  
  else{data$error2[i]<-data$error[i]}
  
}
data$abs<-abs(data$error2)

# ????ss??errot??  
sdata=data.frame()
rr=data.frame()
ss=data.frame()
for (i in 1:10){
  subdata=subset(data,pnum==i)
  ss[i,1]=sum(subdata$abs)
  if (nrow(subdata)>0){
    line.model = lm(error2~study, data=subdata)
    rr[i,1] = round(summary(line.model)$r.squared,3)
    #if (rr[i,1]>0.2){
     # sdata=rbind(data, subdata)}
  }
}


# ɢ??ͼ error
for (i in 1:10){
  subdata=subset(data,pnum==i)
  
  png(file = paste("sub",as.character(i),".png"),width=600,height=450,units="px",bg = "transparent")
  pic=ggplot(subdata,aes(x=study,y=error2))+
    xlim(1,360)+
    ylim(-180,180)+
    #geom_point()+
    geom_point(size = 4, alpha = .4)+
    # geom_smooth(method=lm)+
    # geom_abline(slope = 1)+
    ggtitle(paste("sub",as.character(i)))+
    theme_classic()
  print(pic)
  dev.off();
}

# ɢ??ͼ response
for (i in 1:11){
  subdata=subset(data,pnum==i)
  
  png(file = paste("sub",as.character(i),".png"),width=600,height=450,units="px",bg = "transparent")
  pic=ggplot(subdata,aes(x=study,y=response))+
    xlim(1,360)+
    ylim(1,360)+
    #geom_point()+
    geom_point(size = 4, alpha = .4)+
     #geom_smooth(method=lm)+
     #geom_abline(slope = 1)+
    ggtitle(paste("sub",as.character(i)))+
    theme_classic()
  print(pic)
  dev.off();
}


ggplot(rr, aes(x = V1)) +
  geom_histogram(binwidth = 0.03)+
  theme_classic()
