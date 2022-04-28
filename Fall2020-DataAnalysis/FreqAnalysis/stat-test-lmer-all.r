################################################################################
# This program performs statistical tests for darpa HST overall performance data.
# Statistical results are published to output file created by sink()

################################################################################
# INSTRUCTIONS:
#     Change DIR, data, and the sink file before running in the terminal
#     To run in the terminal
#     $ R
#     Then run the script by executing the following
#     > source('stat-test-glmer.r')
#     If running from console in RStudio on windows computer, use setwd
#     for example - setwd('C:/Users/numur/Desktop/darpa_dataanalysis/src')
################################################################################
# COMMENTS:
# - To test a different metric, Ctrl replace the name of the metric


################################################################################
# useful link https://www.youtube.com/watch?v=wWjVQ8eqR9k
options(contrasts=c("contr.sum","contr.poly"))
remove(list = ls())

# parameters 
DIR = 'C:/Users/milli/OneDrive/Documents/DowntownGames-Biodex/Fall2020-DataAnalysis/FreqAnalysis'
subjects = list("Sub202","Sub203","Sub208","Sub209","Sub211","Sub212","Sub214")
# subjects = list("Sub208","Sub211","Sub212","Sub214")

set = 'all' # 'all' or 'combined'
# loads packages and dependencies
library(ez)
library(car) 
library(lme4)
library(multcomp)

################################################################################
################################################################################
#               Import Data 
################################################################################
################################################################################
if (set=='all') {
  data = read.csv(paste(DIR,paste("stroke-freq-all",".csv",sep=""),sep="/"))
} else {
  data = read.csv(paste(DIR,paste("stroke-freq-average",".csv",sep=""),sep="/"))
}

# data = subset(data, Subject!="Sub202" & Subject!="Sub203" & Subject!="Sub209")
data = subset(data, BallFreq!="0.5Hz")

data_SL1 = subset(data, SL1!=1000)
data_SL0 = subset(data, SL0!=1000)
data_A0 = subset(data, A0!=1000)
data_A1 = subset(data, A1!=1000)

# define file to save data to
sink(paste(DIR,"Stats",paste("percent-tests.txt",sep="-"),sep="/"))

cat("Percent decrease in paretic arm with loading \n")
if (set=='all') {
  lmer_model = lm(SL1 ~  BallFreq*Subject,# + (1|Subject),
                    data = data_SL1)
} else {
  lmer_model = lm(SL1 ~  BallFreq,# + (1|Subject),
                  data = data_SL1)
}
anov = Anova(lmer_model,type="II")
print(anov)
posthoc<-pairwise.t.test(data_SL1$SL1,
                         data_SL1$BallFreq,
                         paired = FALSE,
                         p.adjust.method = "bonferroni",detailed = TRUE)
print(posthoc)

cat("\n\n\n")
cat("########################################################################### \n")
cat("Percent decrease in paretic arm with no loading \n")
if (set=='all') {
  lmer_model = lm(SL0 ~  BallFreq*Subject,# + (1|Subject),
                  data = data_SL0)
} else {
  lmer_model = lm(SL0 ~  BallFreq,# + (1|Subject),
                  data = data_SL0)
}
anov = Anova(lmer_model,type="II")
print(anov)
posthoc<-pairwise.t.test(data_SL0$SL0,
                         data_SL0$BallFreq,
                         paired = FALSE,
                         p.adjust.method = "bonferroni",detailed = TRUE)
print(posthoc)

cat("\n\n\n")
cat("########################################################################### \n")
cat("Percent decrease due to loading in paretic arm \n")
if (set=='all') {
  lmer_model = lm(A0 ~  BallFreq*Subject,# + (1|Subject),
                  data = data_A0)
} else {
  lmer_model = lm(A0 ~  BallFreq,# + (1|Subject),
                  data = data_A0)
}
anov = Anova(lmer_model,type="II")
print(anov)
posthoc<-pairwise.t.test(data_A0$A0,
                         data_A0$BallFreq,
                         paired = FALSE,
                         p.adjust.method = "bonferroni",detailed = TRUE)
print(posthoc)


cat("\n\n\n")
cat("########################################################################### \n")
cat("Percent decrease due to loading in nonparetic arm \n")
if (set=='all') {
  lmer_model = lm(A1 ~  BallFreq*Subject,# + (1|Subject),
                  data = data_A1)
} else {
  lmer_model = lm(A1 ~  BallFreq,# + (1|Subject),
                  data = data_A1)
}
anov = Anova(lmer_model,type="II")
print(anov)
posthoc<-pairwise.t.test(data_A1$A1,
                         data_A1$BallFreq,
                         paired = FALSE,
                         p.adjust.method = "bonferroni",detailed = TRUE)
print(posthoc)

if (set=='all') {
for(i in 1:length(subjects)) {
  # define file to save data to
  sink(paste(DIR,"Stats",paste(subjects[i],"percent-tests.txt",sep="-"),sep="/"))
  
  cat("Percent decrease in paretic arm with loading \n")
  data_SL1_subject = subset(data_SL1, Subject==subjects[i])
  lmer_model = lm(SL1 ~  BallFreq,# + (1|Subject),
                  data = data_SL1_subject)
  anov = Anova(lmer_model,type="II")
  print(anov)
  posthoc<-pairwise.t.test(data_SL1_subject$SL1,
                           data_SL1_subject$BallFreq,
                           paired = FALSE,
                           p.adjust.method = "bonferroni",detailed = TRUE)
  print(posthoc)
  
  cat("\n\n\n")
  cat("Percent decrease in paretic arm with no loading \n")
  data_SL0_subject = subset(data_SL0, Subject==subjects[i])
  lmer_model = lm(SL0 ~  BallFreq,# + (1|Subject),
                  data = data_SL0_subject)
  anov = Anova(lmer_model,type="II")
  print(anov)
  
  cat("\n\n\n")
  cat("Percent decrease due to loading in paretic arm \n")
  data_A0_subject = subset(data_A0, Subject==subjects[i])
  lmer_model = lm(A0 ~  BallFreq,# + (1|Subject),
                  data = data_A0_subject)
  anov = Anova(lmer_model,type="II")
  print(anov)
  posthoc<-pairwise.t.test(data_SL1_subject$A0,
                           data_SL1_subject$BallFreq,
                           paired = FALSE,
                           p.adjust.method = "bonferroni",detailed = TRUE)
  print(posthoc)
  
  cat("\n\n\n")
  cat("Percent decrease due to loading in nonparetic arm \n")
  data_A1_subject = subset(data_A1, Subject==subjects[i])
  lmer_model = lm(A1 ~  BallFreq,# + (1|Subject),
                  data = data_A1_subject)
  anov = Anova(lmer_model,type="II")
  print(anov)
  
}
}
