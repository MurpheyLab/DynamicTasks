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
subjects = list("Sub202","Sub203","Sub208","Sub209","Sub211","Sub212")

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
data = read.csv(paste(DIR,paste("stroke-freq",".csv",sep=""),sep="/"))


for(i in 1:length(subjects)) {
  # define file to save data to
  sink(paste(DIR,"Stats",paste(subjects[i],"tests.txt",sep="-"),sep="/"))
  
  cat("Energy at Resonance Metric  \n")
  data_subject = subset(data, Subject==subjects[i])
  lmer_model = lm(EResonance ~  Arm*Loading*BallFreq,
                      data = data_subject)
  anov = Anova(lmer_model,type="II")
  print(anov)
  
  cat("\n\n\n")
  cat("0.5Hz \n")
  data_freq = subset(data_subject, BallFreq=="0.5Hz")
  lmer_model = lm(EResonance ~  Arm*Loading,
                  data = data_freq)
  anov = Anova(lmer_model,type="II")
  print(anov)
  cat("\n")
  cat("Paretic \n")
  data_A0 = subset(data_freq, Arm=="paretic")
  lmer_model = lm(EResonance ~  Loading,
                    data = data_A0)
  anov = Anova(lmer_model,type="II")
  print(anov)
  cat("\n")
  cat("Non-Paretic \n")
  data_A1 = subset(data_freq, Arm=="nonparetic")
  lmer_model = lm(EResonance ~  Loading,
                    data = data_A1)
  anov = Anova(lmer_model,type="II")
  print(anov)
  
  cat("\n")
  cat("1Hz \n")
  data_freq = subset(data_subject, BallFreq=="1Hz")
  lmer_model = lm(EResonance ~  Arm*Loading,
                  data = data_freq)
  anov = Anova(lmer_model,type="II")
  print(anov)
  cat("\n")
  cat("Paretic \n")
  data_A0 = subset(data_freq, Arm=="paretic")
  lmer_model = lm(EResonance ~  Loading,
                  data = data_A0)
  anov = Anova(lmer_model,type="II")
  print(anov)
  cat("\n")
  cat("Non-Paretic \n")
  data_A1 = subset(data_freq, Arm=="nonparetic")
  lmer_model = lm(EResonance ~  Loading,
                  data = data_A1)
  anov = Anova(lmer_model,type="II")
  print(anov)
  
  cat("\n")
  cat("1.5Hz \n")
  data_freq = subset(data_subject, BallFreq=="1.5Hz")
  lmer_model = lm(EResonance ~  Arm*Loading,
                  data = data_freq)
  anov = Anova(lmer_model,type="II")
  print(anov)
  cat("\n")
  cat("Paretic \n")
  data_A0 = subset(data_freq, Arm=="paretic")
  lmer_model = lm(EResonance ~  Loading,
                  data = data_A0)
  anov = Anova(lmer_model,type="II")
  print(anov)
  cat("\n")
  cat("Non-Paretic \n")
  data_A1 = subset(data_freq, Arm=="nonparetic")
  lmer_model = lm(EResonance ~  Loading,
                  data = data_A1)
  anov = Anova(lmer_model,type="II")
  print(anov)
  
  cat("\n")
  cat("2.5Hz \n")
  data_freq = subset(data_subject, BallFreq=="2.5Hz")
  lmer_model = lm(EResonance ~  Arm*Loading,
                  data = data_freq)
  anov = Anova(lmer_model,type="II")
  print(anov) 
  cat("\n")
  cat("Paretic \n")
  data_A0 = subset(data_freq, Arm=="paretic")
  lmer_model = lm(EResonance ~  Loading,
                  data = data_A0)
  anov = Anova(lmer_model,type="II")
  print(anov)
  cat("\n")
  cat("Non-Paretic \n")
  data_A1 = subset(data_freq, Arm=="nonparetic")
  lmer_model = lm(EResonance ~  Loading,
                  data = data_A1)
  anov = Anova(lmer_model,type="II")
  print(anov)
  
}

# define file to save data to
sink(paste(DIR,"Stats",paste("all-trial-tests.txt",sep="-"),sep="/"))

cat("Energy at Resonance Metric \n")
lmer_model = lmer(EResonance ~  Arm*Loading*BallFreq + (1|Subject),
                  data = data)
#print(summary(lmer_model), correlation=TRUE)
# lmer_model = lmer(EResonance ~  Arm*Loading*BallFreq + ((Subject)/(Arm*Loading*BallFreq)),
#                 data = data)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("0.5Hz \n")
data_freq = subset(data, BallFreq=="0.5Hz")
lmer_model = lmer(EResonance ~  Arm*Loading + (1|Subject),
                data = data_freq)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Paretic \n")
data_A0 = subset(data_freq, Arm=="paretic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                data = data_A0)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Non-Paretic \n")
data_A1 = subset(data_freq, Arm=="nonparetic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                data = data_A1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("1Hz \n")
data_freq = subset(data, BallFreq=="1Hz")
lmer_model = lmer(EResonance ~  Arm*Loading + (1|Subject),
                data = data_freq)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Paretic \n")
data_A0 = subset(data_freq, Arm=="paretic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                  data = data_A0)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Non-Paretic \n")
data_A1 = subset(data_freq, Arm=="nonparetic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                  data = data_A1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("1.5Hz \n")
data_freq = subset(data, BallFreq=="1.5Hz")
lmer_model = lmer(EResonance ~  Arm*Loading + (1|Subject),
                data = data_freq)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Paretic \n")
data_A0 = subset(data_freq, Arm=="paretic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                  data = data_A0)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Non-Paretic \n")
data_A1 = subset(data_freq, Arm=="nonparetic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                  data = data_A1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("2.5Hz \n")
data_freq = subset(data, BallFreq=="2.5Hz")
lmer_model = lmer(EResonance ~  Arm*Loading + (1|Subject),
                data = data_freq)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Paretic \n")
data_A0 = subset(data_freq, Arm=="paretic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                  data = data_A0)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\n")
cat("Non-Paretic \n")
data_A1 = subset(data_freq, Arm=="nonparetic")
lmer_model = lmer(EResonance ~  Loading + (1|Subject),
                  data = data_A1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("Paretic \n")
data_A0 = subset(data, Arm=="paretic")
lmer_model = lmer(EResonance ~  BallFreq*Loading + (1|Subject),
                  data = data_A0)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("Non-Paretic \n")
data_A1 = subset(data, Arm=="nonparetic")
lmer_model = lmer(EResonance ~  BallFreq*Loading + (1|Subject),
                  data = data_A1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("Loading \n")
data_SL1 = subset(data, Loading=="35%")
lmer_model = lmer(EResonance ~  BallFreq*Arm + (1|Subject),
                  data = data_SL1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("No Loading \n")
data_SL0 = subset(data, Loading=="0%")
lmer_model = lmer(EResonance ~  BallFreq*Arm + (1|Subject),
                  data = data_SL0)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("Paretic - Loading \n")
data_A0SL1 = subset(data_A0, Loading=="35%")
lmer_model = lmer(EResonance ~  BallFreq + (1|Subject),
                  data = data_A0SL1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("Paretic - No Loading \n")
data_A0SL0 = subset(data_A0, Loading=="0%")
lmer_model = lmer(EResonance ~  BallFreq + (1|Subject),
                  data = data_A0SL0)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("Non-Paretic - Loading \n")
data_A1SL1 = subset(data_A1, Loading=="35%")
lmer_model = lmer(EResonance ~  BallFreq + (1|Subject),
                  data = data_A1SL1)
anov = Anova(lmer_model,type="II")
print(anov)

cat("\n\n\n")
cat("Non-Paretic - No Loading \n")
data_A1SL0 = subset(data_A1, Loading=="0%")
lmer_model = lmer(EResonance ~  BallFreq + (1|Subject),
                  data = data_A1SL0)
anov = Anova(lmer_model,type="II")
print(anov)