################################################################################
# This program performs statistical tests for stroke paper.
# Statistical results are published to output file created by sink()

################################################################################
# INSTRUCTIONS:
#     Change DIR, data, and the sink file before running in the terminal
#     To run in the terminal
#     $ R
#     Then run the script by executing the following
#     > source('Stroke-lmer.r')
#     If running from console in RStudio on windows computer, use setwd
#     for example - setwd('C:/Users/numur/Desktop/darpa_dataanalysis/src')
################################################################################

options(contrasts=c("contr.sum","contr.poly"))
remove(list = ls())

# parameters 
DIR = '/home/milli/Desktop/DynamicTasks/DataAnalysis'

subjects = list("Sub202","Sub203","Sub208","Sub209","Sub211","Sub212","Sub214","Sub218","Sub219","Sub220")
run_indiv_subjects = 0 # 1 is TRUE and 0 is FALSE

# loads packages and dependencies
library(ez)
library(car)
library(lme4)
library(multcomp)
library(rstatix)

################################################################################
################################################################################
#               Import Data 
################################################################################
################################################################################
data_all = read.csv(paste(DIR,paste("stroke-freq",".csv",sep=""),sep="/"))
data_all = subset(data_all, EResonance!=0)
data_modsev = subset(data_all, FMA<40)
# data_modsev = data_all

# define file to save data to
sink(paste(DIR,"Stats",paste("e-at-res-stats-mod-sev.txt",sep="-"),sep="/"))

cat("\n")
cat("################################################################################ \n")
cat("###############            Energy at Resonance Metric            ############### \n")
cat("################################################################################ \n")

cat("\n")
cat("################################################################################ \n")
cat("Test for normality: Shapiro test\n")
cat("################################################################################ \n")
normality = data_modsev %>%
  group_by(Arm,Loading,BallFreq) %>%
  shapiro_test(EResonance)
print(normality)


cat("\n")
cat("###############            All Experimental Factors            ############### \n")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Arm*Loading*BallFreq + (1+Arm*Loading*BallFreq|Subject),
                  data = data_modsev)
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_modsev,EResonance,
                wid = .(Subject),
                within = .(Arm,Loading,BallFreq),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)


cat("\n\n\n")
cat("###############            0.5Hz Trials           ############### \n")
data_freq = subset(data_modsev, BallFreq=="0.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Arm*Loading + (1+Arm*Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Arm,Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            1Hz Trials           ############### \n")
data_freq = subset(data_modsev, BallFreq=="1Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Arm*Loading + (1+Arm*Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Arm,Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            1.5Hz Trials           ############### \n")
data_freq = subset(data_modsev, BallFreq=="1.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Arm*Loading + (1+Arm*Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Arm,Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            2.5Hz Trials           ############### \n")
data_freq = subset(data_modsev, BallFreq=="2.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Arm*Loading + (1+Arm*Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Arm,Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("################################################################################ \n")
cat("##########    Energy at Resonance Metric Within the Paretic Arm      ########### \n")
cat("################################################################################ \n")
data = subset(data_modsev, Arm=="paretic")

cat("\n")
cat("###############            All Experimental Factors            ############### \n")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading*BallFreq + (1+Loading*BallFreq|Subject),
                  data = data,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_modsev,EResonance,
                wid = .(Subject),
                within = .(Loading,BallFreq),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            0.5Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="0.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            1Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="1Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            1.5Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="1.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            2.5Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="2.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("################################################################################ \n")
cat("##########   Energy at Resonance Metric Within the Non-Paretic Arm   ########### \n")
cat("################################################################################ \n")
data = subset(data_modsev, Arm=="nonparetic")

cat("\n")
cat("###############            All Experimental Factors            ############### \n")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading*BallFreq + (1+Loading*BallFreq|Subject),
                  data = data,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_modsev,EResonance,
                wid = .(Subject),
                within = .(Loading,BallFreq),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            0.5Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="0.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            1Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="1Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            1.5Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="1.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n\n\n")
cat("###############            2.5Hz Trials           ############### \n")
data_freq = subset(data, BallFreq=="2.5Hz")
cat("Linear Mixed Model \n")
lmer_model = lmer(EResonance ~  Loading + (1+Loading|Subject),
                  data = data_freq,
                  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
anov = Anova(lmer_model,type="II")
print(anov)
cat("\nANOVA (although ~23 trials are missing) \n")
mod.ez<-ezANOVA(data_freq,EResonance,
                wid = .(Subject),
                within = .(Loading),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

if (run_indiv_subjects==1) {
  for(i in 1:length(subjects)) {
    data_subject = subset(data_all, Subject==subjects[i])

    # define file to save data to
    sink(paste(DIR,"Stats",paste(subjects[i],"e-at-res-stats.txt",sep="-"),sep="/"))

    cat("\n")
    cat("################################################################################ \n")
    cat("###############            Energy at Resonance Metric            ############### \n")
    cat("################################################################################ \n")

    cat("\n")
    cat("###############            All Experimental Factors            ############### \n")
    attach(data_subject)
    aov_results <- aov(EResonance~Arm*Loading*BallFreq)
    print(summary(aov_results))
    detach(data_subject)

    cat("\n\n\n")
    cat("###############            0.5Hz Trials           ############### \n")
    data_freq = subset(data_subject, BallFreq=="0.5Hz")
    attach(data_freq)
    aov_results <- aov(EResonance~Arm*Loading)
    print(summary(aov_results))
    detach(data_freq)

    cat("\n\n\n")
    cat("###############            1Hz Trials           ############### \n")
    data_freq = subset(data_subject, BallFreq=="1Hz")
    attach(data_freq)
    aov_results <- aov(EResonance~Arm*Loading)
    print(summary(aov_results))
    detach(data_freq)

    cat("\n\n\n")
    cat("###############            1.5Hz Trials           ############### \n")
    data_freq = subset(data_subject, BallFreq=="1.5Hz")
    attach(data_freq)
    aov_results <- aov(EResonance~Arm*Loading)
    print(summary(aov_results))
    detach(data_freq)

    cat("\n\n\n")
    cat("###############            2.5Hz Trials           ############### \n")
    data_freq = subset(data_subject, BallFreq=="2.5Hz")
    attach(data_freq)
    aov_results <- aov(EResonance~Arm*Loading)
    print(summary(aov_results))
    detach(data_freq)

    cat("\n")
    cat("################################################################################ \n")
    cat("##########    Energy at Resonance Metric Within the Paretic Arm      ########### \n")
    cat("################################################################################ \n")
    data_subject = subset(data_subject, Arm=="paretic")

    cat("\n\n\n")
    cat("###############            1.5Hz Trials           ############### \n")
    data_freq = subset(data_subject, BallFreq=="1.5Hz")
    attach(data_freq)
    aov_results <- aov(EResonance~Loading)
    print(summary(aov_results))
    detach(data_freq)
  }
}