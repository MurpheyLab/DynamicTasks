################################################################################
# Hypothesis 2 determines whether haptic forces significantly impacts energy@resonance (comparing F0_B1 and F1_B1)
# Hypothesis 3 determines whether a moving ball significantly impacts energy@resonance (comparing F1_B0 and F1_B1)

# This program performs ANOVA repeated measures statistical tests for darpa HST data.
# Statistical results are published to output file created by sink()

# The exANOVA function in the ez package is used to perform statistical tests.
# This package and dependencies will need to be installed.
# To install run : install.package("ez")
################################################################################
# INSTRUCTIONS:
#     Change DIR, data, and the sink file before running in the terminal
#     To run in the terminal
#     $ R
#     Then run the script by executing the following
#     > source('stat-test-rm.r')
################################################################################
# COMMENTS:
# - To test a different metric, Ctrl replace the name of the metric
# There are two equivalent ways of performing a repeated measures ANOVA:
#   1. aov(Score~(SupportLevel*SL_setnum) + Error(wid = .(Subject)/(SupportLevel*SL_setnum)))
#   2. ezANOVA(data,dv=Score,wid=wid = .(Subject),within = .(SupportLevel,SL_setnum),between = NULL, type = 2, detailed = TRUE)
# The second also tests for sphericity, but sometimes gives errors if nothing is significant.
################################################################################

options(contrasts=c("contr.sum","contr.poly"))
remove(list = ls())

# loads packages and dependencies
library(ez) 
library(rstatix) # for the %>% function

# parameters 
DIR = 'C:/Users/milli/OneDrive/Documents/DowntownGames-Biodex/Fall2020-DataAnalysis/FreqAnalysis'

################################################################################
################################################################################
#               Import Data 
################################################################################
################################################################################
data_original = read.csv(paste(DIR,paste("freq","metrics.csv",sep="-"),sep="/"))

################################################################################
################################################################################
#              Metric: Energy at Resonance for
#                      F0_B1, F1_B1_add, and F1_B1
################################################################################
################################################################################

data = subset(data_original, Force!="F0_B0" & Force!="F1_B0" & Force!="F1_B1_add")
data[] <- lapply(data, function(x) if(is.factor(x)) factor(x) else x)

# define file to save data to
sink(paste(DIR,"Stats","Controls",paste("EnergyatResonance","h2.txt",sep="-"),sep="/"))

cat("\n")
cat("########################################################################### \n")
cat("######### Energy at Resonance for F0_B1 and F1_B1 ############# \n")
cat("########################################################################### \n")

# Check assumptions
# Assumption 1: Independence of samples - satisfied by the experimental set-up
# Assumption 2: Data is normally distributed - check visually using boxplots or histograms
#               AND use the shapiro test (if the null hypothesis is violated p<.05, it is not
#               normally distributed
cat("Test for normality \n")
normality = data %>%
  group_by(Force,BallFreq) %>%
  shapiro_test(Resonance)
print(normality)

# Assumption 3: Sphericity/Equal variances between treatments - ezANOVA performs this test
#               if the null hypothesis is violated p<.05 for any particular factor, use a corrected
#               p-value, possibly the Greenhouse-Geisser p[GG]
mod.ez<-ezANOVA(data,Resonance,
                  wid = .(Subject),
                  within = .(Force,BallFreq),
                  between = NULL, type = 2, detailed = TRUE)

print(mod.ez)

# If a factor is significant, that means that at least one of the groups is different from the others
# But you still do not know which group is different. T-tests can help.
posthoc<-pairwise.t.test(data$Resonance,
                         data$Force,
                         paired = TRUE,
                         p.adjust.method = "bonferroni")
print(posthoc)
# # Compare all group/factor combinations for putting asterisks on plots
# data$combo <- paste(data$Force,data$BallFreq)
# posthoc<-pairwise.t.test(data$Resonance,
#                          data$combo,
#                          paired = TRUE,
#                          p.adjust.method = "bonferroni")
# print(posthoc)

cat("\n")
cat("########################################################################### \n")
cat("######### 0.5Hz:      Energy at Resonance for F0_B1 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="0.5Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("########################################################################### \n")
cat("######### 1Hz:      Energy at Resonance for F0_B1 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="1Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("########################################################################### \n")
cat("######### 1.5Hz:      Energy at Resonance for F0_B1 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="1.5Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("########################################################################### \n")
cat("######### 2.5Hz:      Energy at Resonance for F0_B1 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="2.5Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

################################################################################
################################################################################
#              Metric: Energy at Resonance for
#                      F1_B0 and F1_B1
################################################################################
################################################################################

data = subset(data_original, Force=="F1_B0" | Force=="F1_B1")
data[] <- lapply(data, function(x) if(is.factor(x)) factor(x) else x)

# define file to save data to
sink(paste(DIR,"Stats","Controls",paste("EnergyatResonance","h3.txt",sep="-"),sep="/"))

cat("\n")
cat("########################################################################### \n")
cat("################ Energy at Resonance for F1_B0 and F1_B1 ################## \n")
cat("########################################################################### \n")

# Check assumptions
# Assumption 1: Independence of samples - satisfied by the experimental set-up
# Assumption 2: Data is normally distributed - check visually using boxplots or histograms
#               AND use the shapiro test (if the null hypothesis is violated p<.05, it is not
#               normally distributed
cat("Test for normality \n")
normality = data %>%
  group_by(Force,BallFreq) %>%
  shapiro_test(Resonance)
print(normality)

# Assumption 3: Sphericity/Equal variances between treatments - ezANOVA performs this test
#               if the null hypothesis is violated p<.05 for any particular factor, use a corrected
#               p-value, possibly the Greenhouse-Geisser p[GG]
mod.ez<-ezANOVA(data,Resonance,
                wid = .(Subject),
                within = .(Force,BallFreq),
                between = NULL, type = 2, detailed = TRUE)

print(mod.ez)

# If a factor is significant, that means that at least one of the groups is different from the others
# But you still do not know which group is different. T-tests can help.
posthoc<-pairwise.t.test(data$Resonance,
                         data$Force,
                         paired = TRUE,
                         p.adjust.method = "bonferroni")
print(posthoc)
# # Compare all group/factor combinations for putting asterisks on plots
# data$combo <- paste(data$Force,data$BallFreq)
# posthoc<-pairwise.t.test(data$Resonance,
#                          data$combo,
#                          paired = TRUE,
#                          p.adjust.method = "bonferroni")
# print(posthoc)

cat("\n")
cat("########################################################################### \n")
cat("######### 0.5Hz:      Energy at Resonance for F1_B0 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="0.5Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("########################################################################### \n")
cat("######### 1Hz:      Energy at Resonance for F1_B0 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="1Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("########################################################################### \n")
cat("######### 1.5Hz:      Energy at Resonance for F1_B0 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="1.5Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)

cat("\n")
cat("########################################################################### \n")
cat("######### 2.5Hz:      Energy at Resonance for F1_B0 and F1_B1 ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="2.5Hz")
mod.ez<-ezANOVA(data_freq,Resonance,
                wid = .(Subject),
                within = .(Force),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
