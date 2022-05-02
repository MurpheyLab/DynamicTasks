################################################################################
# Hypothesis 4 determines whether window size significantly impacts energy@resonance 

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
data = read.csv(paste(DIR,paste("freq","metrics","windows.csv",sep="-"),sep="/"))

################################################################################
################################################################################
#              Metric: Energy at Resonance
################################################################################
################################################################################

data[] <- lapply(data, function(x) if(is.factor(x)) factor(x) else x)

# define file to save data to
sink(paste(DIR,"Stats","Controls",paste("EnergyatResonance","windows","h4.txt",sep="-"),sep="/"))

cat("\n")
cat("########################################################################### \n")
cat("######### Energy at Resonance for Different Window Sizes ############# \n")
cat("########################################################################### \n")

# Check assumptions
# Assumption 1: Independence of samples - satisfied by the experimental set-up
# Assumption 2: Data is normally distributed - check visually using boxplots or histograms
#               AND use the shapiro test (if the null hypothesis is violated p<.05, it is not
#               normally distributed
cat("Test for normality \n")
normality = data %>%
  group_by(Window,BallFreq) %>%
  shapiro_test(Energy)
print(normality)

# Assumption 3: Sphericity/Equal variances between treatments - ezANOVA performs this test
#               if the null hypothesis is violated p<.05 for any particular factor, use a corrected
#               p-value, possibly the Greenhouse-Geisser p[GG]
mod.ez<-ezANOVA(data,Energy,
                  wid = .(Subject),
                  within = .(Window,BallFreq),
                  between = NULL, type = 2, detailed = TRUE)

print(mod.ez)

# If a factor is significant, that means that at least one of the groups is different from the others
# But you still do not know which group is different. T-tests can help.
posthoc<-pairwise.t.test(data$Energy,
                         data$Window,
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
cat("######### 0.5Hz:     ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="0.5Hz")
mod.ez<-ezANOVA(data_freq,Energy,
                wid = .(Subject),
                within = .(Window),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
posthoc<-pairwise.t.test(data_freq$Energy,
                         data_freq$Window,
                         paired = TRUE,
                         p.adjust.method = "bonferroni")
print(posthoc)

cat("\n")
cat("########################################################################### \n")
cat("######### 1Hz:      ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="1Hz")
mod.ez<-ezANOVA(data_freq,Energy,
                wid = .(Subject),
                within = .(Window),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
posthoc<-pairwise.t.test(data_freq$Energy,
                         data_freq$Window,
                         paired = TRUE,
                         p.adjust.method = "bonferroni")
print(posthoc)

cat("\n")
cat("########################################################################### \n")
cat("######### 1.5Hz:    ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="1.5Hz")
mod.ez<-ezANOVA(data_freq,Energy,
                wid = .(Subject),
                within = .(Window),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
posthoc<-pairwise.t.test(data_freq$Energy,
                         data_freq$Window,
                         paired = TRUE,
                         p.adjust.method = "bonferroni")
print(posthoc)

cat("\n")
cat("########################################################################### \n")
cat("######### 2.5Hz:      ############# \n")
cat("########################################################################### \n")

data_freq = subset(data, BallFreq=="2.5Hz")
mod.ez<-ezANOVA(data_freq,Energy,
                wid = .(Subject),
                within = .(Window),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
posthoc<-pairwise.t.test(data_freq$Energy,
                         data_freq$Window,
                         paired = TRUE,
                         p.adjust.method = "bonferroni")
print(posthoc)
