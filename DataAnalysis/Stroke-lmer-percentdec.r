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

set = 'combined' # 'all' or 'combined'


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
data = read.csv(paste(DIR,paste("stroke-freq-percent-loss-aggregate",".csv",sep=""),sep="/"))

data_all = subset(data, BallFreq!="0.5Hz")
data_modsev = subset(data_all, FMA<40)

sink(paste(DIR,"Stats",paste("percent-tests-mod-sev.txt",sep="-"),sep="/"))

cat("\n")
cat("################################################################################ \n")
cat("################            Percent Decrease Metrics            ################ \n")
cat("################################################################################ \n")

cat("\n\n\n")
cat("########################################################################### \n")
cat("Percent decrease in paretic arm with and without loading combined\n")
mod.ez<-ezANOVA(data_modsev,A_com,
                wid = .(Subject),
                within = .(BallFreq),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
posthoc<-pairwise.t.test(data_modsev$A_com,
                         data_modsev$BallFreq,
                         paired = TRUE,
                         p.adjust.method = "bonferroni",detailed = TRUE)
print(posthoc)

cat("\n\n\n")
cat("########################################################################### \n")
cat("Percent decrease in paretic arm with loading \n")
mod.ez<-ezANOVA(data_modsev,SL1,
                wid = .(Subject),
                within = .(BallFreq),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
posthoc<-pairwise.t.test(data_modsev$SL1,
                         data_modsev$BallFreq,
                         paired = FALSE,
                         p.adjust.method = "bonferroni",detailed = TRUE)
print(posthoc)

cat("\n\n\n")
cat("########################################################################### \n")
cat("Percent decrease in paretic arm with no loading \n")
mod.ez<-ezANOVA(data_modsev,SL0,
                wid = .(Subject),
                within = .(BallFreq),
                between = NULL, type = 2, detailed = TRUE)
print(mod.ez)
posthoc<-pairwise.t.test(data_modsev$SL0,
                         data_modsev$BallFreq,
                         paired = FALSE,
                         p.adjust.method = "bonferroni",detailed = TRUE)
print(posthoc)

