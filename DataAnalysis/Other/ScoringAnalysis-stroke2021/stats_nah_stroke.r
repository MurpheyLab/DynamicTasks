

options(contrasts=c("contr.sum","contr.poly"))
remove(list = ls())
#library(ez) # loads the ez package and dependencies

col_default = c('steelblue','tomato3')

#DIR = "/media/ola/Data/Research/DowntownData/Fall2020LabData/" # linux
DIR = "E:\\Research\\DowntownWork\\DowntownData\\Stroke2021\\" # windows

data = read.csv(paste(DIR,paste("all","metrics","nah.csv",sep="-"),sep="\\"))
sink(paste(DIR,paste("stroke","comparison","nah.txt",sep="-"),sep="\\"))

type_participants = "N-a-H stroke"

attach(data)

#################################################
cat("\n")
cat("############ Score ############ \n")
# mod.ez<-ezANOVA(data,ModScore,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
# print(mod.ez)
mscoremodel<-aov(Score~(Subject))
print(summary(mscoremodel))

mscoremodel<-aov(Score~(Arm*Suplev) + Error(Subject/(Arm*Suplev)))
print(summary(mscoremodel))

# mscoremodel<-aov(Score~(Dominance*Side) + Error(Subject))
# print(summary(mscoremodel))

jpeg(paste(DIR,'NaH-score-arm-suplev.jpg'))
    boxplot(Score~Arm*Suplev,main = type_participants,xlab = "Arm*SupLev",ylab="Score")
dev.off()
jpeg(paste(DIR,'NaH-score-arm.jpg'))
    boxplot(Score~Arm,main = type_participants,xlab = "Arm",ylab="Score", ylim=c(0,30))
dev.off()
# jpeg(paste(DIR,'NaH-score-subject.jpg'))
#     boxplot(Score~Subject,main = type_participants,xlab = "Subject",ylab="Score", ylim=c(0,35))
# dev.off()

# #################################################
# cat("\n")
# cat("############ Time-Between Average ############ \n")
# #mod.ez<-ezANOVA(data,TIB,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
# #print(mod.ez)
#
# mscoremodel<-aov(TBave~(Subject))
# print(summary(mscoremodel))
#
# mscoremodel<-aov(TBave~(arm) + Error(Subject/(arm)))
# print(summary(mscoremodel))
#
# mscoremodel<-aov(TBave~(Dominance*Side) + Error(Subject))
# print(summary(mscoremodel))
#
# jpeg(paste(DIR,'TBave-arm.jpg'))
#     boxplot(TBave~arm,main = type_participants,xlab = "arm",ylab="Time-Between Average")
# dev.off()
# jpeg(paste(DIR,'TBave-subject.jpg'))
#     boxplot(TBave~Subject,main = type_participants,xlab = "Subject",ylab="Time-Between Average", ylim=c(0,1.5))
# dev.off()
#
#
# #################################################
# cat("\n")
# cat("############ Time-Between Median ############ \n")
# #mod.ez<-ezANOVA(data,TIB,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
# #print(mod.ez)
# mscoremodel<-aov(TBmedian~(Subject))
# print(summary(mscoremodel))
#
# mscoremodel<-aov(TBmedian~(arm) + Error(Subject/(arm)))
# print(summary(mscoremodel))
#
# mscoremodel<-aov(TBmedian~(Dominance*Side) + Error(Subject))
# print(summary(mscoremodel))
#
# jpeg(paste(DIR,'TBmedian-arm.jpg'))
#     boxplot(TBmedian~arm,main = type_participants,xlab = "arm",ylab="Time-Between Median") #ylim=c(0,10)
# dev.off()
# jpeg(paste(DIR,'TBmedian-subject.jpg'))
#     boxplot(TBmedian~Subject,main = type_participants,xlab = "Subject",ylab="Time-Between Median", ylim=c(0,1.5)) #ylim=c(0,10)
# dev.off()
