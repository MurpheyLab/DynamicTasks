

options(contrasts=c("contr.sum","contr.poly"))
remove(list = ls())
#library(ez) # loads the ez package and dependencies

col_default = c('steelblue','tomato3')

#DIR = "/media/ola/Data/Research/DowntownData/Fall2020LabData/" # linux
DIR = "E:\\Research\\DowntownData\\Fall2020LabData\\" # windows

data = read.csv(paste(DIR,paste("all","metrics","nah.csv",sep="-"),sep="\\"))
sink(paste(DIR,paste("labmates","comparison","nah.txt",sep="-"),sep="\\"))

type_participants = "labmates"

attach(data)

#################################################
cat("\n")
cat("############ Score ############ \n")
# mod.ez<-ezANOVA(data,ModScore,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
# print(mod.ez)
mscoremodel<-aov(Score~(Subject))
print(summary(mscoremodel))

mscoremodel<-aov(Score~(category) + Error(Subject/(category)))
print(summary(mscoremodel))

mscoremodel<-aov(Score~(Dominance*Side) + Error(Subject))
print(summary(mscoremodel))

jpeg(paste(DIR,'Score-category.jpg'))
    boxplot(Score~category,main = type_participants,xlab = "category",ylab="Score")
dev.off()
jpeg(paste(DIR,'Score-subject.jpg'))
    boxplot(Score~Subject,main = type_participants,xlab = "Subject",ylab="Score", ylim=c(0,35))
dev.off()

#################################################
cat("\n")
cat("############ Time-Between Average ############ \n")
#mod.ez<-ezANOVA(data,TIB,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
#print(mod.ez)

mscoremodel<-aov(TBave~(Subject))
print(summary(mscoremodel))

mscoremodel<-aov(TBave~(category) + Error(Subject/(category)))
print(summary(mscoremodel))

mscoremodel<-aov(TBave~(Dominance*Side) + Error(Subject))
print(summary(mscoremodel))

jpeg(paste(DIR,'TBave-category.jpg'))
    boxplot(TBave~category,main = type_participants,xlab = "category",ylab="Time-Between Average")
dev.off()
jpeg(paste(DIR,'TBave-subject.jpg'))
    boxplot(TBave~Subject,main = type_participants,xlab = "Subject",ylab="Time-Between Average", ylim=c(0,1.5))
dev.off()


#################################################
cat("\n")
cat("############ Time-Between Median ############ \n")
#mod.ez<-ezANOVA(data,TIB,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
#print(mod.ez)
mscoremodel<-aov(TBmedian~(Subject))
print(summary(mscoremodel))

mscoremodel<-aov(TBmedian~(category) + Error(Subject/(category)))
print(summary(mscoremodel))

mscoremodel<-aov(TBmedian~(Dominance*Side) + Error(Subject))
print(summary(mscoremodel))

jpeg(paste(DIR,'TBmedian-category.jpg'))
    boxplot(TBmedian~category,main = type_participants,xlab = "category",ylab="Time-Between Median") #ylim=c(0,10)
dev.off()
jpeg(paste(DIR,'TBmedian-subject.jpg'))
    boxplot(TBmedian~Subject,main = type_participants,xlab = "Subject",ylab="Time-Between Median", ylim=c(0,1.5)) #ylim=c(0,10)
dev.off()
