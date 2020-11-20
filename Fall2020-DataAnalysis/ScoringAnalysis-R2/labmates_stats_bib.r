

options(contrasts=c("contr.sum","contr.poly"))
remove(list = ls())
#library(ez) # loads the ez package and dependencies

col_default = c('steelblue','tomato3')

#DIR = "/media/ola/Data/Research/DowntownData/Fall2020RoundTwo/"
DIR = "E:\\Research\\DowntownData\\Fall2020RoundTwo\\" # windows

data = read.csv(paste(DIR,paste("all","metrics","bib.csv",sep="-"),sep="\\"))
sink(paste(DIR,paste("labmates","comparison","bib.txt",sep="-"),sep="\\"))

type_participants = "labmates"
freq_labels = c("Freq1","Freq2","Freq3","Freq4")

attach(data)

#################################################
cat("\n")
cat("############ Score ############ \n")
# mod.ez<-ezANOVA(data,ModScore,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
# print(mod.ez)
mscoremodel<-aov(Score~(Freq*Condition) + Error(Subject/(Freq*Condition)))
print(summary(mscoremodel))

mscoremodel<-aov(Score~(Dominance*Side) + Error(Subject))
print(summary(mscoremodel))

jpeg(paste(DIR,'Score_forces.jpg'))
    boxplot(Score~Forces,main = type_participants,xlab = "Forces",ylab="Score")
dev.off()
jpeg(paste(DIR,'Score_freqs.jpg'))
    boxplot(Score~Freq,main = type_participants,xlab = "Freqs",ylab="Score")
           #col=col_default,ylim=c(0,40),names=supportlevel_labels)
    #legend("bottomleft",c("controls","stroke"),cex=1.2,fill=col_default)
dev.off()
# jpeg(paste(DIR,'Score_forces&freqs&ball.jpg'))
#     boxplot(Score~Forces*Freq*BallMov,main = type_participants,xlab = "Forces*Freq*BallMov",ylab="Score")
# dev.off()
jpeg(paste(DIR,'Score_forces&freqs&ball.jpg'))
    boxplot(Score~Freq*Condition,main = type_participants,xlab = "Freq*Condition",ylab="Score")
dev.off()

#################################################
cat("\n")
cat("############ Time ball-in-bowl ############ \n")
#mod.ez<-ezANOVA(data,TIB,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
#print(mod.ez)
mscoremodel<-aov(TIB~(Freq*Condition) + Error(Subject/(Freq*Condition)))
print(summary(mscoremodel))

mscoremodel<-aov(TIB~(Dominance*Side) + Error(Subject))
print(summary(mscoremodel))

jpeg(paste(DIR,'TIB_forces.jpg'))
    boxplot(TIB~Forces,main = type_participants,xlab = "Forces",ylab="Time Ball-in-bowl")
dev.off()
jpeg(paste(DIR,'TIB_freqs.jpg'))
    boxplot(TIB~Freq,main = type_participants,xlab = "Freqs",ylab="Time Ball-in-bowl")
           #col=col_default,ylim=c(0.5,1),names=supportlevel_labels)
    #legend("bottomleft",c("controls","stroke"),cex=1.2,fill=col_default)
dev.off()
# jpeg(paste(DIR,'TIB_forces&freqs&ball.jpg'))
#     boxplot(TIB~Forces*Freq*BallMov,main = type_participants,xlab = "Forces*Freq*BallMov",ylab="Time Ball-in-bowl",ylim=c(5,30))
# dev.off()
jpeg(paste(DIR,'TIB_forces&freqs&ball.jpg'))
    boxplot(TIB~Freq*Condition,main = type_participants,xlab = "Freq*Condition",ylab="Time Ball-in-bowl",ylim=c(5,30))
dev.off()


#################################################
cat("\n")
cat("############ Mean Ball Energy ############ \n")
#mod.ez<-ezANOVA(data,TIB,Subject,within = .(SupportLevel,Impairment),between = NULL, detailed = TRUE)
#print(mod.ez)
mscoremodel<-aov(MeanEnergy2~(Freq*Condition) + Error(Subject/(Freq*Condition)))
print(summary(mscoremodel))

mscoremodel<-aov(MeanEnergy2~(Dominance*Side) + Error(Subject))
print(summary(mscoremodel))

jpeg(paste(DIR,'MeanEnergy_forces.jpg'))
    boxplot(MeanEnergy2~Forces,main = type_participants,xlab = "Forces",ylab="Mean Energy")#,ylim=c(0,10))
           #col=col_default,ylim=c(0.5,1),names=supportlevel_labels)
    #legend("bottomleft",c("controls","stroke"),cex=1.2,fill=col_default)
dev.off()
jpeg(paste(DIR,'MeanEnergy_freq.jpg'))
    boxplot(MeanEnergy2~Freq,main = type_participants,xlab = "Freqs",ylab="Mean Energy")#,ylim=c(0,10))
dev.off()
# jpeg(paste(DIR,'MeanEnergy_forces&freqs&ball.jpg'))
#     boxplot(MeanEnergy~Forces*Freq*BallMov,main = type_participants,xlab = "Forces*Freqs*BallMov",ylab="Mean Energy",ylim=c(0,10))
# dev.off()
jpeg(paste(DIR,'MeanEnergy_forces&freqs&ball.jpg'))
    boxplot(MeanEnergy2~Freq*Condition,main = type_participants,xlab = "Freq*Condition",ylab="Mean Energy")#,ylim=c(0,10))
dev.off()
