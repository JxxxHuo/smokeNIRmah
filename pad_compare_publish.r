## MSC smoke
setwd("G:\\new420_smoke")

las_msc<-read.table('SAS1709-1_SmoothingSGolay_MSC_SGolay.txt',header = TRUE, sep=",", dec = ".")
lty_msc<-read.table('STYX1709-1_SmoothingSGolay_MSC_SGolay.txt',header = TRUE, sep=",", dec = ".")
#llt<-read.table('SLT1801S-1_SmoothingSGolay_MSC_SGolay.txt',header = TRUE, sep=",", dec = ".")

## MSC smoke old
setwd("G:\\new420_smoke")

las_msc_old<-read.table('ASsmoke_old_SmoothingSGolay_MSC_SGolay.txt',header = TRUE, sep=",", dec = ".")
lty_msc_old<-read.table('TYsmoke_old_SmoothingSGolay_MSC_SGolay.txt',header = TRUE, sep=",", dec = ".")
pad<-read.table('emptyfilterpad_SmoothingSGolay_MSC_SGolay.txt',header = TRUE, sep=",", dec = ".")

install.packages("mdatools")
library(mdatools)

NIRdata<-rbind(las_msc,lty_msc)
label<-c(rep('AS',nrow(las_msc)),rep('TY',nrow(lty_msc)))
faclabel<-factor(label)
#---------calculate the result by group-------------
#labelnas<-unlist(la.dataly(las_msc$X.1,substring,0,7))
#labelnty<-unlist(la.dataly(lty_msc$X.1,substring,0,8))
#NIRdata$Xn<-c(labelnas,labelnty)

#Xn<-c(labelnas,labelnty)
#data<-unname(as.matrix(NIRdata[,3:ncol(NIRdata)]))
#Nframe<-tibble(month=Xn,data=(data))
#md<-aggregate(Nframe$data,list(Nframe$month),mean)
#md$Group.1

seqB<-grepl('B',NIRdata$`X.1`)&grepl('TY',NIRdata$`X.1`)
dataB<-NIRdata[seqB,3:ncol(NIRdata)]
source("G:\\maha3.r")

### This is an important step which gets the mahalanobis distance between NIRS ###
mah<-maha3(NIRdata[,3:ncol(NIRdata)],colMeans(dataB),cov(dataB))
### maha3 is a modified function which is different from old traditional mah

co<-sd(mah)/mean(mah)

plot(mah,type='b',col='blue',xlab='sample index',ylab='Mahalanobis Distance')
lines((sum(grepl('AS',label))+1):length(mah),mah[which(grepl('TY',label))],type='b',col='red')
lines(tmp2$mah,type='b',col='green')
legend(100,30000,c("AS(3.3e4,5.2e4)"),text.col=c('blue'),bty="n")
legend(500,10000,c("TY(42,3.5e3)"),text.col=c('red'),bty="n")
legend(100,5000,c("paired NIR(0,207)"),text.col=c('green'),bty="n")
range(mah[which(grepl('TY',label))])  # TY range
range(mah[which(grepl('AS',label))])  # AS range
range(na.omit(tmp2$mah))

plot((sum(grepl('AS',label))+1):length(mah),mah[which(grepl('TY',label))],type='b',col='red',
     xlab='sample index',ylab='Mahalanobis Distance',pch=1)
lines((sum(grepl('AS',label))+1):length(mah),tmp2$mah[which(grepl('TY',label))],type='b',
      col='green',pch=2)
legend(500,2500,c("TY"),text.col=c('red'),pch=1,bty='n',col='red')
legend(500,2300,c("paired NIR"),text.col=c('green'),pch=2,bty='n',col='green')

s1<-which(grepl('AS1709',las_msc_old$`X.1`))
s2<-which(grepl('TY1709',lty_msc_old$`X.1`))
NIRdata2<-rbind(las_msc_old[s1[1]:nrow(las_msc_old),],lty_msc_old[s2[1]:nrow(lty_msc_old),])
label2<-c(rep('AS',nrow(las_msc_old[s1[1]:nrow(las_msc_old),])),rep('TY',nrow(lty_msc_old[s2[1]:nrow(lty_msc_old),])))
faclabel2<-factor(label2)

s1<-which(grepl('AS1709',las_msc_old$`X.1`))
s2<-which(grepl('TY1709',lty_msc_old$`X.1`))
NIRdata2<-rbind(las_msc_old[s1[1]:nrow(las_msc_old),],lty_msc_old[s2[1]:nrow(lty_msc_old),])
seqB2<-grepl('B',NIRdata2$`X.1`)&grepl('TY',NIRdata2$`X.1`)
dataB2<-NIRdata2[seqB2,3:ncol(NIRdata2)]
source("G:\\maha3.r")
mah2<-maha3(NIRdata2[,3:ncol(NIRdata2)],colMeans(dataB2),cov(dataB2))
co2<-sd(mah2)/mean(mah2)

TYmah2<-mah2[grepl("TY",label2)]
TYmah<-mah[grepl("TY",label)]
plot(TYmah,type="b")
lines(TYmah2,col='red',type="b")
cco<-sd(TYmah)/mean(TYmah)
cco2<-sd(TYmah2)/mean(TYmah2)

library(rmatio)
setwd("G:\\")
tmp2<-read.mat('pad_compare.mat')

plot(tmp2$mah,type='b',col='green',xlab='sample index',ylab='Mahalanobis Distance')
legend(200,150,c("AS"),text.col=c('blue'),bty="n")
legend(600,60,c("TY"),text.col=c('red'),bty="n")
range(na.omit(tmp2$mah))

source("G:\\maha3.r")
mahtyty<-rep('NA',nrow(NIRdata) )
matched<-list()


library(ggplot2)
library(dplyr)
library(hrbrthemes)

data<-data.frame(type=c(label[which(grepl('TY',label))], rep('Paired NIR',length(which(grepl('TY',label)))),
                        label[which(grepl('AS',label))]),
                  value=c(mah[which(grepl('TY',label))],tmp2$mah[which(grepl('TY',label))],
                          mah[which(grepl('AS',label))]))

data%>% ggplot(aes(x=value,fill=type))+geom_histogram(binwidth=1000,alpha=0.5,color="#e9ecef",position='identity')+
  scale_fill_manual(values=c("red","green","blue"))+theme_ipsum()+labs(fill="")



data<-data.frame(type=c(label[which(grepl('TY',label))], rep('Paired NIR',length(which(grepl('TY',label))))),
                 value=c(mah[which(grepl('TY',label))],tmp2$mah[which(grepl('TY',label))]))

data%>% ggplot(aes(x=value,fill=type))+geom_histogram(binwidth=100,alpha=0.5,color="#e9ecef",position='identity')+
  scale_fill_manual(values=c("red","green"))+theme_ipsum()+labs(fill="")



##########################smv TEST###############################################
