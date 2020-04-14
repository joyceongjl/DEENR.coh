#script for ANE data
#Q1 code
ane.sp6.ev4<-read.csv("D:/Rutgers_postdoc/wavelet_lecture_workshop/DEENR.coh/ANE_sp6_ev4.csv")
rownames(ane.sp6.ev4)<-ane.sp6.ev4[,1]
ane.sp6.ev4<-ane.sp6.ev4[,-1]
colnames(ane.sp6.ev4)<-seq(1955,2014,1)
str(ane.sp6.ev4)
ane.sp6.ev4<-as.matrix(ane.sp6.ev4)

#Q2 code (one way, many alternatives)
plot(1955:2014, ane.sp6.ev4[1,]/6+1, type="l", xlab="time", ylab="time series index", ylim=c(0,7))
for (i in 2:6){
  lines(1955:2014, ane.sp6.ev4[i,]/6+i)
}

plot(1955:2014, ane.sp6.ev4[7,]/4+1, type="l", xlab="time", ylab="time series index", ylim=c(0,5))
lines(1955:2014, ane.sp6.ev4[8,]/4+2)
lines(1955:2014, ane.sp6.ev4[9,]/4+3)
lines(1955:2014, ane.sp6.ev4[10,]/4+4)