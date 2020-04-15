###script for ANE data
library(wsyn)
library(corrplot)#to plot correlation matrix

##Q1 code
ane.sp6.ev4<-read.csv("D:/Rutgers_postdoc/wavelet_lecture_workshop/DEENR.coh/ANE_sp6_ev4.csv")
rownames(ane.sp6.ev4)<-ane.sp6.ev4[,1]
ane.sp6.ev4<-ane.sp6.ev4[,-1]
colnames(ane.sp6.ev4)<-seq(1955,2014,1)
str(ane.sp6.ev4)
ane.sp6.ev4<-as.matrix(ane.sp6.ev4)

##Q2 code (one way, many alternatives)
plot(1955:2014, ane.sp6.ev4[1,]/6+1, type="l", xlab="time", ylab="time series index", ylim=c(0,7))
for (i in 2:6){
  lines(1955:2014, ane.sp6.ev4[i,]/6+i)
}

plot(1955:2014, ane.sp6.ev4[7,]/4+1, type="l", xlab="time", ylab="time series index", ylim=c(0,5))
lines(1955:2014, ane.sp6.ev4[8,]/4+2)
lines(1955:2014, ane.sp6.ev4[9,]/4+3)
lines(1955:2014, ane.sp6.ev4[10,]/4+4)

##Q3 code (wmf and wpmf plots)
anesp6wmf<-wmf(ane.sp6.ev4[1:6,], times=1955:2014, scale.min=2, scale.max.input = 20)
plotmag(anesp6wmf)
#seems like these 6 species were only synchronous at 2-8 yr timescales from 1970s-1980s
anesp6wpmf<-wpmf(ane.sp6.ev4[1:6,], times=1955:2014, scale.min=2, scale.max.input = 20, sigmethod="fft", nrand=1000)
plotmag(anesp6wpmf, colorbar=T,neat=T, title="ANE 6 species wpmf", sigthresh=0.95)
#strong phase synchrony around 3-8 yr timescales from 60s to mid 80s.
sp6timeavg<-colMeans(Mod(anesp6wpmf$values), na.rm=T)
plot(anesp6wpmf$timescales, sp6timeavg, type="l", xlab="Timescale", ylab="Power")
#5-8 yr timescales important
#I might try to do timescale bands of 2-4, 4-8, 8-16
#below is a plot of winter NAO to see the dominant timescales
wt.NAO<-wt(ane.sp6.ev4[7,], times=1955:2014, scale.min=2, scale.max.input = 20)
plotmag(wt.NAO, title="winter NAO")
timeavg.nao<-colMeans(Mod(wt.NAO$values), na.rm=T)
plot(wt.NAO$timescales, timeavg.nao, type="l", xlab="Timescale", ylab="Power")
#strong 8 yr timescale, also 4-6 year timescale

##Q4 code (wavelet coherence matrices)
anecoh2.4<-synmat(ane.sp6.ev4, times=1955:2014, method="coh", scale.max.input = 16, tsrange=c(2,4))# coherence value
rownames(ane.sp6.ev4)#full names of species and env vars
rownames(anecoh2.4)<-c("T.tra", "P.spp", "B.boo", "D.mac", "P.ery", "S.aur", "NAO", "MEI", "AMO", "Oil")
colnames(anecoh2.4)<-c("T.tra", "P.spp", "B.boo", "D.mac", "P.ery", "S.aur", "NAO", "MEI", "AMO", "Oil")
anecohsig2.4<-synmat(ane.sp6.ev4, times=1955:2014, method="coh.sig.fast", scale.max.input = 16, tsrange=c(2,4), nsurrogs=1000) #significance of coherence based on 1000 surrogates
anecoh2.4pvals<-1-anecohsig2.4 #because the default values for synmat gives 1-p-values
colbwr<-colorRampPalette(c("blue", "white", "red"))#to specify colour palette
corrplot(anecoh2.4, method="number", type="lower", tl.pos="ld", tl.srt=0, tl.offset=0.7, 
         col=colbwr(10), is.corr=TRUE, diag=F, tl.col="black", p.mat=anecoh2.4pvals, 
         sig.level=0.01, insig="blank")
mtext("Coherence at timescale 2-4 years (p<0.01)", side=3, line=3)
#5 of the 6 species have strong coherences with each other at very short timescales of 2-4 years.
#for further species analyses, choose B. boops and P. erythrinus
#MEI and crude oil seem to have strong coherences with Penaeus spp (shrimps) and dentex at 2-4yr timescales
#for further env analyses, choose crude oil and shrimps

anecoh4.8<-synmat(ane.sp6.ev4, times=1955:2014, method="coh", scale.max.input = 16, tsrange=c(4,8))
rownames(anecoh4.8)<-c("T.tra", "P.spp", "B.boo", "D.mac", "P.ery", "S.aur", "NAO", "MEI", "AMO", "Oil")
colnames(anecoh4.8)<-c("T.tra", "P.spp", "B.boo", "D.mac", "P.ery", "S.aur", "NAO", "MEI", "AMO", "Oil")
anecohsig4.8<-synmat(ane.sp6.ev4, times=1955:2014, method="coh.sig.fast", scale.max.input = 16, 
                    tsrange=c(4,8), nsurrogs=1000)
anecoh4.8pvals<-1-anecohsig4.8 
corrplot(anecoh4.8, method="number", type="lower", tl.pos="ld", tl.srt=0, tl.offset=0.7, 
         col=colbwr(10), is.corr=TRUE, diag=F, tl.col="black", p.mat=anecoh4.8pvals, 
         sig.level=0.01, insig="blank")
mtext("Coherence at timescale 4-8 years (p<0.01)", side=3, line=3)
#all 6 species is strongly coherent with another species at 4-8 yr timescales
#for further env analyses try crude oil and dentex.

anecoh8.16<-synmat(ane.sp6.ev4, times=1955:2014, method="coh", scale.max.input = 16, tsrange=c(8,16))
rownames(anecoh8.16)<-c("T.tra", "P.spp", "B.boo", "D.mac", "P.ery", "S.aur", "NAO", "MEI", "AMO", "Oil")
colnames(anecoh8.16)<-c("T.tra", "P.spp", "B.boo", "D.mac", "P.ery", "S.aur", "NAO", "MEI", "AMO", "Oil")
anecohsig8.16<-synmat(ane.sp6.ev4, times=1955:2014, method="coh.sig.fast", scale.max.input = 16, 
                     tsrange=c(8,16), nsurrogs=1000)
anecoh8.16pvals<-1-anecohsig8.16 
corrplot(anecoh8.16, method="number", type="lower", tl.pos="ld", tl.srt=0, tl.offset=0.7, 
         col=colbwr(10), is.corr=TRUE, diag=F, tl.col="black", p.mat=anecoh8.16pvals, 
         sig.level=0.01, insig="blank")
mtext("Coherence at timescale 8-16 years (p<0.01)", side=3, line=3)
#only P.ery and T.tra, and P.ery and shrimps are coherent at 8-16yr timescales.
#for further species analyses, try P.ery and T.tra. 

##Q5 code (specific examples)
#for species, choose P.ery (row 5) and P.spp (row 2)
#for env, choose cruide oil (row 10) and dentex (row 4)
