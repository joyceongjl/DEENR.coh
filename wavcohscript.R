#install packages
install.packages("wsyn")
library("wsyn")

#simulate biological timeseries
#assume 50 timesteps (years) in simulated data
#sp1 and sp2 should be synchronous at 4 years, sp2 and sp3 should be synchronous at 8 year timescales with a lag
times<-1:50
sp1<-sin(2*pi*times/4)+rnorm(50, mean=0, sd=0.5)+1.1
sp2<-sp1 + sin(2*pi*times/8)+rnorm(50, mean=0, sd=0.5)
sp3<-3*sin(2*pi*times/8+2*pi*runif(1))+3.1 +rnorm(50, mean=0, sd=0.5)

dat<-rbind(sp1, sp2, sp3)

#box-cox transformation on data, for functionsin wsyn to run, data need to be de-meaned, detrended (especially if using fourier surrogates - these also require ~normal marginals)
bcdat<-cleandat(dat, times,5)#results in a list, different options for 1-5
clndat<-bcdat$cdat
str(clndat)#check that its in the right format, data matrix
clndat[,1:5]#checking of the first few columns of all 3 rows

#see what the 3 timeseries look like:
plot(times, clndat[1,]/3+1, type="l", xlab="time", ylab="time series index", ylim=c(0,4))
for (i in 2:dim(clndat)[1]){
  lines(times, clndat[i,]/3+i)
}

#simulate 2 environmental timeseries, ev1 with cycles at 4 years and ev2 with cycles at 8 years
ev1<-sin(2*pi*times/4)+rnorm(50, mean=0, sd=0.5)+1.1
ev2<-sin(2*pi*times/8)+rnorm(50, mean=0, sd=0.5)+1.1
ev12<-rbind(ev1,ev2)
clnev<-cleandat(ev12, times,5)$cdat

#plot environmental timeseries
plot(times, clnev[1,], type="l", xlab="time", ylab="env var index", ylim=c(-3,3))
lines(times, clnev[2,], col="blue")

# species wavelet mean field plot, to see if there is overall synchrony among species
spwmf<-wmf(clndat, times=1:50, scale.min=2, scale.max.input = 16)#takes into account both phase and magnitude
plotmag(spwmf)
spwpmf<-wpmf(clndat, times=1:50, scale.min=2, scale.max.input = 16, sigmethod="fft", nrand=1000)#only phase synchrony taken into account
plotmag(spwpmf, colorbar=T,neat=T, title="All 3 species wpmf", sigthresh=0.95)
timeavg<-colMeans(Mod(spwpmf$values), na.rm=T)
plot(spwpmf$timescales, timeavg, type="l", xlab="Timescale", ylab="Power")

#environmental wavelet transforms
wt.ev1<-wt(clnev[1,], times=1:50, scale.min=2, scale.max.input = 16)
plotmag(wt.ev1, title="ev1")
wt.ev2<-wt(clnev[2,], times=1:50, scale.min=2, scale.max.input = 16)
plotmag(wt.ev2, title="ev2")
timeavg.ev2<-colMeans(Mod(wt.ev2$values), na.rm=T)
plot(wt.ev2$timescales, timeavg.ev2, type="l", xlab="Timescale", ylab="Power")

#The various plots above would help us to determine the timescale bands that are likely to be important
#we already know that timescales ~4 and ~8 are important because these are simulated data.

#wavelet coherence between sp 1 and sp2
sp1.2coh<-coh(dat1=clndat[1,], dat2=clndat[2,], times=1:50, norm="powall", 
              sigmethod="fast", nrand=1000, scale.min=2, scale.max.input = 16)
sp1.2coh<-bandtest(sp1.2coh, c(3,5))
sp1.2coh<-bandtest(sp1.2coh, c(6,10))
get_bandp(sp1.2coh)
plotmag(sp1.2coh)#black lines are the 96th and 99th quantiles of coherence of surrogates
#solid red line is actual value of coherence, dashed red line is for significance of coherence
plotphase(sp1.2coh)
#negative phase indicates dat2 is leading dat 1, positive phase indicates dat 1 is leading dat 2.
#result is that as predicted, sp 1 is synchronous (in-phase) with sp2 at around 3-5 year timescales.
#general rules are : 
#if -0.25pi < phase < 0.25pi, this is approximately in-phase
#if -0.75pi <= phase <= -0.25pi, or 0.25pi <= phase <= 0.75pi, this is lagged synchrony
#if phase < -0.75pi or phase > 0.75pi, this is approximately anti-phase. 

#wavelet coherence between sp2 and sp3
sp2.3coh<-coh(dat1=clndat[2,], dat2=clndat[3,], times=1:50, norm="powall", 
              sigmethod="fast", nrand=1000, scale.min=2, scale.max.input = 16)
sp2.3coh<-bandtest(sp2.3coh, c(3,5))
sp2.3coh<-bandtest(sp2.3coh, c(6,10))
get_bandp(sp2.3coh)
plotmag(sp2.3coh)
plotphase(sp2.3coh)
#note that long timescales are a bit harder to detect compared to short timescales, just because there are more short cycles, so sometimes the long timescale coherence may not show up
#assuming that there was an 8 year cycle, the plotphase graph tells us that sp2 was lagging behind sp3, such that they were almost anti-phase, which is what we see in the plot of the timeseries

#wavelet coherence between sp1 and ev1 (4 year cycles)

#wavelet coherence between sp3 and ev2 (8 year cycles)