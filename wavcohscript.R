##install packages
install.packages("wsyn")
library("wsyn")

##simulate biological timeseries
##assume 50 timesteps (years) in simulated data
##sp1 and sp2 should be synchronous at 4 years, sp2 and sp3 should be synchronous at 8 year timescales with a lag
times<-1:50
sp1<-sin(2*pi*times/4)+rnorm(50, mean=0, sd=0.5)+1.1
sp2<-sp1 + sin(2*pi*times/8)+rnorm(50, mean=0, sd=0.5)
sp3<-3*sin(2*pi*times/8+2*pi*runif(1))+3.1 +rnorm(50, mean=0, sd=0.5)

dat<-rbind(sp1, sp2, sp3)

##box-cox transformation on data, for functionsin wsyn to run, data need to be de-meaned, detrended (especially if using fourier surrogates - these also require ~normal marginals)
bcdat<-cleandat(dat, times,5)#results in a list, different options for 1-5
clndat<-bcdat$cdat
str(clndat)#check that its in the right format, data matrix
clndat[,1:5]#checking of the first few columns of all 3 rows

##see what the 3 timeseries look like:
plot(times, clndat[1,]/3+1, type="l", xlab="time", ylab="time series index", ylim=c(0,4))
for (i in 2:dim(clndat)[1]){
  lines(times, clndat[i,]/3+i)
}

##simulate 2 environmental timeseries, ev1 with cycles at 4 years and ev2 with cycles at 8 years
ev1<-sin(2*pi*times/4)+rnorm(50, mean=0, sd=0.5)+1.1
ev2<-sin(2*pi*times/8)+rnorm(50, mean=0, sd=0.5)+1.1
ev12<-rbind(ev1,ev2)
clnev<-cleandat(ev12, times,5)$cdat

##plot environmental timeseries
plot(times, clnev[1,], type="l", xlab="time", ylab="env var index", ylim=c(-3,3))
lines(times, clnev[2,], col="blue")

## species wavelet mean field plot, to see if there is overall synchrony among species
spwmf<-wmf(clndat, times=1:50, scale.min=2, scale.max.input = 16)#takes into account both phase and magnitude
plotmag(spwmf)
spwpmf<-wpmf(clndat, times=1:50, scale.min=2, scale.max.input = 16, sigmethod="fft", nrand=1000)#only phase synchrony taken into account
plotmag(spwpmf, colorbar=T,neat=T, title="All 3 species wpmf", sigthresh=0.95)
timeavg<-colMeans(Mod(spwpmf$values), na.rm=T)
plot(spwpmf$timescales, timeavg, type="l", xlab="Timescale", ylab="Power")

##environmental wavelet transforms
wt.ev1<-wt(clnev[1,], times=1:50, scale.min=2, scale.max.input = 16)
plotmag(wt.ev1, title="ev1")
wt.ev2<-wt(clnev[2,], times=1:50, scale.min=2, scale.max.input = 16)
plotmag(wt.ev2, title="ev2")
timeavg.ev2<-colMeans(Mod(wt.ev2$values), na.rm=T)
plot(wt.ev2$timescales, timeavg.ev2, type="l", xlab="Timescale", ylab="Power")

#The various plots above would help us to determine the timescale bands that are likely to be important
#we already know that timescales ~4 and ~8 are important because these are simulated data.

##wavelet coherence between sp 1 and sp2
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
plot(times, clndat[1,], type="l", xlab="time", ylab="Transformed index", ylim=c(-3,3))
lines(times, clndat[2,], col="blue")
#plot of timeseries shows the in-phase relationship between these two species

##wavelet coherence between sp2 and sp3
#sigmethod="fast" is preferred because more surrogates can be used in the same computation time
#for the different significance methods, see the optional reading Sheppard et al. 2017 for fast coherence and the wsyn vignette for details on fft (fourier surrogates) and aaft (amplitude-adjusted fourier surrogates)
sp2.3coh<-coh(dat1=clndat[2,], dat2=clndat[3,], times=1:50, norm="powall", 
              sigmethod="fast", nrand=1000, scale.min=2, scale.max.input = 16)
sp2.3coh<-bandtest(sp2.3coh, c(3,5))
sp2.3coh<-bandtest(sp2.3coh, c(6,10))
get_bandp(sp2.3coh)
plotmag(sp2.3coh)
plotphase(sp2.3coh)
#note that long timescales are a bit harder to detect compared to short timescales, just because there are more short cycles, so sometimes the long timescale coherence may not show up
#assuming that there was an 8 year cycle, the plotphase graph tells us that sp2 was lagging behind sp3, such that they were almost anti-phase, which is what we see in the plot of the timeseries
plot(times, clndat[2,], type="l", xlab="time", ylab="Transformed index", ylim=c(-3,3))
lines(times, clndat[3,], col="blue")
#plot of timeseries shows the anti-phase relationship between these two species


##wavelet coherence between sp1 and ev1 (4 year cycles)
sp1ev1coh<-coh(dat1=clnev[1,], dat2=clndat[1,], times=1:50, norm="powall", 
              sigmethod="fast", nrand=1000, scale.min=2, scale.max.input = 16)
sp1ev1coh<-bandtest(sp1ev1coh, c(3,5))
sp1ev1coh<-bandtest(sp1ev1coh, c(6,10))
get_bandp(sp1ev1coh)# as expected, sig p-value at 3-5yr timescales
plotmag(sp1ev1coh)#dotted red line higher than black lines at timescale=4yrs
plotphase(sp1ev1coh)
#result shows that sp1 and ev1 are synchronous (in-phase) at 4 year timescales, as expected.
plot(times, clnev[1,], type="l", xlab="time", ylab="Transformed index", ylim=c(-3,3))
lines(times, clndat[1,], col="blue")
#plot of timeseries shows the in-phase relationship between sp1 and ev1, which have 4 year cycles


##wavelet coherence between sp3 and ev2 (8 year cycles)
sp3ev2coh<-coh(dat1=clnev[2,], dat2=clndat[3,], times=1:50, norm="powall", 
               sigmethod="fast", nrand=1000, scale.min=2, scale.max.input = 16)
sp3ev2coh<-bandtest(sp3ev2coh, c(3,5))
sp3ev2coh<-bandtest(sp3ev2coh, c(6,10))
get_bandp(sp3ev2coh)# as expected, sig p-value at 6-10yr timescales, mean phase value is closer to anti-phase
plotmag(sp3ev2coh)#dotted red line higher than black lines at timescales 6-10yrs
plotphase(sp3ev2coh)
#result shows that sp3 is leading ev2 at 6-10 year timescales, such that these two timeseries are anti-phase or anti-synchronous
plot(times, clnev[2,], type="l", xlab="time", ylab="Transformed index", ylim=c(-3,3))
lines(times, clndat[3,], col="blue")
#plot of timeseries shows the anti-phase relationship between sp3 and ev2 at ~8 year timescales

##upload some real ANE FAO catch data, ANE_sp6_ev4.csv
#note that these timeseries data have already been optimal box-cox transformed
#csv file consists of transformed catches of 6 species from ANE (region 27) and also 4 environmental variables
#NAOdm = North Atlantic Oscillation from dec-mar, MEIjj = Multivariate ENSO index from july to june of the next year,
#AMOyr = Atlantic Multidecadal Oscillation annual value, CrudeOil = price of annual crude oil

#Q1. Get csv data into suitable format to run analyses in package "wsyn", ie. matrix with proper rownames (labels) and column names (years from 1955-2014)
#note, answer to str(data) should be num [1:10, 1:60] ... ...

#Q2. Do plots of the 6 species timeseries and the 4 env vars separatedly (or together) to see what the data looks like
#note, there are multiple ways to plot this

#Q3. Do wavelet mean field and wavelet phasor mean field plots for the 6 species