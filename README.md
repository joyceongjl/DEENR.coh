# DEENR.coh
Scripts for Rutgers DEENR wavelet coherence workshop
Data files and R code will be uploaded here for the Advanced Ecological Data Analysis course

## Notes about data
1. I normally just use a .csv file with the species abundance (or some other measure) in rows and time in columns. 
2. Make sure its read in as a data matrix.
3. Make sure no missing values (NAs), recommendation in vignette is to replace NAs (if very little) with the median of the non-missing values in the timeseries, since less likely to produce significant synchrony or coherence.
4. Time component should be at regular intervals, e.g. daily, weekly, monthly, yearly, etc.
5. Timeseries should be long enough for wavelets to pick up the timescales of interest, ie. 2x smallest timestep and 1/3 entire timeseries 
6. For this workshop, I'm assuming that the data have been checked beforehand already, I'm not going over data checking steps 
