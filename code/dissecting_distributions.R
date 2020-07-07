# notes ----
# dissecting distributions
# sarah.power@alaska.gov
# 8/13/2018

# load ----
library(fpp2)
library(RDS)
library(tidyverse)
library(mixdist)
library(lubridate)
library(zoo) # to convert numeric date back to a number
source('code/distribution_functions.R')

# data ----
# Chig_data only uses Chignik Bay District harvest
chig_data <- read_csv('data/ChigISGrunappt2006-2017catch.by.district.csv') %>% 
  dplyr::select(-Earlycatch, -Latecatch, -Earlytotal, -Latetotal) %>%
  dplyr::rename(#date = Date,
         prop_early_genetics = Propotionearly, 
         early_esc_genetics = Earlyesc, 
         late_esc_genetics = Lateesc,
         #early_catch_genetics = Earlycatch,
         #late_catch_genetics = Latecatch,
         #early_total_genetics = Earlytotal,
         #late_total_genetics = Latetotal,
         chignik = Chignik_Hook_Kujulik) %>%
  mutate(esc = early_esc_genetics + late_esc_genetics,
         #catch = early_catch_genetics + late_catch_genetics,
         catch_chignik = Lagoon, # + chignik + Kumlik,
         run = esc + catch_chignik, #catch is our estimate of what would have occured at the wier had there been no fishing, based on migration timing studies.
         run_early_gen = prop_early_genetics*run,
         date = mdy(Date),
         year = year(date),
         day_of_year = yday(date)) #-> chig_data # convert the Date to its numeric equivalent


# weir_data uses harvest data from farther flung areas
weir_data <- read_csv('data/ChigISGrunappt2006-2017.csv') %>% 
  dplyr::select(-X9) %>%
  dplyr::rename(prop_early_genetics = Propotionearly, 
         early_esc_genetics = Earlyesc, 
         late_esc_genetics = Lateesc,
         early_catch_genetics = Earlycatch,
         late_catch_genetics = Latecatch,
         early_total_genetics = Earlytotal,
         late_total_genetics = Latetotal) %>%
  mutate(esc = early_esc_genetics + late_esc_genetics,
         catch = early_catch_genetics + late_catch_genetics,
         run = esc + catch, #catch is our estimate of what would have occured at the wier had there been no fishing, based on migration timing studies.
         run_early_gen = prop_early_genetics*run,
         date = mdy(Date),
         year = year(date),
         day_of_year = yday(date)) # %>% # convert the Date to its numeric equivalent

year_vector <- c(2006:2008,2010:2017)

fix_date_est <-weir_data %>%
  dplyr::select(run, day_of_year, year) %>%
  dplyr::filter(day_of_year < 186) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(sum_july4 = sum(run))

e06 <- early_look(weir_data, 2006) 
e07 <- early_look(weir_data, 2007) 
e08 <- early_look(weir_data, 2008) 
e10 <- early_look(weir_data, 2010) 
e11 <- early_look(weir_data, 2011) 
e12 <- early_look(weir_data, 2012) 
e13 <- early_look(weir_data, 2013) 
e14 <- early_look(weir_data, 2014) 
e15 <- early_look(weir_data, 2015) 
e16 <- early_look(weir_data, 2016) 
e17 <- early_look(weir_data, 2017) 

p_df <- rbind(e06, e07, e08, e10, e11, e12, e13, e14, e15, e16, e17)
genetics_early_summary <- cbind(year_vector, p_df)

e06 <- early_look(chig_data, 2006) 
e07 <- early_look(chig_data, 2007) 
e08 <- early_look(chig_data, 2008) 
e10 <- early_look(chig_data, 2010) 
e11 <- early_look(chig_data, 2011) 
e12 <- early_look(chig_data, 2012) 
e13 <- early_look(chig_data, 2013) 
e14 <- early_look(chig_data, 2014) 
e15 <- early_look(chig_data, 2015) 
e16 <- early_look(chig_data, 2016) 
e17 <- early_look(chig_data, 2017) 

p_df <- rbind(e06, e07, e08, e10, e11, e12, e13, e14, e15, e16, e17)
genetics_early_summary <- cbind(year_vector, p_df)

for(i in year_vector ){
  distribution_estimation_norms_SEQ_n1(chig_data, i)
}

for(i in year_vector ){
  early_look(chig_data, i)
}

for(i in year_vector ){
  auto_year(chig_data, i)
}
df_year <- data_prep(weir_data, 2007) 

run <- as.ts(df_year$run)
autoplot(run, ts.geom = 'bar', fill = 'grey') +
  ggtitle("runtiming") +
  ylab('sockeye') +
  ylim(0, 200000)

ggplot(chig_data, aes(day_of_year, run)) + geom_line() + 
<<<<<<< HEAD
  ylim(0, 200000) +
=======
  ylim(0, 175000)# +
>>>>>>> b8c303971f3ac98b868f79c4b721ce05d1b06709
  facet_wrap(~year)

ggplot(df_year, aes(day_of_year, run)) + geom_line() +
  ylab("sockeye run per day") +
  xlab("day of year") +
  ggtitle("2007 run timing")

max(weir_data$run)

auto_year(weir_data, 2006)

auto_year(weir_data, 2007)
auto_year(weir_data, 2008)
auto_year(weir_data, 2010)
auto_year(weir_data, 2011)
auto_year(weir_data, 2012)
auto_year(weir_data, 2013)
auto_year(weir_data, 2014)
auto_year(weir_data, 2015)
auto_year(weir_data, 2016)
auto_year(weir_data, 2017)

#July 4th was a historic cut off between the early and later runs. 
#Which is day 185 or 186 if it is a leap year.
yday("2017-07-04")



for(i in year_vector){
  year_stats(chig_data, i)
}

year_stats(weir_data, 2006)
year_stats(weir_data, 2007)
year_stats(weir_data, 2008)
year_stats(weir_data, 2010)
year_stats(weir_data, 2011)
year_stats(weir_data, 2012)
year_stats(weir_data, 2013)
year_stats(weir_data, 2014)
year_stats(weir_data, 2015)
year_stats(weir_data, 2016)
year_stats(weir_data, 2017)


auto_year(chig_data, 2006)
auto_year(chig_data, 2007)
auto_year(chig_data, 2008)
auto_year(chig_data, 2010)
auto_year(chig_data, 2011)
auto_year(chig_data, 2012)
auto_year(chig_data, 2013)
auto_year(chig_data, 2014)
auto_year(chig_data, 2015)
auto_year(chig_data, 2016)
auto_year(chig_data, 2017)

year_stats(chig_data, 2006)
year_stats(chig_data, 2007)
year_stats(chig_data, 2008)
year_stats(chig_data, 2010)
year_stats(chig_data, 2011)
year_stats(chig_data, 2012)
year_stats(chig_data, 2013)
year_stats(chig_data, 2014)
year_stats(chig_data, 2015)
year_stats(chig_data, 2016)
year_stats(chig_data, 2017)



#2008 in detail
df08 <- data_prep(weir_data, 2008)
#write.csv(df08, file = "data/df08.csv")

df08<- as.mixdata(df08)
(fit <- mix(as.mixdata(df08), mixparam(mu=c(175, 205, 240), sigma= c(10,10,10)),constr = mixconstr(consigma="SEQ"), dist= 'gamma', iterlim=5000)) 
(fit<- mix(as.mixdata(df08), mixparam(mu=c(175, 205, 240), sigma= c(10,10,10)), constr = mixconstr(consigma="SEQ"), dist='weibull'))  #, iterlim=5000
mean_guess <- c(fit$parameters$mu[1],fit$parameters$mu[2],fit$parameters$mu[3])
plot(fit)
fit$parameters$mu[1]-1

dist_plot(fit, 2006)
percent_dist(fit, 191) #date when fifty percent is early run. According to distributions: July 10 

year_stats(weir_data, 2006)
df08$dist_percent <- percent_dist(fit, df08$day_of_year)
df08$run_dis <- df08$dist_percent*df08$run
df08$cum <- cumsum(df08$run_dis)


percent_dist(fit, 178) #date when fifty percent is early run. According to distributions? 
pnormfit(fit, 178)


# Delay/Lag from Witteveen and Botz 2007

#Location		                Delay 		Stat Area		              Time Period	

#SEDM		                     5	6	  28115-28190		              65% June, 55% July, 50% Aug., 35% Sept. 	
#SEDM		                     5	6	  28115-28130, 28170-28190		65% June, 55% July, 50% Aug., 35% Sept. 
#Perryville District		     3	4	  27540-27560		              50% June, 60% July, 50% Aug., 35% Sept.
#Western District		         2	3	  27330-27394		              50% June, 60% July, 50% Aug., 35% Sept.
#Outer Chignik Bay/Kujulik	 1	2	  27220-27250		              90% June, 95% rest of Season	
#Cape Kumlik		             2	3	  27262-27264		              90% June, 95% rest of Season
#Chignik Bay District        0	1	  27110		                    100% All Season	
#Weir Count		              -1	0				                        100% All Season
#Eastern District		         3	4   27260, 27270-27296		      75% June, 20% July - September	
#Cape Igvak		               5	6	  26275-26295		              75% June, 20% July - September

# %ages to approximated from GSI
#
#1. Chignik managment harvest (yellow in the spreadsheet)
#2. 90% of Igvak harvest (through 7/25)
#3. and 80% of SEDM harvest through 7/25 (excludes NWSS 7/1-7/25)

#outline
# For each year 2006 - 2008 
#Start by getting weir data.
#For each fishing district starting with closest districts and working further away:
#   Add % of harvest to weir and previously added districts.
#   Estimate early run
#   See how well early run estimation matches up with genetics. 
#If addtional harvest makes early run estimation less similar - add on next district - if it still makes early run estiation less similar to genetics estimation then don't add any more areas. 


#Trying to plot dist_plot(fitpro, year_wanted ) as a ggplot instead of a base plot:

#auto_year<- function (df, year_wanted) {
  df <- chig_data 
  year_wanted <- 2006
  df_year <- data_prep(df, year_wanted) 
  #graph_year(df_year)
  #fitpro <- distribution_estimation_norms(df_year) 
  #dist_plot(fitpro, year_wanted )
  #fitpro <- distribution_estimation_norms_MFX(df_year) 
  #dist_plot(fitpro, year_wanted ) 
  fitpro <- distribution_estimation_norms_SEQ(df_year) 
  dist_plot(fitpro, year_wanted )
#}

fitpro$vmat
fitpro$mixdata
fitpro$se
fitpro$parameters
fitpro$parameters$mu[1]
minx <- min(fitpro$mixdata$day_of_year)
maxx <- max(fitpro$mixdata$day_of_year)
ggplot(data = fitpro$mixdata, aes(fitpro$mixdata$day_of_year))+
  #geom_bar(aes(y = fitpro$mixdata$run))+
  #need to incorporate fitpro$pi which is the percenatge of the whole that distribution represents.
  stat_function(fun =dnorm, n = 101, args = list(mean = fitpro$parameters$mu[1], sd = fitpro$parameters$sigma[1]), col = "green") +
  stat_function(fun =dnorm, n = 101, args = list(mean = fitpro$parameters$mu[2], sd = fitpro$parameters$sigma[2]), col = "green") +
  stat_function(fun =dnorm, n = 101, args = list(mean = fitpro$parameters$mu[3], sd = fitpro$parameters$sigma[3]), col = "green") +
  ylab("")+
  scale_y_continuous(breaks = NULL)


#possible help here:
#https://stackoverflow.com/questions/5688082/overlay-histogram-with-density-curve
# create some data to work with
x = rnorm(1000);

# overlay histogram, empirical density and normal density
p0 = qplot(x, geom = 'blank') +   
  geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_color_manual(name = 'Density', values = c('red', 'blue')) + 
  theme(legend.position = c(0.85, 0.85))

print(p0)
ggplot(data = fitpro$mixdata, aes(day_of_year))+
  geom_density()
geom_freqpoly( aes(run), freq =TRUE)
hist(data = fitpro$mixdata, freq = TRUE)
