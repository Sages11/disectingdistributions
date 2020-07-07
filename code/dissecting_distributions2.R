# notes ----
# dissecting distributions2
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
source('code/weighted_logistic_regression_genetics.R')

# data ----
# Chig_data only uses Chignik Bay District harvest
chig_data <- read_csv('data/ChigISGrunappt2006-2017catch.by.district.csv') %>% 
  dplyr::select(-Earlycatch, -Latecatch, -Earlytotal, -Latetotal, -Propotionearly) %>%
  dplyr::rename(
    chignik = Chignik_Hook_Kujulik) %>%
  mutate(esc = Earlyesc + Lateesc,
         catch_chignik = Lagoon + chignik,
         run = esc + catch_chignik, #catch is our estimate of what would have occured at the wier had there been no fishing, based on migration timing studies.
         date = mdy(Date),
         year = year(date),
         day_of_year = yday(date)) #-> chig_data 

year_vector <- c(2006:2008,2010:2017)

#This call on the yearly logistic eq.s from weighted_logistic_regression_genetics.R
gen_calc <- function(df, ayear){
  df <- df %>%
    filter(year == ayear) %>%
    select(day_of_year, run, date, year)
  new_data <- df %>%
    select(day_of_year)
  gen <- augment(eval(as.name(paste0("bl",ayear))), newdata= new_data, type.predict = "response") 
  df <- left_join(df, gen, by = "day_of_year")  %>% 
    rename(prop_early_genetics = .fitted, prop_early_genetics_se = .se.fit) %>% 
    mutate(run_early_gen = prop_early_genetics*run) 
}

gen2006 <- gen_calc(chig_data, 2006) # calculates the genetic percentage per day based on logistic regression
gen2007 <- gen_calc(chig_data, 2007)
gen2008 <- gen_calc(chig_data, 2008)
gen2010 <- gen_calc(chig_data, 2010)
gen2011 <- gen_calc(chig_data, 2011)
gen2012 <- gen_calc(chig_data, 2012)
gen2013 <- gen_calc(chig_data, 2013)
gen2014 <- gen_calc(chig_data, 2014)
gen2015 <- gen_calc(chig_data, 2015)
gen2016 <- gen_calc(chig_data, 2016)
gen2017 <- gen_calc(chig_data, 2017)

#puts all the data back together after figuring out genetics %age for each day.
chig_data <- bind_rows(gen2006, gen2007, gen2008, gen2010, gen2011, gen2012, 
                       gen2013, gen2014, gen2015, gen2016, gen2017)

#Maybe this is not needed, it just allows you look and isn't in the calculation.
early_look <- function(df, year_wanted, mu_guess =175, sigma_guess = 9){
  df <- data_prep_early(df, year_wanted)
  fit <- mix(as.mixdata(df), mixparam(mu = mu_guess, sigma = sigma_guess), dist="gamma", iterlim=5000)
  dist_plot(fit, year_wanted)
}

early_look(chig_data, 2006, 175)#, 8)
early_look(chig_data, 2007, 175)
early_look(chig_data, 2008, 177, 8)
early_look(chig_data, 2010, 177, 9)
early_look(chig_data, 2011, 169, 9)
early_look(chig_data, 2012, 175, 9)
early_look(chig_data, 2013, 174, 9)
early_look(chig_data, 2014, 174, 9)
early_look(chig_data, 2015, 177, 10)
early_look(chig_data, 2016, 174, 9)
early_look(chig_data, 2017, 174, 9)


for(i in year_vector ){
  #auto_year(chig_data, i)
  year_stats(chig_data, i)
}

cdf06 <- year_stats(chig_data, 2006)
cdf07 <- year_stats(chig_data, 2007)
cdf08 <- year_stats(chig_data, 2008)
cdf10 <- year_stats(chig_data, 2010)
cdf11 <- year_stats(chig_data, 2011)
cdf12 <- year_stats(chig_data, 2012)
cdf13 <- year_stats(chig_data, 2013)
cdf14 <- year_stats(chig_data, 2014)
cdf15 <- year_stats(chig_data, 2015)
cdf16 <- year_stats(chig_data, 2016)
cdf17 <- year_stats(chig_data, 2017)
#combines individual plots together
(x <- cowplot::plot_grid(cdf06, cdf07, cdf08, cdf10, cdf11, cdf12, 
                        ncol = 3, align = "hv") )                       
(y <- cowplot::plot_grid(cdf13, cdf14, cdf15, cdf16, cdf17,  
                        ncol = 3, align = "hv"))


pdf06 <- auto_year(chig_data, 2006)
pdf07 <- auto_year(chig_data, 2007)
pdf08 <- auto_year(chig_data, 2008)
pdf10 <- auto_year(chig_data, 2010)
pdf11 <- auto_year(chig_data, 2011)
pdf12 <- auto_year(chig_data, 2012)
pdf13 <- auto_year(chig_data, 2013)
pdf14 <- auto_year(chig_data, 2014)
pdf15 <- auto_year(chig_data, 2015)
pdf16 <- auto_year(chig_data, 2016)
pdf17 <- auto_year(chig_data, 2017)

plot(pdf06)

z <- cowplot::plot_grid(pdf06, pdf07, pdf08, pdf10, pdf11, pdf12, 
                        ncol = 3, align = "hv")                        
w <- cowplot::plot_grid(pdf13, pdf14, pdf15, pdf16, pdf17,  
                        ncol = 3, align = "hv")
z
w
