# notes ----
# functions for dissecting distributions
# sarah.power@alaska.gov
# 8/15/2018

# load ----
library(tidyverse)
library(mixdist)
library(data.table)
library(lubridate)
library(gridExtra)
library(sfsmisc)
citation("mixdist")
#library(here)
#?ecdf.ksCI()
windowsFonts(Times=windowsFont("Times New Roman"))
options(scipen = 999)

theme_sleek <- function(base_size = 12, base_family = "Times") {
  half_line <- base_size/2
  theme_light(base_size = 12, base_family = "Times") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      panel.border = element_rect(fill = NA),
      legend.key.size = unit(0.9, "lines"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}

#theme_set(theme_sleek())

# functions ----
data_prep <- function(df, year_wanted){
  #df <- df_data
  #year_wanted <- 2006
  df %>% 
    dplyr::filter(year(date) == year_wanted) %>%
    # Create a data frame whose first column are the dates in numeric format
    # and whose second column are the frequencies. 
    # This is required for fitting the mixture. See mixdata {mixdist}
    dplyr::select(day_of_year, run) -> df
}

#This function is for preping early run data as defined by genetics for distribution fitting
data_prep_early <- function(df, year_wanted){
  df %>% 
    dplyr::filter(year(date)==year_wanted) %>%
    # Create a data frame whose first column are the dates in numeric format
    # and whose second column are the frequencies. 
    # This is required for fitting the mixture. See mixdata {mixdist}
    dplyr::select(day_of_year, run_early_gen) -> df
}

year_stats <- function (df, year_wanted){
  #df <- df_data
  #year_wanted <- 2013
  df %>%
    filter(year(date)==year_wanted)-> df
  #run_size <-sum(df$run)
  #print("Run")
  #print(run_size)
  df_fit <- data_prep(df, year_wanted)
  fit <- distribution_estimation_norms_SEQ(df_fit) 
  
  #fit <- (fit <- mix(as.mixdata(df08), mixparam(mu=c(177, 205, 245), sigma= c(10,10,8.6)), constr = mixconstr(conmu="MFX", fixmu= c(TRUE, FALSE, FALSE)), dist= 'gamma', iterlim=5000)) #constr = mixconstr(consigma="SEQ"),
  #dist_plot (fit, year_wanted)
  #dist_plot (fit)
  df_fit_early <- data_prep_early(df, year_wanted)
  fit_early <- distribution_estimation_norms(df_fit_early)
  
  #print("Mean day-of-year & sd based on runtiming distributions alone")
  #print(c(yday(floor(fit$parameters$mu[1])), fit$parameters$sigma[1]))
  #print("Mean day-of-year & sd based on additional genetics information")
  #print(c(yday(fit_early$parameters$mu[1]), fit_early$parameters$sigma[1]))
  
  #Must use group_by(day_of_year) %>% to get correct calculations
  df %>%
    dplyr::group_by(day_of_year) %>% 
    dplyr::mutate(dist_percent = percent_dist(fit,day_of_year),
                  run_early_dis = dist_percent*run,
                  run_early_dis_all = dist_percent*run_all,
                  run_early_gen_all = prop_early_genetics*run_all) -> df
  
  df %>%
    dplyr::group_by(year) %>% 
    dplyr::mutate(cum_run_dis = cumsum(replace_na(run_early_dis, 0)),
                  cum_run_dis_all = cumsum(replace_na(run_early_dis_all, 0)),
                  cum_run_gen = cumsum(replace_na(run_early_gen, 0)),
                  cum_run_gen_all = cumsum(replace_na(run_early_gen_all, 0)),
                  cum_run_all = cumsum(replace_na(run_all,0))) %>%
    dplyr::ungroup(year) -> df #   %>% View(cum_run_gen) #
  
  p_early_dis <- sum(df$run_early_dis)/sum(df$run)
  p_early_gen <- sum(df$run_early_gen)/sum(df$run)
  
  
  min(df$day_of_year, na.rm = TRUE)
  
  # this is for simultaneous graphing. 
  df_long <- df %>%
    gather(model_type, proportions, dist_percent, prop_early_genetics)
  length(df_long$proportions)
  
  
  ks.test(df$dist_percent, df$prop_early_genetics)
  
  #ecdf.ksCI(df$dist_percent)
  
  ggplot(df, aes(day_of_year, prop_early_genetics)) +
    geom_point(size=2) + theme_light()
  ggplot(df, aes(day_of_year, dist_percent)) +
    geom_point(size=2) + theme_light()
  
  log_curve <- ggplot(df_long, aes(x = day_of_year, y = proportions), color = model_type) +
    geom_point(aes(pch = model_type)) + theme_light() +
    theme(legend.justification = c(.5,0), legend.position = "bottom")
  
  ggsave(filename = paste0("figures/log_curv", year_wanted, ".png", sep = ""), device = png(), width = 7, height = 9, units = "in", dpi = 300)
  
  #print("Number of early run by runtiming distribution")
  #print(max(df$cum_run_dis))
  #print("Number of early run by genetics runtiming distribution")
  #print(max(df$cum_run_gen))
  #print("runtiming distribution/genetics runtiming distribution")
  #print(max(df$cum_run_dis)/max(df$cum_run_gen))
  #xaxis <- tickr(df, day_of_year, 5)
  
  #df %>%
  #  dplyr::select(day_of_year, dist_percent, prop_early_genetics) %>% 
  #  melt(id = "day_of_year") -> df2
  #ggplot(df2, aes(day_of_year, value, colour = variable))+
  #  geom_line(size = 3)+
  #  scale_colour_manual(name = "Modeled by",
  #                      labels = c("Distribution only", "Genetics"), 
  #                      values=c("green", "blue"))+
  #  labs(y = "Proportion of run", x= "Day of the year")+
  #  theme(legend.justification = c(1,1), legend.position = c(1,1))+
  #  ggtitle(paste0(year_wanted, " Runtiming Assignment "))
  
  df %>%
    dplyr::select(day_of_year, cum_run_dis, cum_run_gen) %>% 
    melt(id = "day_of_year") -> df3
  unique(df3$variable)
  #yaxis <- tickr(df3, value, 10000)  
  runCDF <- ggplot(df3, aes(day_of_year, value, group = variable)) +
    geom_line(size = 1.5, aes(linetype = variable), show.legend = FALSE) +
    scale_linetype_manual(#name = "Modeled by",
      #labels = c("Distribution only", "Genetics"),
      values = c("solid", "dotted")) +
    labs(y = " ", x= " ") +
    #scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels)+
    #scale_y_continuous(breaks = yaxis$breaks, labels = yaxis$labels)+
    coord_cartesian(xlim = c(150, 220)) +
    ggtitle(paste0(year_wanted)) + #, " Early Run Estimation")) + 
    theme_bw() +
    theme(legend.justification = c(.5,0), legend.position = "bottom")
  my_list <- list(df = df, "logistic" = log_curve, "runCDF" = runCDF)
  return(my_list) 
  
}

graph_year <- function(df){
  #Graph daily weir run = escapement + catch for a year ----
  ggplot(df, aes(day_of_year, run_early)) + 
    geom_line() +
    labs(title = "Early Run Daily weir passage", x = "date in number format", y = "Number of fish")
}

graph_year_early <- function(df){
  #Graph daily weir run = escapement + catch for a year ----
  ggplot(df, aes(day_of_year, run_early)) + 
    geom_line() +
    labs(title = "Early Run Daily weir passage", x = "date in number format", y = "Number of fish")
}

#This function takes the genetically defined early run and using it's mean and standard deviation fits a "gamma"
# (normal) distribution. This steps helps us to see what the genetically defined early run parameters are.
#This give us a basis for our starting points when fitting the distributions without genetics
early_look <- function(df, year_wanted){
  #df <- df_data
  #year_wanted <- 2010
  df <- data_prep_early(df, year_wanted)
  fit <- mix(as.mixdata(df), mixparam(mu=mean(df$day_of_year), sigma=sd(df$day_of_year)), dist="gamma", iterlim=5000)
  dist_plot(fit, year_wanted)
  output  <- cbind(fit$parameters[1], fit$parameters[2], fit$parameters[3],fit$se[1], fit$se[2], fit$se[3])
  return(output)
}

#The following functions starting with distribution_estimation make use of the Expectation Maximization algorithms
#found int he mixdist package. Some amount of coding is in place to "automate" the functions by using the mean and sigmas
#from the data set, but that is only one way to set the starting values, one can also guess at the means and sigmas of the distributions

distribution_estimation_norms <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  
  if(missing(mean_guess_given)) {
    mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  } else {
    mean_guess = mean_guess_given
  }
  if(missing(sigma_guess_given)) {
    sigma_guess = rep(9, each = num_of_distributions)
  } else {
    sigma_guess = sigma_guess_given
  }
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess), dist=distibution_guess, iterlim=5000))  
  
}

#This fits with the additional constraint that means are fit with equal spacing.
distribution_estimation_norms_MES <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  
  if(missing(mean_guess_given)) {
    mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  } else {
    mean_guess = mean_guess_given
  }
  if(missing(sigma_guess_given)) {
    sigma_guess = rep(9, each = num_of_distributions)
  } else {
    sigma_guess = sigma_guess_given
  }
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(conmu="MES"), dist=distibution_guess, iterlim=5000))  
  
}

#This fits with the additional constraint that standard deviations are equal.
distribution_estimation_norms_SEQ <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  #df<-df_fit 
  if(missing(mean_guess_given)) {
    mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  } else {
    mean_guess = mean_guess_given
  }
  if(missing(sigma_guess_given)) {
    sigma_guess = rep(9, each = num_of_distributions)
  } else {
    sigma_guess = sigma_guess_given
  }
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(consigma="SEQ"), dist=distibution_guess, iterlim=5000))  
  
}

#This fits with the additional constraint that the first mean is fixed.
distribution_estimation_norms_SEQ_n1 <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  
  if(missing(mean_guess_given)) {
    mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  } else {
    mean_guess = mean_guess_given
  }
  if(missing(sigma_guess_given)) {
    sigma_guess = rep(9, each = num_of_distributions)
  } else {
    sigma_guess = sigma_guess_given
  }
  
  fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(consigma="SEQ"), dist=distibution_guess, iterlim=5000)
  
  mean_guess <- c(fitpro$parameters$mu[1],fitpro$parameters$mu[2],fitpro$parameters$mu[3])
  fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(conmu= "MFX", fixmu = c(TRUE, FALSE, FALSE)), dist=distibution_guess, iterlim=5000)
  
}

#Fits with the additional constraint that the sigma on the first distribution is fixed
distribution_estimation_norms_SFX <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  
  mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  sigma_guess = c(8.6, 8.6, 8.6)
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(consigma="SFX", fixsigma= c(TRUE, FALSE, FALSE)), dist=distibution_guess, iterlim=5000))  
}

#Fits with the additional constraint that the first mean is fixed, in this case at 174, and sigma is guessed to by 8.6 for each of the 3 distributions.
distribution_estimation_norms_MFX <-  function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  
  mean_guess = c(174, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  sigma_guess = c(8.6, 8.6, 8.6)
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(conmu="MFX", fixmu= c(TRUE, FALSE, FALSE)), dist=distibution_guess, iterlim=5000))  
}

#Fits with the constraint that the Weibull distribution is used with the guess of the sigmas  = 9 
distribution_estimation_weibull <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'weibull'){
  
  if(missing(mean_guess_given)) {
    mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  } else {
    mean_guess = mean_guess_given
  }
  if(missing(sigma_guess_given)) {
    sigma_guess = rep(9, each = num_of_distributions)
  } else {
    sigma_guess = sigma_guess_given
  }
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess), dist=distibution_guess))  #, iterlim=5000
  
}

#Fits with the constraint that the Weibull distribution is used with the sigmas equal. 
distribution_estimation_weibull_SEQ <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'weibull'){
  
  if(missing(mean_guess_given)) {
    mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  } else {
    mean_guess = mean_guess_given
  }
  if(missing(sigma_guess_given)) {
    sigma_guess = rep(9, each = num_of_distributions)
  } else {
    sigma_guess = sigma_guess_given
  }
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess), constr = mixconstr(consigma="SEQ"), dist=distibution_guess))  #, iterlim=5000
  
}

#plots the distributions
dist_plot <- function (fitpro, year_wanted){
  #Plot the results  
  # pdf(NULL)
  # dev.control(displaylist="enable")
  plot(fitpro, main=year_wanted) 
  grid()  
  legend("topright", lty=1, lwd=c(1, 1, 2), c("Original Distribution to be Fit", "Individual Fitted Distributions", "Fitted Distributions Combined"), col=c("blue", "red", rgb(0.2, 0.7, 0.2)), bg="white")  
  # p1 <- recordPlot()
  # invisible(dev.off())
  # return(p1)
  #Estimated mean date and sigmas.
  #summary(fitpro) 
}

auto_year<- function (df, year_wanted) {
  # df <- weir_data
  # year_wanted <- 2013
  
  df <- data_prep(df, year_wanted) 
  #graph_year(df_year)
  #fitpro <- distribution_estimation_norms(df_year) 
  #dist_plot(fitpro, year_wanted )
  #fitpro <- distribution_estimation_norms_MFX(df_year) 
  #dist_plot(fitpro, year_wanted ) 
  fitpro <- distribution_estimation_norms_SEQ(df) 
  dist_plot(fitpro, year_wanted)
}

#The density for the specified distribution (dist_num)
dnormfit <- function(fit, x, dist_num = 1){
  #fit$parameters$pi[dist_num]*dnorm(x, fit$parameters$mu[dist_num],fit$parameters$sigma[dist_num])
  dnorm(x, fit$parameters$mu[dist_num],fit$parameters$sigma[dist_num])
}

#This is the percent for the specified distribution (dist_num) compared to all the other distributions. 
# x = day of the year
#If unspecified the distribution number is 1, representing the first, or black lake distribution.
percent_dist <- function(fit, x){
  dnormfit(fit, x, 1)/sum(dnormfit(fit, x, 1), dnormfit(fit, x, 2), dnormfit(fit, x, 3), na.rm =TRUE) 
}

#This bootstraps the estimated SD for percent of the first distribution on day x, for the "fit" model. 
#this isn't working
percent_dist_se <- function(fit, x, dist_num = 1, sim = 1000){
  m1 <- rnorm(sim, fit$parameters$mu[1], fit$se$mu.se[1])
  m1[is.na(m1)] <- 0
  mu1 <- dnorm(x, m1, fit$parameters$sigma[1])
  dnorm(x, 171, fit$parameters$sigma[1])
  
  m2 <- rnorm(sim, fit$parameters$mu[2], fit$se$mu.se[2])
  m2[is.na(m2)] <- 0
  mu2 <- dnorm(x, m2, fit$parameters$sigma[2])
  
  m3 <- rnorm(sim, fit$parameters$mu[3], fit$se$mu.se[3])
  m3[is.na(m3)] <- 0
  mu3 <- dnorm(x, m3, fit$parameters$sigma[3])
  
  vector <- mu1/sum(mu1, mu2, mu3)
  dnorm(x, mu1, fit$parameters$sigma[1])/sum(dnorm(x, mu1, fit$parameters$sigma[1]), dnorm(x, mu2, fit$parameters$sigma[2]), dnorm(x, mu3, fit$parameters$sigma[3]), na.rm =TRUE) 
  sd(vector)
}

pnormfit <- function(fit, x, dist_num =1){
  fit$parameters$pi[dist_num]*pnorm(x, fit$parameters$mu[dist_num],fit$parameters$sigma[dist_num])
}

weir_date <- function(area, harvest_date){
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
  weir_arrival_date <-
    case_when(
      area == "Unimak" ~ harvest_date + 12,
      area == "Ikatan" ~ harvest_date + 11,
      area == "Dolgoi" ~ harvest_date + 9,
      area == "Shumagins" ~ harvest_date + 7, #This line & above a guestimate based on other values and distances.
      area == "SEDM" ~ harvest_date + 6,
      area == "Perryville" ~ harvest_date + 4,
      area == "Western" ~ harvest_date + 3,
      area == "Chignik" ~ harvest_date + 1,
      area == "Central" ~ harvest_date + 2,
      area == "Eastern" ~ harvest_date + 4,
      area == "Igvak" ~ harvest_date + 6
    )
}
#unique(harvest_data$area)
allocation <- function(area, harvest_date){
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
  if(month(harvest_date) < 7){
    percent_allocated <-
      case_when(
        area == "Unimak" ~ 0,
        area == "Ikatan" ~ 0,
        area == "Dolgoi" ~ .0,
        area == "Shumagins" ~ .0, #This line & above a guestimate based on other values and distances.
        area == "SEDM" ~ .65,
        area == "Perryville" ~ .50,
        area == "Western" ~ .50,
        area == "Chignik" ~ 1,
        area == "Central" ~ .90,
        area == "Eastern" ~ .75,
        area == "Igvak" ~ .75
      )
  }else if(month(harvest_date) == 7){
    percent_allocated <-
      case_when(
        area == "Unimak" ~ 0,
        area == "Ikatan" ~ 0,
        area == "Dolgoi" ~ .0,
        area == "Shumagins" ~ .0, #This line & above a guestimate based on other values and distances.
        area == "SEDM" ~ .55,
        area == "Perryville" ~ .60,
        area == "Western" ~ .60,
        area == "Chignik" ~ 1,
        area == "Central" ~ .95,
        area == "Eastern" ~ .20,
        area == "Igvak" ~ .20
      )
  }else if(month(harvest_date) == 8){
    percent_allocated <-
      case_when(
        area == "Unimak" ~ 0,
        area == "Ikatan" ~ 0,
        area == "Dolgoi" ~ .00,
        area == "Shumagins" ~ .00, #This line & above a guestimate based on other values and distances.
        area == "SEDM" ~ .50,
        area == "Perryville" ~ .50,
        area == "Western" ~ .50,
        area == "Chignik" ~ 1,
        area == "Central" ~ .95,
        area == "Eastern" ~ .20,
        area == "Igvak" ~ .20
      )
  }else if(month(harvest_date) > 8){
    percent_allocated <-
      case_when(
        area == "Unimak" ~ 0,
        area == "Ikatan" ~ 0,
        area == "Dolgoi" ~ .05,
        area == "Shumagins" ~ .05, #This line & above a guestimate based on other values and distances.
        area == "SEDM" ~ .35,
        area == "Perryville" ~ .35,
        area == "Western" ~ .35,
        area == "Chignik" ~ 1,
        area == "Central" ~ .95,
        area == "Eastern" ~ .20,
        area == "Igvak" ~ .20
      )
  }
  return(percent_allocated)
}

#This function still needs work. 
tails_difference <- function(fit, x, dist_a =1, dist_b =2){
  difference <- abs(pnormfit(fit, x, 1) + 1 - pnormfit(fit, x, 2))
}

#This function still needs work. 
tails_equal_date <- function(fit, dist_a =1, dist_b =2){
  # find X (date) when pnorm(dist 1) = 1 - pnorm(dist 2) taking into account the proporiton (pi)
  # each distribution is of the whole multinomial distribution
  
  #find x where
  pnormfit(fitpro, 13338, 2) #is about equal to  
  1 - pnormfit(fitpro, 13338+1, 1) 
  # AKA 
  pnormfit(fitpro, 13338+1, 2) - 1 + pnormfit(fitpro, 13338+1, 1) # is as small as possible. 
  
  # start with x (date) halfway between 2 means
  current_x <- floor((fit$parameters$mu[dist_a]+fit$parameters$mu[dist_b])/2)
  difference_current <- tails_difference(fit, current_x)
  difference_new <- tails_difference(fit, current_x+1)
  
  # if the next day the difference is greater 
  if(difference_current < difference_new ){
    # then consider the day before
    new_x = current_x - 1
    difference_new = tails_difference(fit, new_x)
    while (difference_new < difference_current){ #keep going back a day until the smallest difference is found
      difference_current <- difference_new
      current_x <- new_x
      new_x <- new_x - 1
      difference_new <- tails_difference(fit, new_x)
    } } else{ # difference_new <= difference_current # keep going forward a day until the smallest difference is found
      while (difference_new < difference_current){
        difference_current <- difference_new
        current_x <- new_x
        new_x = new_x + 1
        difference_new = tails_difference(fit, new_x)
      }
    }
  
  return(current_x) #(as.Date.numeric(current_x)) # date 
}