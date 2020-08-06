# notes ----
# functions for dissecting distributions
# sarah.power@alaska.gov
# 7/7/2020


# libraries ----

if(!require("fpp2"))   install.packages("fpp2") 
if(!require("RDS"))   install.packages("RDS") 
if(!require("tidyverse"))   install.packages("tidyverse") # dplyr, ggplot, etc.
if(!require("mixdist"))   install.packages("mixdist") 
if(!require("lubridate"))   install.packages("lubridate") 
if(!require("grid"))   install.packages("grid") 
if(!require("data.table"))   install.packages("data.table") 
if(!require("lubridate"))   install.packages("lubridate") 
if(!require("gridExtra"))   install.packages("gridExtra") 
if(!require("cowplot"))   install.packages("cowplot") 
if(!require("zoo"))   install.packages("zoo")  # to convert numeric date back to a number can conflict with lubridate.

citation("mixdist")
#library(here)

# display settings
windowsFonts(Times=windowsFont("Times New Roman"))
windowsFonts("helvetica" = windowsFont("helvetica"))
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
#theme_set(theme_cowplot())

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# functions ----


# Create a data frame whose first column are the dates in numeric format
# and whose second column are the frequencies. 
# This is required for fitting the mixture. See mixdata {mixdist}
data_prep <- function(df, year_wanted){
  df %>% 
    filter(year(date)==year_wanted) %>%
    dplyr::select(day_of_year, run) -> df
}

# Create a data frame whose first column are the dates in numeric format
# and whose second column are the frequencies. 
# This is required for fitting the mixture. See mixdata {mixdist}
#This using genetics methods to find the early run component. 
#This was helpful when initially figuring out early. 
data_prep_early <- function(df, year_wanted){
  df %>% 
    filter(year(date)==year_wanted) %>%
    dplyr::select(day_of_year, run_early_gen) -> df
}

#Input data frame of runtiming and genetic testing
#output graphs and tables comparing two methods for each year
year_stats <- function (df, year_wanted){
  #df <- df_data
  #year_wanted <- 2007
  df %>%
    filter(year(date)==year_wanted)-> df
  #run_size <-sum(df$run)
  df_fit <- data_prep(df, year_wanted)
  fit <- distribution_estimation_norms_SEQ(df_fit) 
  #dist_plot (fit, year_wanted)
  df_fit_early <- data_prep_early(df, year_wanted)
  fit_early <- distribution_estimation_norms(df_fit_early)
  
  #print("Mean day-of-year & sd based on runtiming distributions alone")
  #print(c(yday(fit$parameters$mu[1]), fit$parameters$sigma[1]))
  #print("Mean day-of-year & sd based on additional genetics information")
  #print(c(yday(fit_early$parameters$mu[1]), fit_early$parameters$sigma[1]))
  
  #this first graph is not displayed, but can be used in years with genetics 
  # to find starting values for the fit for distribution_estimation_norms_SEQ(df_fit) 
  
  #Must use group_by(day_of_year) %>% to get correct calculations
  
  df <- df %>%
    dplyr::group_by(day_of_year) %>% 
    dplyr::mutate(dist_percent = percent_dist(fit,day_of_year),
                  run_early_dis = dist_percent*run,
                  run_early_dis_all = dist_percent*run_all,
                  run_early_gen_all = prop_early_genetics*run_all)

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
  
  
  
  # this "long" data set needed for simultaneous graphing of both methods
  df_long <- df %>%
    gather(model_type, proportions, dist_percent, prop_early_genetics)
  length(df_long$proportions)

  #ks.test(df$dist_percent, df$prop_early_genetics)
  #ecdf.ksCI(df$dist_percent)
  
  #this is to plot vertical lines in plots 
  jun1 <- as.integer(yday(as.Date(paste0(year_wanted, "-06-01"))))
  jul1 <- as.integer(yday(as.Date(paste0(year_wanted, "-07-01"))))
  aug1 <- as.integer(yday(as.Date(paste0(year_wanted, "-08-01"))))
  d <- data.frame(day_of_year = c(jun1, jul1, aug1), 
                  event = c("June 1", "July 1", "August 1"),
                  ltype = c("dotted","dotdash","dotted"))
  
  #Plot just early proportions determined with genetics (not saved)
  ggplot(df, aes(day_of_year, prop_early_genetics)) +
    geom_point(size=2) + 
    ggtitle(label = paste0(year_wanted)) +
    geom_vline(data = d, mapping = aes(xintercept = day_of_year, linetype = ltype), show.legend = FALSE) +
    geom_text(data=d, mapping=aes(x=day_of_year, y=0, label= event), size=4, angle=90, vjust= -0.4, hjust = -0.5) +
    scale_shape_manual(name = "Modeled by",
                       labels = c("Run Timing", "Genetics"), 
                       values = c(19, 17)) +
    theme(plot.title = element_text(vjust = -8, hjust = 1), legend.justification = c(.5,0), legend.position = "bottom") +
    labs(y = "Proportion of early sockeye determined with genetics", x= "Day of the year") +
    coord_cartesian(xlim = c(150, 220))
  #dev.off()
  
  #Plot just early proportions determined with run timing (not saved)
  ggplot(df, aes(day_of_year, dist_percent)) +
    geom_point(size=2) + 
    ggtitle(label = paste0(year_wanted)) +
    geom_vline(data = d, mapping = aes(xintercept = day_of_year, linetype = ltype), show.legend = FALSE) +
    geom_text(data=d, mapping=aes(x=day_of_year, y=0, label= event), size=4, angle=90, vjust= -0.4, hjust = -0.5) +
    scale_shape_manual(name = "Modeled by",
                       labels = c("Run Timing", "Genetics"), 
                       values = c(19, 17)) +
    theme(plot.title = element_text(vjust = -8, hjust = 1), legend.justification = c(.5,0), legend.position = "bottom") +
    labs(y = "Proportion of early sockeye determined with run timing", x= "Day of the year") +
    coord_cartesian(xlim = c(150, 220))
  #dev.off()
  
  #Plot just early proportions both methods (saved)
  early_prop <- ggplot(df_long, aes(x = day_of_year, y = proportions), color = model_type) +
    geom_point(aes(pch = model_type)) +
    ggtitle(label = paste0(year_wanted, " Proportions of Early Run")) +
    #theme_Publication() +
    #theme(legend.justification = c(.5,0), legend.position = "bottom") +
    geom_vline(data = d, mapping = aes(xintercept = day_of_year, linetype = ltype), show.legend = FALSE) +
    geom_text(data=d, mapping=aes(x=day_of_year, y=0, label= event), size=6, angle=90, vjust= -0.4, hjust = -0.5) +
    scale_shape_manual(name = "Modeled by",
                       labels = c("Run Timing", "Genetics"), 
                       values = c(19, 17)) +
    theme(text = element_text(size =16), plot.title = element_text(vjust = -8, hjust = 1, size = 20), legend.justification = c(.5,0), legend.position = "bottom") +
    labs(y = "Proportion of early sockeye run", x= "Day of the year") +
    coord_cartesian(xlim = c(150, 220))
  #early_prop
  ggsave(filename = paste0("figures/early_prop", year_wanted, ".png", sep = ""), device = png(), width = 7, height = 9, units = "in", dpi = 300)
  save_plot(filename = paste0("figures/early_prop", year_wanted, ".png", sep = ""), early_prop, ncol = 1, nrow = 1, base_height =3.71, base_asp = 1.618)
  
  #dev.off()

  #df %>%
  #  mutate(dist_percent = percent_dist(fit, df$day_of_year),
  #         run_early_dis = dist_percent*run,
  #         cum_run_dis = cumsum(run_early_dis),
  #         cum_run_gen = cumsum(run_early_gen),
  #         cum_run_all = cumsum(replace_na(run_all,0))) ->df #cum_run_all
  #min(df$day_of_year, na.rm = TRUE)
  #print("Number of early run by runtiming distribution")
  #print(max(df$cum_run_dis))
  #print("Number of early run by genetics runtiming distribution")
  #print(max(df$cum_run_gen))
  #print("runtiming distribution/genetics runtiming distribution")
  #print(max(df$cum_run_dis)/max(df$cum_run_gen))
  #xaxis <- tickr(df, day_of_year, 5)
  
  #For ploting the cumulative run via both methods 
  df_long_cum <- df %>%
    gather(model_type, cumulative_run, cum_run_dis, cum_run_gen)
  length(df_long_cum$cumulative_run)
  
  #Plot of cumulative nun via both methods
  early_cum <- ggplot(df_long_cum, aes(x = day_of_year, y = cumulative_run), color = model_type) +
    geom_point(aes(pch = model_type)) + 
    ggtitle(label = paste0(year_wanted, " Early Run Estimation")) +
    #theme_Publication() +
    geom_vline(data = d, mapping = aes(xintercept = day_of_year, linetype = ltype), show.legend = FALSE) +
    geom_text(data=d, mapping=aes(x=day_of_year, y=0, label= event), size=5, angle=90, vjust=-0.4, hjust =-0.5) +
    scale_shape_manual(name = "Modeled by",
                        labels = c("Run Timing", "Genetics"), 
                        values = c(19, 17)) +
    theme(text = element_text(size =16), plot.title = element_text(vjust = -8, hjust = 0.05, size = 20), legend.justification = c(.5,0), legend.position = "bottom") +
    labs(y = "Cumulative run", x= "Day of the year") +
    coord_cartesian(xlim = c(150, 220))
  #early_cum
  ggsave(filename = paste0("figures/early_cum", year_wanted, ".png", sep = ""), device = png(), width = 7, height = 9, units = "in", dpi = 300)
  
  my_list <- list(df = df,"early_prop" = early_prop, "early_cum" = early_cum)
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


#This function looks at what genetics says the early run should look like. 
#This might be helpful if someone is trying to figure out starting values for the run timing method. 
early_look <- function(df, year_wanted){
  df <- data_prep_early(df, year_wanted)
  fit <- mix(as.mixdata(df), mixparam(mu=mean(df$day_of_year), sigma=sd(df$day_of_year)), dist="gamma", iterlim=5000)
  dist_plot(fit, year_wanted)
}

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

distribution_estimation_norms_SFX <- function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  
  mean_guess = c(mean(df$day_of_year) -30, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  sigma_guess = c(8.6, 8.6, 8.6)
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(consigma="SFX", fixsigma= c(TRUE, FALSE, FALSE)), dist=distibution_guess, iterlim=5000))  
}

distribution_estimation_norms_MFX <-  function(df, num_of_distributions = 3, mean_guess_given , sigma_guess_given, distibution_guess = 'gamma'){
  
  mean_guess = c(174, mean(df$day_of_year), mean(df$day_of_year)+30) #currently the default is for three distributions.
  sigma_guess = c(8.6, 8.6, 8.6)
  
  (fitpro <- mix(as.mixdata(df), mixparam(mu=mean_guess, sigma=sigma_guess),constr = mixconstr(conmu="MFX", fixmu= c(TRUE, FALSE, FALSE)), dist=distibution_guess, iterlim=5000))  
}

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
  df_year <- data_prep(df, year_wanted) 
  #graph_year(df_year)
  #fitpro <- distribution_estimation_norms(df_year) 
  #dist_plot(fitpro, year_wanted )
  #fitpro <- distribution_estimation_norms_MFX(df_year) 
  #dist_plot(fitpro, year_wanted ) 
  fitpro <- distribution_estimation_norms_SEQ(df_year) 
  dist_plot(fitpro, year_wanted )
}


dnormfit <- function(fit, x, dist_num = 1){
  #fit$parameters$pi[dist_num]*dnorm(x, fit$parameters$mu[dist_num],fit$parameters$sigma[dist_num])
  dnorm(x, fit$parameters$mu[dist_num],fit$parameters$sigma[dist_num])
}

percent_dist <- function(fit, x, dist_num = 1){
  dnormfit(fit, x, dist_num)/(dnormfit(fit, x, 1)+dnormfit(fit, x, 2)+dnormfit(fit, x, 3)) 
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
