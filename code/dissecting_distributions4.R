# notes ----
# dissecting distributions
# sarah.power@alaska.gov
# 7/7/2020

# load ----
source('code/distribution_functions4.R')

# data ----
data2019 <- read_csv('data/CMA_sockeye_2006_2019.csv') %>% 
  dplyr::rename(
    early_esc_genetics = early_esc, 
    late_esc_genetics = late_esc) %>%
  mutate(prop_early_genetics = 1 - prop_late, 
         esc = early_esc_genetics + late_esc_genetics,
         date = mdy(Date),
         year = year(date),
         day_of_year = yday(date)) # %>% #-> chig_data # convert the Date to its numeric equivalent 

datayear <- data2019

nas <- map_df(datayear, function(x) sum(is.na(x)))# figure out home many columns have NA's need purrr package

datayear[!complete.cases(datayear),]$date #print out lines of non complete data to check. This should be 0.

#datayear[is.na(datayear)] <- 0 # replace NAs with 0.

# Defining harvest Change df_data to the one you want to run through.
# df_data <- dfweironly <- datayear %>% mutate(harvest = 0)
#h_name = "h_none"
# df_data <- datayear %>% mutate(harvest = area_27110)
#h_name = "h_to_272_10_lagoon"
# df_data <- df18ocdb_272_20 <- datayear %>% mutate(harvest = area_27110 + area_27220)
#h_name = "h_to_272_20"
# df_data <- df18ocdb_272_30 <- datayear %>% mutate(harvest = area_27110 + area_27220 + area_27230)
#h_name = "h_to_272_30"
#df_data <- outercb <- datayear %>% mutate(harvest = area_27110 + area_27220 + area_27230 + area_27240)
#h_name = "h_to_272_40_OuterCB"
# df_data <- all_h <- datayear %>% mutate(harvest = datayear %>% select(area_27110:area_27560) %>% rowSums())
#h_name = "h_all"
# df_data <- kujulik <- datayear %>% mutate(harvest = datayear %>% select(area_27110:area_27253) %>% rowSums())
#h_name = "h_to_272_53_kujulik"
# df_data <- kumlik <- datayear %>% mutate(harvest = datayear %>% select(area_27110:area_27253, area_27262, area_27264) %>% rowSums())
#h_name = "h_to_272_64_kumlik"
# df_data <- west_kuj <- datayear %>% mutate(harvest = datayear %>% select(area_27110:area_27253, area_27390:area_27395) %>% rowSums())
#h_name = "h_west_to_kuj"
df_data <- miles_out45 <- datayear %>% mutate(harvest = datayear %>% select(
            area_27110, 
            area_27220, area_27230, area_27240, area_27251, area_27253, area_27262, area_27264,
            area_27390, area_27393, area_27395) %>% rowSums())
h_name = "miles_out45"

#So the stat areas are

#272-10	thru 272-40
#272-51, 272-53, 272-62, 272-64
#273-90 thru 273-95

#Checking remove later
sum(df_data$harvest)

# More data processing.
post_harvest <- function(df = datayear){
  df <- df %>% 
    mutate(run = esc + harvest, #catch is our estimate of what would have occured at the wier had there been no fishing, based on migration timing studies.
           run_early_gen = prop_early_genetics*run, 
           harvest_all = datayear %>% select(area_27110:area_27560) %>% rowSums(),
           #harvest_all == harvest currently attributed to Chignik - although not all of that harvest is used when we develop our runtiming estimates.
           run_all = esc + harvest_all,
           run_early_gen_all = prop_early_genetics*run) %>%
    select(-c(area_27110:area_27560))
  return(df)
}

df_data <- post_harvest(df_data)

#Checking remove later
sum(df_data$harvest_all)

year_vector <- c(2006:2008,2010:2019)

dgen <- read_csv('data/chig_genetics_by_weir_date.csv') %>%
  mutate(date = as.Date(paste(year, month, day, sep = '-')),
         day_of_year = yday(date))
  

# figures ----
out <- list()
#p <-par(mfrow = c(2,5))
for(i in 1:length(year_vector) ){
  png(file = paste0("figures/early_genetics", year_vector[i],h_name, ".png"), height = 4, width = 6, units = "in", res = 300)
  out[[i]] <- early_look(df_data, year_vector[i])
  dev.off()
}
#par(p)
eout <- do.call("rbind", out)

for(i in 1:length(year_vector) ){
  png(file = paste0("figures/mixture_auto", year_vector[i], h_name, ".png"), height = 4, width = 6, units = "in", res = 300)
  auto_year(df_data, year_vector[i])
  dev.off()
}

#plots <- list()
#for(i in 1:length(year_vector) ) {
#  plots[[i]] <- year_stats(df_data, year_vector[i])
#  ggsave(filename = paste0("figures/CDF",year_vector[i], h_name, ".png", sep = ""), device = png(), width = 6, height = 4, units = "in", dpi = 300)
#}

#do.call(grid.arrange, plots)
#plots[[2]]

#plots = lapply(1:9, function(.x) year_stats(,year_vector[i]))
#require(gridExtra)
#do.call(grid.arrange,  plots)

y06 <- year_stats(df_data, 2006) # 
y07 <- year_stats(df_data, 2007)
y08 <- year_stats(df_data, 2008)
y10 <- year_stats(df_data, 2010)
y11 <- year_stats(df_data, 2011)
y12 <- year_stats(df_data, 2012)
y13 <- year_stats(df_data, 2013)
y14 <- year_stats(df_data, 2014)
y15 <- year_stats(df_data, 2015)
y16 <- year_stats(df_data, 2016)
y17 <- year_stats(df_data, 2017)
y18 <- year_stats(df_data, 2018)
y19 <- year_stats(df_data, 2019)

dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()

fig <- cowplot::plot_grid(y06$early_cum, y07$early_cum, y08$early_cum, y10$early_cum, y11$early_cum, y12$early_cum, 
                          y13$early_cum, y14$early_cum, y15$early_cum, y16$early_cum, y17$early_cum, y18$early_cum, y19$early_cum, ncol = 3)

#multiplot(y06$early_cum, y07$early_cum, y08$early_cum, y10$early_cum, cols =2)

dev.off()
#Resources for grobs (graphical objects) 
#https://stackoverflow.com/questions/33114380/centered-x-axis-label-for-muliplot-using-cowplot-package
#https://stackoverflow.com/questions/17059099/saving-grid-arrange-plot-to-file
#http://www.sthda.com/english/wiki/wiki.php?id_contents=7930

y.grob <- textGrob(paste0("Count of early run sockeye using harvest ", h_name), gp=gpar(col="black", fontsize=15), rot=90)
x.grob <- textGrob("Day of year", gp=gpar(col="black", fontsize=15))
fig <- grid.arrange(arrangeGrob(fig, left = y.grob, bottom = x.grob))
ggsave(filename = paste0("figures/fig_year_stats_", h_name, ".png", sep = ""), device = png(), width = 9, height = 16, units = "in", dpi = 300, fig)


fig6 <- cowplot::plot_grid(y15$early_cum, y18$early_cum, ncol = 1)
#leg.grob <- legendGrob(c("Run Timing", "Genetics"), pch = c(19, 17), gp=gpar(col = 2:3, fill = "gray"))
y.grob <- textGrob(paste0("Cumulative count of early run sockeye"), gp=gpar(col="black", fontsize=15), rot=90)
x.grob <- textGrob(paste0("Day of year"), gp=gpar(col="black", fontsize=15))
fig6 <- grid.arrange(arrangeGrob(fig6, left = y.grob, bottom = x.grob))
ggsave(filename = paste0("figures/fig_6_15_18", h_name, ".png", sep = ""), device = png(), width = 4, height = 6, units = "in", dpi = 300, fig6)





#analysis ----

#grab the values from yYY$df that match those in dgen

m <- rbind(y06$df, y07$df, y08$df, y10$df, y11$df, y12$df, y13$df, y14$df, y15$df, y16$df, y17$df, y18$df, y19$df) %>%
  dplyr::select(date, year, day_of_year, dist_percent, prop_early_genetics, run_early_dis_all, run_early_gen_all, run_all, cum_run_all, 
                cum_run_dis_all, cum_run_gen_all) 
m <- m[,c("date", "year", "day_of_year", "dist_percent", "prop_early_genetics", "run_early_dis_all", "run_early_gen_all", "run_all", 
          "cum_run_dis_all", "cum_run_gen_all" , "cum_run_all")]

write.csv(m, file = paste0("figures/p_of_run_by_day_", h_name, ".csv", sep = ""))

p_of_run_by_year_ <- m %>%
  group_by(year) %>%
  summarise(dis_p = max(cum_run_dis_all, na.rm = TRUE)/max(cum_run_all, na.rm = TRUE),
            gen_p = max(cum_run_gen_all, na.rm = TRUE)/max(cum_run_all, na.rm = TRUE),
            diff = dis_p - gen_p)
write.csv(p_of_run_by_year_, file = paste0("figures/p_of_run_by_year_", h_name, ".csv", sep = ""))


m1 <- m %>%
  rename(genetics = cum_run_gen_all, runtiming = cum_run_dis_all) %>% 
  gather(key = "method", value = "cum_fish", runtiming, genetics)%>% 
  ggplot(aes(day_of_year, cum_fish, group = method)) +
  geom_line(aes(linetype = method), size = 1) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  geom_vline(data = d, mapping = aes(xintercept = day_of_year, linetype = ltype), show.legend = FALSE) +
  geom_text(data=d, mapping=aes(x=day_of_year, y=0, label= event), size=4, angle=90, vjust= -0.4, hjust = -0.5) +
  theme_bw() +
  facet_wrap(~year, ncol = 3) +
  theme(legend.position = "bottom") +
  labs(y = "Cumulative run", x= "Day of the year")
m1

ggsave(filename = paste0("figures/cum_early_run ", h_name, ".png", sep = ""), device = png(), width = 6, height = 9, units = "in", dpi = 300, m1)



new <- dplyr::inner_join(dgen, m, by = c("year", "day_of_year")) %>%
  dplyr::select(-month, -day, -chig_proportion, -date.y) %>%
  #black_proportion is the rawer genetics sample estimates, pro_early_genetics is the proportion after the black_proportion data has been fit to a logistic regresion. 
  rename(genetics_p = prop_early_genetics, genetics_sd = sd_bayes, runtiming_p = dist_percent) %>% 
  mutate(dif = genetics_p - runtiming_p) #%>%
  
  #select(date.x, year, day_of_year, sample_size, genetics_p, genetics_sd, dif,  runtiming_p)

ks.test(new$genetics_p, new$runtiming_p)

new1 <- new %>%
  #split(.$year) %>% 
  dplyr::filter(year == 2015) %>% 
  dplyr::select(genetics_p, runtiming_p)# %>% #, day_of_year)

ks.test(new1$genetics_p, new1$runtiming_p)
shapiro.test(new1$dif) # difference is normally distributed therefore we can do paired sample T-tests
t.test(new1$genetics_p, new1$runtiming_p, alternative = "two.sided")
#even 2015 and 2018 which are the years where the percentages differ the most the paired-t-test-sample pvlaue is over .20.


# statistical testing ----

data <- new
testing <- data %>%
  group_by(year) %>%
  nest() %>% 
  mutate(ks = purrr::map_dbl(data, ~ ks.test(.$genetics_p, .$runtiming_p)$p.value, data = .),
         shapiro = purrr::map_dbl(data, ~ shapiro.test(.$dif)$p.value, data = .),
         paired_t = purrr::map_dbl(data, ~ t.test(.$genetics_p, .$runtiming_p, alternative = "two.sided")$p.value, data = .)) %>% 
  select(-data)

write.csv(testing, file = paste0("figures/testing_", h_name, ".csv", sep = ""))


#proportion graph ----  
new <- new %>% 
  #gather(key = "method", value = "p", black_proportion, dist_percent)
  rename(genetics = genetics_p, runtiming = runtiming_p) %>% 
  gather(key = "method", value = "p", genetics, runtiming)

newplot <- new %>%
  ggplot(aes(day_of_year, p)) +
  geom_point(aes(shape = method, color = method), size = 2) +
  scale_shape_manual(values = c(1,2)) +
  theme_bw() +
  facet_wrap(~year, ncol = 3) +
  theme(legend.position = "bottom") +
  xlab(paste0("harvest = ", h_name))

ggsave(filename = paste0("figures/earlyrun_p_", h_name, ".png", sep = ""), device = png(), width = 6, height = 9, units = "in", dpi = 300)

h_name

###End of current code.


# Logistics Regression
glm.fit <- glm(p ~ day_of_year, data = new, family = binomial)

# might need to remove coding after this line

dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()

data <- df_data %>%
  group_by(year) %>%
  nest()

df06 <- data.frame(data[[2]][1])
df07 <- data.frame(data[[2]][2])
df08 <- data.frame(data[[2]][3])
df10 <- data.frame(data[[2]][4])
df11 <- data.frame(data[[2]][5])
df12 <- data.frame(data[[2]][6])
df13 <- data.frame(data[[2]][7])
df14 <- data.frame(data[[2]][8])
df15 <- data.frame(data[[2]][9])
df16 <- data.frame(data[[2]][10])
df17 <- data.frame(data[[2]][11])
df18 <- data.frame(data[[2]][12])

dfa <- df_data %>%
  select(year, day_of_year, run) %>%
  group_by(year) %>%
  nest() %>%
  mutate(fit = map(data, ~ distribution_estimation_norms_SEQ(df = .))) #,

fit06 <- dfa[[3]][1]
fit07 <- dfa[[3]][2]
fit08 <- dfa[[3]][3]
fit10 <- dfa[[3]][4]
fit11 <- dfa[[3]][5]
fit12 <- dfa[[3]][6]
fit13 <- dfa[[3]][7]
fit14 <- dfa[[3]][8]
fit15 <- dfa[[3]][9]
fit16 <- dfa[[3]][10]
fit17 <- dfa[[3]][11]
fit18 <- dfa[[3]][12]





         

