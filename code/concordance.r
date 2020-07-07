# notes ----
# Determining possible harvest allocation schemes which should help with 
# run timing estimates of early run transitional timing.
# sarah.power@alaska.gov
# 10/29/2018

# load ----
library(here)

# data ----
# load ----
source('code/distribution_functions.R')
#library(lubridate)
getwd()

# data ----
chig_data <-read_csv('data/ChigISGrunappt2006-2017catch.by.district.csv') %>% 
  dplyr::select( Date, Propotionearly, Earlyesc, Lateesc) %>%
  dplyr::rename(
    date = Date,
    prop_early_genetics = Propotionearly, 
    early_esc_genetics = Earlyesc, 
    late_esc_genetics = Lateesc
  ) %>%
  mutate(esc = early_esc_genetics + late_esc_genetics,
         date_weir = lubridate::mdy(date),
         date = date_weir,
         day_of_year = lubridate::yday(date_weir))%>% # convert the Date to its numeric equivalent
  filter(year(date_weir) %in% c(2006:2016))

harvest_data <- read_csv('data/harvest19692017.csv') %>% 
  filter(year %in% c(2006:2016)) %>%
  mutate(date_harvest = mdy(date_begin),
         day_of_year_harvest = yday(date_harvest),
         date_weir = weir_date(area, date_harvest), 
         allocation_percent = allocation(area, date_harvest),
         chig_catch = round(fish_count*allocation_percent, 0),
         month_of_weir = month(date_weir)
  )


# This section to figure out patterns of harvest an what might be different about 2008 and 2015 in particular
#with(harvest_data, table(month_of_weir,year, area ))
unique(harvest_data$area)

xyz <- harvest_data %>% 
  #filter(area %in% c("SEDM")) %>% #, "Western")
  group_by(year, month_of_weir) %>%
  summarise(month_harv = sum(chig_catch)) %>%
  mutate(month_of_weir=paste('m', month_of_weir, sep="_")) 

year_harv <- harvest_data %>%
  group_by(year) %>%
  summarise(year_harv = sum(chig_catch))

xyz <- merge(xyz, year_harv, by = "year") 

#class(xyz)
xyz <- xyz %>%
  mutate(permonth_harv = month_harv/year_harv) %>%
  select(year, month_of_weir, permonth_harv) %>%
  spread(month_of_weir, permonth_harv)


xyz2 <-harvest_data %>% 
  group_by(year, month_of_weir) %>%
  summarise(month_harv = sum(fish_count))
###end section

#create a wide data set
harv <- harvest_data %>% 
  dplyr::select(date_weir, year, month_of_weir, area, chig_catch) %>% 
  tidyr::spread(area, chig_catch)

harv2 <- harv %>% 
  group_by(year, month_of_weir) %>%
  summarise( Central= sum(Central, na.rm = TRUE),
             Chignik= sum(Chignik, na.rm = TRUE),
             Dolgoi = sum(Dolgoi, na.rm = TRUE ),
             Eastern= sum(Eastern, na.rm = TRUE),
             Igvak  = sum(Igvak, na.rm = TRUE),
             Ikatan = sum(Ikatan, na.rm = TRUE),
             Perryville= sum(Perryville, na.rm = TRUE),
             SEDM   = sum(SEDM, na.rm = TRUE),
             Shumagins = sum(Shumagins, na.rm = TRUE),
             Unimak = sum(Unimak, na.rm = TRUE),
             Western= sum(Western, na.rm = TRUE)
             )

harv3 <- harv %>% 
  group_by(year) %>%
  summarise( Central= sum(Central, na.rm = TRUE),
             Chignik= sum(Chignik, na.rm = TRUE),
             Dolgoi = sum(Dolgoi, na.rm = TRUE ),
             Eastern= sum(Eastern, na.rm = TRUE),
             Igvak  = sum(Igvak, na.rm = TRUE),
             Ikatan = sum(Ikatan, na.rm = TRUE),
             Perryville= sum(Perryville, na.rm = TRUE),
             SEDM   = sum(SEDM, na.rm = TRUE),
             Shumagins = sum(Shumagins, na.rm = TRUE),
             Unimak = sum(Unimak, na.rm = TRUE),
             Western= sum(Western, na.rm = TRUE)
  )

harv4 <- harv3 %>% 
  group_by(year) %>%
  mutate(yearly_harv = sum(Central,Chignik, Dolgoi, Eastern, Igvak, Ikatan, Perryville, SEDM, Shumagins, Unimak, Western, na.rm = TRUE),
         Central= round(Central/yearly_harv*100),
         Chignik= round(Chignik/yearly_harv*100),
         Dolgoi = round(Dolgoi/yearly_harv*100),
             Eastern= round(Eastern/yearly_harv*100),
             Igvak  = round(Igvak/yearly_harv*100),
             Ikatan = round(Ikatan/yearly_harv*100),
             Perryville= round(Perryville/yearly_harv*100),
             SEDM   = round(SEDM/yearly_harv*100),
             Shumagins = round(Shumagins/yearly_harv*100),
             Unimak = round(Unimak/yearly_harv*100),
             Western= round(Western/yearly_harv*100))
# We can see in 2015 there was quite a bit of Chignik attributed harvest caught in the SEDM & Western region, 
# with relatively less caught in Chignik and Central.  
  

catch <-aggregate( chig_catch ~ date_weir, data=harvest_data, FUN=sum)

harv <- harv %>% 
  full_join(catch, by = "date_weir") #catch is our estimate of what would have occured at the wier had there been no fishing, based on migration timing studies.

harv$chig_catch%>% 
  replace_na(0)

#Some harvest would have arrived at the weir after the weir is removed, or when the weir was down.
something <-anti_join(harv, chig_data, by ="date_weir") #most of this is 2009 data since we don't have 2009 weir data

all <- full_join(chig_data, harv, by = "date_weir") %>%
  mutate(day_of_year = yday(date_weir))

#need to replace na with 0 so we can sum the columns
all$chig_catch <- all$chig_catch %>% 
  replace_na(0)

all$esc <- all$esc %>% 
  replace_na(0)

all <- all%>%
  mutate(run = chig_catch + esc,
         run_early_gen = prop_early_genetics*run)

p6 <- year_stats(all, 2006)
p7 <- year_stats(all, 2007)
p8 <- year_stats(all, 2008)
p10 <- year_stats(all, 2010)
p11 <- year_stats(all, 2011)
p12 <- year_stats(all, 2012)
p13 <- year_stats(all, 2013)
p14 <- year_stats(all, 2014)
p15 <- year_stats(all, 2015)
p16 <- year_stats(all, 2016)

grid.arrange(p6, p7, p8, p10, p11, p12, nrow = 3)
grid.arrange(p11, p12, p13, p14, p15, p16, nrow = 3)

year_vector <- c(2006:2008,2010:2017)

par(mfrow = c(3,2))
auto_year(all, 2006)
auto_year(all, 2007)
auto_year(all, 2008)
auto_year(all, 2010)
auto_year(all, 2011)
auto_year(all, 2012)
auto_year(all, 2014)
auto_year(all, 2015)
auto_year(all, 2016)
auto_year(all, 2017)

auto_year(all, 2008)
auto_year(all, 2012)
auto_year(all, 2014)
auto_year(all, 2015)

all %>% select()

bar <- ggplot(all, aes(x = date, y =)) +
  geom_bar(aes(fill = run), position = position_stack())

