load("~/Downloads/6month.rdata")
install.packages("usmap")
install.packages("sf")
library(sf)
library(ggplot2)
library(dplyr)
library(lubridate)
usa <- usmap::us_map()

locs <- train.df %>% select(lat, long) %>% unique
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
locs_sf <-  st_as_sf(x = locs, coords = c("long", "lat"), crs = projcrs)
ggplot() + geom_sf(data = usa)    +
  geom_sf(data = locs_sf) + xlab("") +
  geom_sf(data = locs_sf[1:2,], col = 2)

(locs[1,])

data_for_first_loc <- inner_join(train.df, locs[1,])

data_for_first_loc <- data_for_first_loc %>% 
  mutate(Month = month(DateTime))

ggplot(data = data_for_first_loc, aes(x = DateTime, y = GHI_Meas)) + geom_line()

table(data_for_first_loc$Month)

ggplot(data = data_for_first_loc %>% filter(Month == 9 ), aes(x = DateTime, y = GHI_Meas)) + geom_line()


