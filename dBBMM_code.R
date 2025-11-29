

# Load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(amt)


#load packages
pacman::p_load(tidyverse, janitor, lubridate, SDLfilter, adehabitatHR, 
               dplyr, sf, ggplot2, ggmap, readxl)

#from tutorial:
# library(ggplot2)
library(ggpubr)
library(ggspatial)
library(viridis)
library(scales)
library(glue)
library(grid)
# library(lubridate) 
# library(dplyr)
library(beepr)
library(reader)

library(move)
library(adehabitatHR)

# library(rgeos)
library(sf)
library(sp)
library(stringr)
library(terra)
library(cleangeo) 
library(plotly)

library(rnaturalearth)
library(tictoc)
library(MetBrewer)

library(tidyr)
library(purrr)
library(progressr) #progress bar in dbbmm function
library(readxl)
library(flextable)
library(jtools)



# Data prep --------------------------------------------------------------------

#dingo_data = main df

dingo_data <- read_csv("data/tanami_collars.csv") %>% filter(x > 0)

#create new combined datetime column
dingo_data <- dingo_data %>%
  mutate(
    datetime = dmy_hms(paste(cst_date, cst_time), tz = "Australia/Darwin")
  )


#check for duplicate timestamps

#one row per duplicate
dup <- dingo_data %>%
  group_by(dogname) %>%
  filter(duplicated(datetime)) %>%
  select(x, y, datetime, dogname, attempt_id, fix_id)
dup

#all duplicates
dup_all <- dingo_data %>%
  group_by(dogname, datetime) %>%
  filter(n() > 1) %>%
  arrange(dogname, datetime, attempt_id, fix_id)

dup_all


# take only the record with the most satellites
# Make sure it's sorted and de-duplicated per dog-time
dingo_clean <- dingo_data %>%
  arrange(dogname, datetime, hdop, desc(satellites)) %>%
  group_by(dogname, datetime) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(dogname, datetime)

# Drop tibble class
dingo_df <- as.data.frame(dingo_clean)

str(dingo_df$datetime)  # should be POSIXct

# 3. Minimal move() call (no data= yet)
dingo_movedata <- move(
  x      = dingo_df$x,
  y      = dingo_df$y,
  time   = dingo_df$datetime,       # already POSIXct, don’t wrap in as.POSIXct again
  proj   = CRS("+proj=longlat +datum=WGS84"),
  animal = dingo_df$dogname         # character is fine
)

dingo_movedata@data <- dingo_df

str(dingo_movedata)

# Was having issues with order ob objects not the same as trackID; fixed after doing the steps above through chatgpt
# any(duplicated(dingo_clean$datetime))                         # across whole dataset
# any(duplicated(dingo_clean[, c("dogname", "datetime")]))      # per dog



# model parameters --------------------------------------------------------

# According to Kranstauber et al 2012: over a range of possible margin window
# size combinations, margins of 9–11 locations and window sizes of around 30
# seemed to perform best (Default is margin = 11, window = 31).
#
# The window size should be >2x margin, and the margin be at least 3 pts.
#
# Increasing window size increases the reliability of the motion variance
# estimate, but at the cost of missing rapid changes in behaviour. Increasing
# the margin increases the likelihood of identifying weak breakpoints, at the
# cost of not detecting breakpoints due to limited observations in window.
#
# The window and margin should be chosen based on biologically relevant time
# intervals that you expect behavioural changes to occur. The values will also
# depend on the interval of tracking data.



# Prep constant parameters and output locations

## constants that stay the same across snakes
raster_res  <- 30
margin      <- 3
window_size <- 12   

## output folders
# one base folder for everything from this run
output_base <- file.path("dBBMMs_pipeline_output")

out_dir_dbb <- file.path(output_base, "dBBMMs_saved")
out_dir_var <- file.path(output_base, "Motion_variance_saved")
out_dir_log <- file.path(output_base, "logs")


# Shapefiles --------------------------------------------------------------

# Get the target projection from the spatial data provided
#load shapefiles

# load distance to artificial food layer
dist_artificial_food <- terra::rast("spatial_layers/artificial_food/dist_artificial_food.tif")

# load distance to artificial water layer
dist_artificial_water <- terra::rast("spatial_layers/artificial_water/dist_artificial_water.tif")

# load distance to roads layer
dist_roads <- terra::rast("spatial_layers/roads/dist_roads.tif")

#check projection
st_crs(dist_roads)

#   st_transform(crs = st_crs(dist_roads))


#52K is the UTM zone
#code 32753

## target projection: UTM Zone 55S
utm52 <- sf::st_crs(32753)         
utm52_sp <- as(utm52, "CRS")       # convert to sp::CRS for move/spTransform



# dingo_movedata must already exist and have valid WGS84 lon/lat coordinates

dingo_movedata <- spTransform(dingo_movedata, utm52_sp)


# Test one dBBMM models --------------------------------------------------------

# Pull out two individuals to compare

# shaun = intermediate example
# boss = mine example

#extract shaun only
movedata_shaun <- dingo_movedata[[c("shaun")]]

# dBBMM - don't run unless necessary (can take a while)
dBBMM_shaun <- brownian.bridge.dyn(movedata_shaun, 
                                     ext = 0.75,
                                     raster = 30, #raster pixel size, should always be bigger than error
                                     # dimSize = 10, #only used if raster is not set
                                     margin = 3, #previously 11
                                     window.size = 13, #previously 31
                                     location.error = 10,
                                     timestep = median(timeLag(movedata_shaun,"hours")/15)
                                   )

#the output is a raster, where each cell has a probability of the animal occurring there

#plot
plot(sqrt(dBBMM_shaun), main = "Shaun UD") #quick hack that changes color scale to highlight the colors better
# points(dBBMM_shaun, col="red",  cex=.5, pch=20)
contour(dBBMM_shaun, add = TRUE,
        levels = c(0.5, 0.75, 0.95)) #add contour lines to plot


# Motion variance

#extract the variance of the movement
dataVar_shaun <- brownian.motion.variance.dyn(movedata_shaun,
                                               margin = 3,
                                               window.size = 13,
                                               location.error = 10) #5 m

getMotionVariance(dataVar_shaun)

#plot with base R
plot(timestamps(dataVar_shaun),
     getMotionVariance(dataVar_shaun),
     type = "s",
     xlab = "Time",
     ylab = "Variance",
     main = "Variance of movement"
)


##try plotting with ggplot instead

# Extract timestamps and variance
df_shaun <- data.frame(
  time = as.POSIXct(timestamps(dataVar_shaun)), # Convert timestamps to POSIXct for proper time handling
  variance = getMotionVariance(dataVar_shaun) # Extract motion variance
)



# Mine example ------------------------------------------------------------


#extract boss only
movedata_boss <- dingo_movedata[[c("boss")]]

# dBBMM - don't run unless necessary (can take a while)
dBBMM_boss <- brownian.bridge.dyn(movedata_boss, 
                                   ext = 0.75,
                                   raster = 30, #raster pixel size, should always be bigger than error
                                   # dimSize = 10, #only used if raster is not set
                                   margin = 3, #previously 11
                                   window.size = 13, #previously 31
                                   location.error = 10,
                                   timestep = median(timeLag(movedata_boss,"hours")/15)
)

#the output is a raster, where each cell has a probability of the animal occurring there

#plot
plot(sqrt(dBBMM_boss), main = "boss UD") #quick hack that changes color scale to highlight the colors better
# points(dBBMM_boss, col="red",  cex=.5, pch=20)
contour(dBBMM_boss, add = TRUE,
        levels = c(0.5, 0.75, 0.95)) #add contour lines to plot


# Motion variance

#extract the variance of the movement
dataVar_boss <- brownian.motion.variance.dyn(movedata_boss,
                                              margin = 3,
                                              window.size = 13,
                                              location.error = 10) #5 m

getMotionVariance(dataVar_boss)

#plot with base R
plot(timestamps(dataVar_boss),
     getMotionVariance(dataVar_boss),
     type = "s",
     xlab = "Time",
     ylab = "Variance",
     main = "Variance of movement"
)


##try plotting with ggplot instead

# Extract timestamps and variance
df_boss <- data.frame(
  time = as.POSIXct(timestamps(dataVar_boss)), # Convert timestamps to POSIXct for proper time handling
  variance = getMotionVariance(dataVar_boss) # Extract motion variance
)



# Plot motion variance ----------------------------------------------------

# combine into one df
# Add individual labels
df_boss  <- df_boss  %>% mutate(individual = "boss")
df_shaun <- df_shaun %>% mutate(individual = "shaun")

# Combine the two data frames
df_combined <- bind_rows(df_boss, df_shaun)

# plot both individuals on the same plot
mv_combined <- ggplot(df_combined, aes(x = datetime, 
                        y = motion_variance, 
                        color = dogname)) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Dynamic Brownian Motion Variance",
    x = "Time",
    y = "Motion variance",
    color = "Individual"
  )

mv_combined

# seems like boss might have an error in his data

# Plot mv AND distance ---------------------------

#rename for joining
df_combined <- df_combined %>%
  rename(dogname = individual,
         datetime = time,
         motion_variance = variance)

dingo_data_sel <- dingo_data %>%
  filter(dogname %in% c("shaun", "boss")) %>%
  select(datetime, dogname, mine_away, mintemp, maxtemp, landunit, disttip, distwater, distprimroad)

# add distance to mine to Shaun's data
df_combined <- df_combined  %>%
  left_join(., dingo_data_sel, by = c("dogname", "datetime"))
  


# plot shaun --------------------------------------------------------------


# create another plot with distance to mine over time, to compare the motion variance to mine distance

#plot Shaun's mv
shaun_mv <- ggplot(df_combined %>% filter(dogname == "shaun"), 
                   aes(x = datetime, 
                        y = motion_variance)) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Dynamic Brownian Motion Variance",
    x = "Time",
    y = "Motion variance")

shaun_mv

# Plot shaun's distance to tip
shaun_disttip <- ggplot(df_combined %>% filter(dogname == "shaun"), 
                  aes(x = datetime, 
                      y = disttip)) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distance to tip over time",
    x = "Time",
    y = "Distance to tip")

shaun_disttip


library(cowplot)

shaun_grid <- plot_grid(
  shaun_mv,
  shaun_disttip,
  ncol = 1,   # 1 column
  nrow = 2    # 2 rows
)


# plot boss ---------------------------------------------------------------



#plot boss's mv
boss_mv <- ggplot(df_combined %>% filter(dogname == "boss"), 
                   aes(x = datetime, 
                       y = motion_variance)) +
  geom_line(linewidth = 1) +
  # scale_y_continuous(
  #   breaks = seq(0, 400, 50),
  #   labels = seq(0, 400, 50),
  #   limits = c(0, 400),
  # ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Dynamic Brownian Motion Variance",
    x = "Time",
    y = "Motion variance")

boss_mv

# Plot boss's distance to tip
boss_disttip <- ggplot(df_combined %>% filter(dogname == "boss"), 
                        aes(x = datetime, 
                            y = disttip)) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distance to tip over time",
    x = "Time",
    y = "Distance to tip")

boss_disttip


boss_grid <- plot_grid(
  boss_mv,
  boss_disttip,
  ncol = 1,   # 1 column
  nrow = 2    # 2 rows
)

boss_grid





# Plots -------------------------------------------------------------------

mv_combined

boss_grid

shaun_grid

