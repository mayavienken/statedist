## Prepare case study data

library(dplyr)
library(LaMa)
library(leaflet)
library(zoo)
library(moveHMM)
library(scales)

## load data + overview

# the dataset can be found here (size of 416.2 MB): 

# Bastille-Rousseau G, Yackulic CB, Gibbs J, Frair JL, Cabrera F, Blake S. 2019. 
# Data from: Migration triggers in a large herbivore: Galápagos giant tortoises navigating resource gradients on volcanoes.
# Movebank Data Repository. https://doi.org/10.5441/001/1.6gr485fk

# data <- read.csv('case_study_tortoises/data/GalapagosTortoise.csv', header = TRUE) 
head(data)
dim(data) # 1841253 x 24

table(data$individual.local.identifier)
length(table(data$individual.local.identifier)) # 96 individual animals

## prep data

# check for GPS-Status (eobs.status), remove everything where status == 'D'
# definition: 
# e-obs status: The record status, from e-obs GPS/accelerometer tags; Allowed values are
# A = position and time within accuracy masks
# B = only time of week and weeknumber valid
# C = only weeknumber valid
# D = no valid data
table(data$eobs.status)
data.1 <- data[data$eobs.status == "A", ]
table(data.1$eobs.status)
data$timestamp

# remove manually marked outliers
table(data.1$manually.marked.outlier)
data.1 <- data.1[data.1$manually.marked.outlier != "true", ]

# only keep selected columns
data.1 <- data.1[, c("timestamp", "location.long", "location.lat", "individual.local.identifier", "eobs.temperature", "heading", "height.above.ellipsoid", "ground.speed")]

dim(data)
dim(data.1)

sum(is.na(data[data$individual.local.identifier=="Carolina",]))
sum(is.na(data.1[data.1$individual.local.identifier=="Carolina",]))

data.1$timestamp <- substr(data.1$timestamp,1,13)

## selecting one specific tortoise "Carolina"
caro <- data.1[data.1$individual.local.identifier == "Carolina", ]
caro <- caro[, c("timestamp", "location.long", "location.lat", "eobs.temperature", "heading", "height.above.ellipsoid", "ground.speed")]

dim(caro)
plot(caro$eobs.temperature, type="l")
sum(is.na(caro$timestamp))
caro <- aggregate(caro[, 2:7],
                    by = list(timestamp = caro$timestamp),
                    FUN = function(x) mean(x, trim = 0.2, na.rm = TRUE))

caro$timestamp <- as.POSIXct(caro$timestamp, format="%Y-%m-%d %H")
sum(duplicated(caro$timestamp))
caro <- caro[!duplicated(caro$timestamp), ]
#remove na
sum(is.na(caro$timestamp))
caro <- caro[!is.na(caro$timestamp), ]
caro.zoo <- zoo(caro[,-1], caro[,1])
caro <- merge(caro.zoo,zoo(,seq(start(caro.zoo),end(caro.zoo),by="hour")),all=TRUE)
rm(caro.zoo)

plot(caro$eobs.temperature, type="l")

caro = data.frame(index(caro), as.data.frame(caro))
names(caro)[1] = "timestamp"
rownames(caro) <- NULL

coords = caro[,2:3]
df.caro <- prepData(coords, type="LL", coordNames = c("location.long", "location.lat"))

df.caro$step <- df.caro$step*1000
max(df.caro$step, na.rm = T)
quantile(df.caro$step, 0.996, na.rm = T)
df.caro$step[df.caro$step > 150] <- NA
idx=which(df.caro$step<=0)
df.caro$temperature <- caro$eobs.temperature
df.caro$temperature <- na.approx(df.caro$temperature)
df.caro$timestamp = caro$timestamp

sum(is.na(df.caro$temperature))

location <- leaflet(df.caro)
location <- addTiles(location)
location <- addPolylines(location, lng = ~x, lat = ~y)
location


# Get blocks of NAs of covariate
ts_index <- index(caro)
temp_vals <- caro$eobs.temperature
na_idx <- is.na(temp_vals)
rle_na <- rle(na_idx)
ends <- cumsum(rle_na$lengths)
starts <- c(1, head(ends, -1) + 1)
na_blocks <- data.frame(
  start = ts_index[starts[rle_na$values]],
  end   = ts_index[ends[rle_na$values]],
  length_hours = rle_na$lengths[rle_na$values]
)
View(na_blocks)
sort(na_blocks$length_hours, decreasing = TRUE)
summary(na_blocks$length_hours)

# 2009-05-14 17:00:00 - 2009-08-11 21:00:00 # first track
# delete in between
# 2009-09-01 18:00:00 - 2009-10-14 01:00:00 # second track 
# delete in between
# 2009-12-17 16:00:00 - 2011-10-29 23:00:00 # third track

block_1_start <- as.POSIXct("2009-05-14 17:00:00")
block_1_end   <- as.POSIXct("2009-08-11 21:00:00")

block_2_start <- as.POSIXct("2009-09-01 18:00:00")
block_2_end   <- as.POSIXct("2009-10-14 01:00:00")

block_3_start <- as.POSIXct("2009-12-17 16:00:00")
block_3_end   <- as.POSIXct("2011-10-29 23:00:00")

df.caro <- df.caro[!(df.caro$timestamp >= block_1_end & df.caro$timestamp <= block_2_start), ]
df.caro <- df.caro[!(df.caro$timestamp >= block_2_end & df.caro$timestamp <= block_3_start), ]

df.caro$track_id <- NA
df.caro$track_id[df.caro$timestamp <= block_1_end] <- 1
df.caro$track_id[df.caro$timestamp >= block_2_start & df.caro$timestamp <= block_2_end] <- 2
df.caro$track_id[df.caro$timestamp >= block_3_start] <- 3

table(df.caro$track_id, useNA = "ifany")

dim(df.caro) # 19504 x 8

#save(df.caro, file = "df_caro2.RData")