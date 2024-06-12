# analysis of movement data from GPS collared muskoxen from Greenland (https://zenodo.org/records/3768080)
rm(list = ls())

library(ctmm)
library(sf)
library(moveHMM)
# import data
DAT <- read.delim("data/muskox_tracking_data_summer_20200417.txt", colClasses = c(datetime = "POSIXct"))
# select columns to create telemetry object, rename them
telemDAT <- DAT[,c(1,6,5,4)]
names(telemDAT) <- c("individual.local.identifier", "timestamp", "location.lat", "location.long")
# create telemetries as list
telems <- as.telemetry(telemDAT, datum = "+proj=utm +zone=27 +datum=WGS84 +units=m +no_defs")
# saveRDS(telems, "outputs/greenland_telemetries.rds")

# Plot positions
plotcols <- hcl.colors(23,"Dynamic")
plot(telems, col = plotcols, error = F)

# # fit ctmm. This takes a while, uncomment below to run
# GUESS <- lapply(telems, ctmm.guess, interactive=F)
# FIT <- mapply(ctmm.select, telems, GUESS, verbose = T, trace = 1, cores = 4)
# saveRDS(FIT, "outputs/greenland_ctmm.rds")
# See summary
lapply(FIT, "[[", 1)|>lapply(summary)

# Get utilization distribution
# import polygon boundaries
bounds <- st_read("data/gadm41_GRL_0.json") |> as_Spatial() |> 
  sp::spTransform(CRSobj = "+proj=utm +zone=27 +datum=WGS84 +units=m +no_defs")
raster::plot(bounds)
FIT1 <- lapply(FIT, "[[", 1)
# # utilization distribution through AKDE. un comment to run
# UD <- akde(telems, FIT1, SP = bounds)
saveRDS(UD, "outputs/greenland_UDs.rds")

# Hidden Markov Model to assess behavior. The goal is to assign behaviors
# (resting, foraging, moving), and estimate the distance travelled on average
# when foraging. We assume that transmission occurs when ingesting contaminated
# soil, so is likely during foraging bouts.
# 