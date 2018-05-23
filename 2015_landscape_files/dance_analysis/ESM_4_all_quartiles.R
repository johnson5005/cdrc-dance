######## This script has been adapted for a new dataset by Douglas Sponsler. Unused portions of the original script have been commented out but not deleted. Additions are followed by the tag # Sponsler:

# Electronic Supplementary Material 4 - R script to simulate dances from known waggle
# dance durations and headings.
# ------------------------------------------------------------------------------------

# Article title: Incorporating variability in honey bee waggle dance decoding improves
# the mapping of communicated resource locations

# Journal: Journal of Comparative Physiology A

# Authors: Roger Sch?rch, Margaret J. Couvillon, Dominic D. R. Burns, Kiah
# Tasman, David Waxman and Francis L. W. Ratnieks

# Corresponding author: Roger Sch?rch, Evolution, Behaviour and
# Environment, School of Life Sciences, University of Sussex, Brighton,
# BN1 9QG, United Kingdom, R.Schuerch@sussex.ac.uk

# Last revised: 2013-08-12
# Older revision: 2013-05-27

# script files that go with this script:
# ESM_3.jag

# data files that go with this script:
# ESM_5.csv

# the file will output a comma separated value file and an ASC raster file
# you must create a "data" folder within the folder from which you are running the scripts
# or adapt the paths for the instructions below

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg
library('googlesheets')  # Read from GoogleDocs
library('magrittr')        # Pipes
library('oce')  # Calculation of azimuth

# Set working directory and create "data" subdirectory
setwd("~/Downloads/version-control/cdrc-dance/2015_landscape_files/dance_analysis") # path to the working directory
dir.create("data", showWarnings = FALSE)

## Download data from GoogleSheets using 'googlesheets' library.  Store as csv file specified in "wagglefile"
waggleFile <- "2015_CDRC_dance.csv" # Set the name for the file to use
# "2015 CDRC Dance Analysis" https://docs.google.com/spreadsheets/d/15EX-j0CAlhNjhP0baeBAw4FuZGjYdRkKQ1PZoj7Q_Qg/edit?usp=sharing
gs_key("15EX-j0CAlhNjhP0baeBAw4FuZGjYdRkKQ1PZoj7Q_Qg") %>%
  gs_download(ws = "Data", to = waggleFile, overwrite = TRUE)

## Load and process dance data
waggleData <- read.csv(waggleFile) # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1) # Sponsler: a flag field removes empty or incomplete lines

# Divide into groups
 #q <- quantile(waggleData$day)
 #waggleData$quantile <- 0
 #waggleData[waggleData$day <= q[5],]$quantile <- 4
 #waggleData[waggleData$day <= q[4],]$quantile <- 3
 #waggleData[waggleData$day <= q[3],]$quantile <- 2
 #waggleData[waggleData$day <= q[2],]$quantile <- 1

# Divide into 3 groups
# q <- quantile(waggleData$day, probs = seq(0,1,(1/3)))
# waggleData$quantile <- 0
# waggleData[waggleData$day <= q[4],]$quantile <- 3
# waggleData[waggleData$day <= q[3],]$quantile <- 2
# waggleData[waggleData$day <= q[2],]$quantile <- 1

# Divide into 2 groups
#q <- quantile(waggleData$day, probs = seq(0,1,0.5))
#waggleData$quantile <- 0
#waggleData[waggleData$day <= q[3],]$quantile <- 2
#waggleData[waggleData$day <= q[2],]$quantile <- 1

# Divide based on days recommended by Rodney: May 8 and May 19
waggleData$quantile <- 3
waggleData[waggleData$day <= 19,]$quantile <- 2
waggleData[waggleData$day <= 8,]$quantile <- 1

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle", 
								 unit = "radian", 
								 rotation = "clock", 
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dancer.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)
   
# store the subset of dances going to that feeder in waggleData
# Establish locations 
# the UTM 17N (EPSG:26917) northing/easting of the hives in meters
location <- data.frame(
  site= c("FSR", "HONEYRUN", "MOOREMAN", "MECHANICSBURG"),
  easting = c(292263.656365, 317334.853833, 295906.044842, 278907.984269 ),
  northing = c(4426271.356893, 4391383.561963, 4406432.40089, 4430377.392885) )

# make the data properly circular
waggleData$heading <- circular(waggleData$heading.radians,
                               type = "angle", 
							   unit = "radian", 
							   rotation = "clock", 
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

## Set up list of analyses to run
waggleData$runGroup <- paste(waggleData$hive, waggleData$quantile, sep=".")
toRun <- unique(waggleData$runGroup)

## Iterate through list 
for (i in 1:length(toRun)){
  print(toRun[i])
# Subset to get just the data for each site and quantile
wD <- waggleData[waggleData$runGroup == toRun[i],]

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- location[location$site == site,]$easting  # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- location[location$site == site,]$northing # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(wD$dancer.id)){
  cat(paste(i, "of", length(wD$dancer.id), "\n"))
  # choose only the i^th dance
  tempData <- wD[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$mean_duration.sec)
  x2 <- rep(NA, length(tempData$mean_duration.sec))
  y2 <- tempData$mean_duration.sec

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]
  
  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dancer.id = rep(tempData$dancer.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dancer.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }
  
  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler: 
  
  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dancer.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)
  
  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dancer.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 5000,
                                 hiveEasting + 5000),
                               northing = c(hiveNorthing, hiveNorthing - 5000,
                                 hiveNorthing + 5000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
# prepare a file name for the raster
rasterName <- paste(i, "May2015.tif", sep="_")
writeRaster(new.data.rast, filename = rasterName, format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer

}



### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)
