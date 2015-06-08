## list_files lists all the files in one of the data folders and counts the number of files. Default time1 = "present"
list_files <- function(time){
  path <- paste("data/", time, "/", sep="")
  list <- list.files(path = path, pattern="*.asc")
  return(list)
}

count_list <- function(list){
  n <- base::length(list)
  n
}

## list_to_names takes the list of files and returns a list of species "namelist"
list_to_names <- function(list){
  parts <- strsplit(list, "_")
  genus <- unlist(lapply(parts,"[",1))
  species <- unlist(lapply(parts,"[",2))
  namelist <- paste(genus, species, sep="_")
  return(namelist)
}

## set_LATLON sets the CRS system. default is WGS84
set_LATLON <- function(CRS = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"){
  LATLON<-CRS(CRS)
  return (LATLON)
}

##raster_to_points takes an .asc file, reads it as a raster, transforms it to points, produces a table of points coordinates, and count the number of points
time_raster <- function(time){
  file_name <- paste("data/", time, "/", spname, "_", time, ".asc", sep="")
  timeraster <-raster(file_name)
  return(timeraster)
}
#convert raster to points where value is 1
raster_to_points <- function(timeraster){
  timeraster2 <- reclassify(timeraster, c(-Inf,0.66,0,0.66,Inf,1))
  timePoints <- rasterToPoints(timeraster, fun=function(timeraster){timeraster==1})
  return(timePoints)
}
#convert raster to points where value is 1 or 0 (used for coastal distance)
raster_to_points_coast <- function(timeraster){
  timeraster2 <- reclassify(timeraster, c(-Inf,0.66,0,0.66,Inf,1))
  timePoints <- rasterToPoints(timeraster, fun=function(timeraster){timeraster==1 | timeraster==0})
  return(timePoints)
}


number_points <- function(timePoints){
  z <- nrow(timePoints)
  return(z)
}







points_table <- function(timePoints){
  time_points_table<-cbind(timePoints[,1],timePoints[,2])
  return(time_points_table)
  }

#read shape file for individual species graph
read_shapefile <- function(shape){
  x <- paste("data/", shape, ".shp", sep="")
  shapefile <- readShapeLines(x,  proj4string=LATLON)
  return(shapefile)
  }

## make table of distance of 20 furthest points from centroid
get_20furthest_points_in_wedge <- function(n, t, time, function_raspts = "raster_to_points", ..., verbose=TRUE){
    if (verbose) {
        message("Going through points in individual wedge...")
                  }
  time_ras <- time_raster(time)
  time_points <- function_raspts(time_ras)
  
 for (p in 1:n){
    x <- number[p,1]
    c <- which(t[,1]==x)
    coord <-  time_points[c,1:2]
    matrixWedge <- rbind(matrixWedge, coord)
    
    dist<- distMeeus(combined_centroid, coord, a=6378137, f=1/298.257223563)
    distWedge <- rbind(distWedge, dist)
  } #close loop to gather all the points in distances in that wedge
  return(distWedge)
}

get_20furthest_coordinates_in_wedge <- function(n, t, time, function_raspts = "raster_to_points",..., verbose=TRUE){
  time_ras <- time_raster(time)
  time_points <- function_raspts(time_ras)
  
  for (p in 1:n){
    x <- number[p,1]
    c <- which(t[,1]==x)
    coord <-  time_points[c,1:2]
    matrixWedge <- rbind(matrixWedge, coord)

  } #close loop to gather all the points in distances in that wedge
  return(matrixWedge)
}


mean_dist_20furthest_points <- function(distWedge,..., verbose=TRUE){
  #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
  if (verbose) {
    message("Finding 20 furthest points...")
  }
  
  order <- distWedge[order(distWedge, decreasing=TRUE)]
  furthest <- order[1:20]
  return(furthest)
}




## take the 20 furthest points and calculate the point at wedge for the mean distance from the 20 furthest
mean_points_table_for_wedge <- function(furthest, distWedge, matrixWedge){
  furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
  for (x in 1:20){
    y <- which(distWedge == furthest[x])
    y <- y[1]
    z <- matrixWedge[y,1:2]
    furthestMatrix <- rbind(furthestMatrix,z)
  }
  #meanmin is the centroid for the min-max wedge
  points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
  nP <- nrow(points)
  #Here we add a condition in case n=1. In this case, geomean chokes and stops the script. If there is only 1 point in that wedge, it is its own mean. 
  if(nP>1){
    mean<-geomean(points)
  }
  if(nP<2) {
    mean <- points
  }
  return(mean)
}

##summarize table data by wedge (rotation column)
summary_by_wedge <- function(datum, var_sum){
summary <- summaryBy(datum[,var_sum]~rotation, datum, FUN=function(x){c(m=mean(x), s=sd(x), c=length(x[which(x>0)]))})
return(summary)
}

#make individual map of shift
make_individual_shift_map <- function(){
namePDF <- paste("change/", spname, "_plot.pdf", sep="")
pdf(namePDF)
plot(shapefile)
points(spatial_points_time2, col=rgb(0, 0, 1, 0.15), pch="+")
points(time1_points, col=" dim grey", pch=".")
points(MeansTime1, col="black", pch="+", cex=2)
points(MeansTime2, col="dark blue", pch="o", cex=2)
points(centroid_time1, col="black", pch=16, cex=1)
points(centroid_time2, col="blue", pch=16, cex=1)
points(combined_centroid, col="red", pch=17, cex=2)
dev.off()
}

##This function displays just the points from the coordinates. For use in QGIS through LifeMapper
make_individual_shift_map_QGIS <- function(){
  ggplot(points(spatial_points_time2, col=rgb(0, 0, 1, 0.15), pch="+", xlim=c(-110, -70), ylim=c(25, 50)))
  points(time1_points, col=" dim grey", pch=".")
  points(MeansTime1, col="black", pch="+", cex=2)
  points(MeansTime2, col="dark blue", pch="o", cex=2)
  points(centroid_time1, col="black", pch=16, cex=1)
  points(centroid_time2, col="blue", pch=16, cex=1)
  points(combined_centroid, col="red", pch=17, cex=2)
}

## Make Polar Plot of a certain table, displaying a certain variable at a certain time (of time diff)
make_polar_plot <- function(data_climate, time){
  data_climate <- data.frame(data_climate)
  colnames(data_climate) <- c("species", "rotation", "distance")
  #data_climate[, "distance"] <- as.numeric(as.factor(data_climate[, "distance"]))
  data_climate[, "distance"] <- as.numeric(data_climate[, "distance"])
  data_climate[, "distance"] <- round(data_climate[, "distance"])
  #Here we create a table for the mean var per wedge for all species
  test <- summaryBy(distance~rotation, data_climate, FUN=function(x){c(m=mean(x), s=sd(x), c=length(x[which(x>0)]))})
  test[,"rotation"] <- as.numeric(as.character(test[,"rotation"]))
  test <- test[order(test[,"rotation"]),] 
  main <- paste("Mean distance from centroid to", time)
  plot <- polar.plot(as.numeric((test[,"distance.m"]/1000)),as.numeric(as.character(test[,"rotation"])),
                     main=main, clockwise=TRUE, lwd=3,line.col=4, start=90, 
                     rp.type="p", label.pos=seq(0,350, 30), labels=seq(0,350, 30))
  return (plot)
  }
