## load all the functions in memory
invisible(sapply(list.files(path = "R", pattern = "R$", full.names = TRUE), source))

folders <- list.dirs("data/")
folders <- folders[-1]

folder1 <- folders[1]
folder1 <- strsplit(folders, "//")
time1 <- folder1[[2]][2]
time2 <- folder1[[1]][2]
TimeDiff <- paste(time2, "-", time1)

#list all present files in time1 data folder
listfiles <- list_files(time1)
nspecies <- count_list(listfiles)

#make list of species names
namelist <- list_to_names(listfiles)

#set coordinates system
LATLON <- set_LATLON()

#set tables for later
TotalDiffDistTable <- matrix(data=NA, nrow=0, ncol=3)
TotalMeansTime1<- matrix(data=NA, nrow=0, ncol=4)
TotalMeansTime2 <- matrix(data=NA, nrow=0, ncol=4)
TotalDistTime1 <- matrix(data=NA, nrow = 0, ncol = 3)
TotalDistTime2 <- matrix(data=NA, nrow = 0, ncol = 3)
CoastDistTime1 <- matrix(data=NA, nrow=0, ncol=3)
CoastDistTime2 <- matrix(data=NA, nrow=0, ncol=3)
DistCoastDiffTable <- matrix(data=NA, nrow=0, ncol=3)

#make loop for future-present niche shift

for (i in 1:nspecies){# for each species in the folder
  
  spname <- namelist[i]
  
  time1_raster <- time_raster(time1)
  time1_points <- raster_to_points(time1_raster)
  nPoints1 <- number_points(time1_points)
  time1_points_table <- points_table(time1_points)
  
    
  if(nPoints1>1){
    centroid_time1 <- geomean(time1_points_table)
  }
  
  if (nPoints1 == 1)
  {
    centroid_time1 <- time1_points_table
  }
  
  time2_raster <- time_raster(time2)
  time2_points <- raster_to_points(time2_raster)
  nPoints2 <- number_points(time2_points)
  time2_points_table <- points_table(time2_points)


  ################################################################################
  #First calculate the distance to the mean 20 furthest points along each wedge. 
  ################################################################################
  
  
  ################################################################################
  #If there are some present cells in the time2 layer 
  ################################################################################
  if (nPoints2>1){
    time2_points_table <- points_table(time2_points)
    centroid_time2 <- geomean(time2_points_table)
 
    if(nPoints1==0){
      combined_centroid <- centroid_time2
    }
    
    if(nPoints1>0){
    combined_centroid<-cbind((centroid_time1[1]+centroid_time2[1])/2,(centroid_time1[2]+centroid_time2[2])/2)
    }
    
    spatial_points_time2 <- SpatialPoints(time2_points_table, proj4string=LATLON)
    spatial_n_time2 <- length(spatial_points_time2)
    
    
    
    table1 <- matrix(data=NA, nrow=0, ncol=1)
    
    for (d in 1:spatial_n_time2){
      dpt <- spatial_points_time2[d,1:2]
      bearing <- bearingRhumb(combined_centroid, dpt)
      table1 <- rbind(table1, bearing)
    }
    
    #create Final tables with no rows
    MeansTime2 <- matrix(data=NA, ncol=2, nrow=0)
    BearingTime2 <- matrix(data=NA, ncol=1, nrow=0)
    DistTime2 <- matrix(data=NA, ncol=1, nrow=0)
    
    
    #lets do the loop
    
    a <- seq(0,360,by=30)
    na <- length(a)-1
    
    
    for (b in 1:na){
      min <- a[b]
      max <- a[b+1]
      
      table1 <- matrix(as.numeric(table1))
      number <- subset(table1, table1[,1] > min & table1[,1] < max)  
      n <- length(number[,1])
      
      distWedge <- matrix(data=NA, ncol=1, nrow=0)
      matrixWedge <- matrix(data=NA, ncol=2, nrow=0)
      #look for original point coordinates to make distance matrix. 
      
      ################################################################################
      #CALCULATIONS FOR wedge IF several cells with 1 present
      
      if(n>1){
        distWedge <- get_20furthest_points_in_wedge(n, table1, time2, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, table1, time2, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
  
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
        
      } 
      
      ################################################################################
      #CALCULATIONS FOR wedge IF only 1 cell with 1 present
      
      if(n == 1){
        distWedge <- get_20furthest_points_in_wedge(n, table1, time2, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, table1, time2, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      ################################################################################
      #CALCULATIONS FOR wedge IF 0 only
      
      if(n<1){
        mean<- 0
        meanBearing<- 0
        meanDist<- 0
      }
      BearingTime2 <- rbind(BearingTime2, meanBearing)
      DistTime2 <- rbind(DistTime2, meanDist)
      MeansTime2 <- rbind(MeansTime2, mean)  
    }
    
    
    
    #SAME THING FOR PRESENT LAYER
    ################################################################################
    
    
    spatial_points_time1 <- SpatialPoints(time1_points_table, proj4string=LATLON)
    spatial_n_time1 <- length(spatial_points_time1)
    
    tableA <- matrix(data=NA, nrow=0, ncol=1)
    
    for (d in 1:spatial_n_time1){
      dpt <- spatial_points_time1[d,1:2]
      bearing <- bearingRhumb(combined_centroid, dpt)
      tableA <- rbind(tableA, bearing)
    }
    #create Final tables with no rows
    MeansTime1 <- matrix(data=NA, ncol=2, nrow=0)
    BearingTime1 <- matrix(data=NA, ncol=1, nrow=0)
    DistTime1 <- matrix(data=NA, ncol=1, nrow=0)
    
    #lets do the loop
    a <- seq(0,360,by=30)
    na <- length(a)-1
    
    
    for (d in 1:na){
      min <- a[d]
      max <- a[d+1]
      
      tableA <- matrix(as.numeric(tableA))
      number <- subset(tableA, tableA[,1] > min & tableA[,1] < max)  
      n <- length(number[,1])
      
      distWedge <- matrix(data=NA, ncol=1, nrow=0)
      matrixWedge <- matrix(data=NA, ncol=2, nrow=0)
      #look for original point coordinates to make distance matrix. 
      
      ################################################################################
      #CALCULATIONS FOR PRESENT IF more than 1 cell with value of 1
      
      if(n>1){
        distWedge <- get_20furthest_points_in_wedge(n, tableA, time1, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, tableA, time1, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      
      ################################################################################
      #CALCULATIONS FOR PRESENT IF only 1 cell with value of 1 present
      
      if(n == 1){
        distWedge <- get_20furthest_points_in_wedge(n, tableA, time1, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, tableA, time1, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      
      ################################################################################
      #CALCULATIONS FOR PRESENT IF 0 only
      
      if(n < 1){
        mean<- 0
        meanBearing<- 0
        meanDist<- 0
      }
      BearingTime1 <- rbind(BearingTime1, meanBearing)
      DistTime1 <- rbind(DistTime1, meanDist)
      MeansTime1 <- rbind(MeansTime1, mean)  
    }
    
    
    #NOW COMPARE BOTH TIMES
    ################################################################################
    
    DiffBearings <- BearingTime2 - BearingTime1
    DiffDist <- DistTime2 - DistTime1
    NewDiffDist <- cbind(spname,seq(30, 360, 30),DiffDist)
    NewMeansTime1 <- cbind(spname,seq(30, 360, 30),MeansTime1)
    NewMeansTime2 <- cbind(spname,seq(30, 360, 30),MeansTime2)
    NewDistTime1 <- cbind(spname, seq(30, 360, 30),DistTime1)
    NewDistTime2 <- cbind(spname, seq(30, 360, 30),DistTime2)
    
    ## make PDF individual map of shift
    shapefile <- readShapeLines("data/FLstate.shp", proj4string = LATLON)
    make_individual_shift_map()    
    ## display shift on map for QGIS
   # make_individual_shift_map_QGIS()
  }#close if loop of future layer contains 1+ cells with value 1
  
  
  
  ################################################################################
  #If there is 1 cell in the FUTURE layer 
  ################################################################################
  if (nPoints2 == 1){
    time2_points_table <- points_table(time2_points)
    centroid_time2 <- time2_points_table
    
    if(nPoints1==0){
      combined_centroid <- centroid_time2
    }
    
    if(nPoints1>0){
      combined_centroid<-cbind((centroid_time1[1]+centroid_time2[1])/2,(centroid_time1[2]+centroid_time2[2])/2)
    }
    
    
    
    spatial_points_time2 <- SpatialPoints(time2_points_table, proj4string=LATLON)
    spatial_n_time2 <- length(spatial_points_time2)
    
    
    
    table1 <- matrix(data=NA, nrow=0, ncol=1)
    
    for (d in 1:spatial_n_time2){
      dpt <- spatial_points_time2[d,1:2]
      #  bearing <- bearingRhumb(PresPastCentroid, dpt)
      bearing <- bearingRhumb(combined_centroid, dpt)
      table1 <- rbind(table1, bearing)
    }
    
    #create Final tables with no rows
    MeansTime2 <- matrix(data=NA, ncol=2, nrow=0)
    BearingTime2 <- matrix(data=NA, ncol=1, nrow=0)
    DistTime2 <- matrix(data=NA, ncol=1, nrow=0)
    
    
    #lets do the loop
    
    a <- seq(0,360,by=30)
    na <- length(a)-1
    
    
    for (b in 1:na){
      min <- a[b]
      max <- a[b+1]
      
      table1 <- matrix(as.numeric(table1))
      number <- subset(table1, table1[,1] > min & table1[,1] < max)  
      n <- length(number[,1])
      
      distWedge <- matrix(data=NA, ncol=1, nrow=0)
      matrixWedge <- matrix(data=NA, ncol=2, nrow=0)
      #look for original point coordinates to make distance matrix. 
      
      ################################################################################
      #CALCULATIONS FOR wedge IF several cells with 1 present
      
      if(n>1){
        distWedge <- get_20furthest_points_in_wedge(n, table1, time2, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, table1, time2, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
        
      } 
      
      ################################################################################
      #CALCULATIONS FOR wedge IF only 1 cell with 1 present
      
      if(n == 1){
        distWedge <- get_20furthest_points_in_wedge(n, table1, time2, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, table1, time2, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      ################################################################################
      #CALCULATIONS FOR wedge IF 0 only
      
      if(n<1){
        mean<- 0
        meanBearing<- 0
        meanDist<- 0
      }
      BearingTime2 <- rbind(BearingTime2, meanBearing)
      DistTime2 <- rbind(DistTime2, meanDist)
      MeansTime2 <- rbind(MeansTime2, mean)  
    }
    
    
    
    #SAME THING FOR PRESENT LAYER
    ################################################################################
    
    
    spatial_points_time1 <- SpatialPoints(time1_points_table, proj4string=LATLON)
    spatial_n_time1 <- length(spatial_points_time1)
    
    tableA <- matrix(data=NA, nrow=0, ncol=1)
    
    for (d in 1:spatial_n_time1){
      dpt <- time1_points[d,1:2]
      #  bearing <- bearingRhumb(PresPastCentroid, dpt)
      bearing <- bearingRhumb(combined_centroid, dpt)
      tableA <- rbind(tableA, bearing)
    }
    #create Final tables with no rows
    MeansTime1 <- matrix(data=NA, ncol=2, nrow=0)
    BearingTime1 <- matrix(data=NA, ncol=1, nrow=0)
    DistTime1 <- matrix(data=NA, ncol=1, nrow=0)
    
    #lets do the loop
    a <- seq(0,360,by=30)
    na <- length(a)-1
    
    
    for (d in 1:na){
      min <- a[d]
      max <- a[d+1]
      
      tableA <- matrix(as.numeric(tableA))
      number <- subset(tableA, tableA[,1] > min & tableA[,1] < max)  
      n <- length(number[,1])
      
      distWedge <- matrix(data=NA, ncol=1, nrow=0)
      matrixWedge <- matrix(data=NA, ncol=2, nrow=0)
      #look for original point coordinates to make distance matrix. 
      
      ################################################################################
      #CALCULATIONS FOR PRESENT IF more than 1 cell with value of 1
      
      if(n>1){
        distWedge <- get_20furthest_points_in_wedge(n, tableA, time1, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, tableA, time1, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      
      ################################################################################
      #CALCULATIONS FOR PRESENT IF only 1 cell with value of 1 present
      
      if(n == 1){
        distWedge <- get_20furthest_points_in_wedge(n, tableA, time1, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, tableA, time1, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      
      ################################################################################
      #CALCULATIONS FOR PRESENT IF 0 only
      
      if(n < 1){
        mean<- 0
        meanBearing<- 0
        meanDist<- 0
      }
      BearingTime1 <- rbind(BearingTime1, meanBearing)
      DistTime1 <- rbind(DistTime1, meanDist)
      MeansTime1 <- rbind(MeansTime1, mean)  
    }
    
    
    #NOW COMPARE BOTH TIMES
    ################################################################################
    
    DiffBearings <- BearingTime2 - BearingTime1
    DiffDist <- DistTime2 - DistTime1
    NewDiffDist <- cbind(spname,seq(30, 360, 30),DiffDist)
    NewMeansTime1 <- cbind(spname,seq(30, 360, 30),MeansTime1)
    NewMeansTime2 <- cbind(spname,seq(30, 360, 30),MeansTime2)
    NewDistTime1 <- cbind(spname, seq(30, 360, 30),DistTime1)
    NewDistTime2 <- cbind(spname, seq(30, 360, 30),DistTime2)
    
    make_individual_shift_map()   
  }#close if loop of future layer contains 1 cell with value 1
  
 
  ################################################################################
  #CALCULATIONS FOR FUTURE IF only 0 
  ################################################################################
  
  if (nPoints2<1){
    
    # Get results of 0 for future layer as it is empty
    ################################################################################
    
    MeansTime2 <- matrix(data=0, ncol=2, nrow=12)
    BearingTime2 <- matrix(data=0, ncol=1, nrow=12)
    DistTime2 <- matrix(data=0, ncol=1, nrow=12)
    
    combined_centroid <- centroid_time1
    
    #SAME THING FOR PRESENT LAYER
    ################################################################################
    
    spatial_points_time1 <- SpatialPoints(time1_points_table, proj4string=LATLON)
    spatial_n_time1 <- length(spatial_points_time1)
    
    table1 <- matrix(data=NA, nrow=0, ncol=1)
    
    for (d in 1:spatial_n_time1){
      dpt <- time1_points[d,1:2]
      #  bearing <- bearingRhumb(PresPastCentroid, dpt)
      bearing <- bearingRhumb(combined_centroid, dpt)
      table1 <- rbind(table1, bearing)
    }
    #create Final tables with no rows
    MeansTime1 <- matrix(data=NA, ncol=2, nrow=0)
    BearingTime1 <- matrix(data=NA, ncol=1, nrow=0)
    DistTime1 <- matrix(data=NA, ncol=1, nrow=0)
    
    #lets do the loop
    a <- seq(0,360,by=30)
    na <- length(a)-1
    
    
    for (b in 1:na){
      min <- a[b]
      max <- a[b+1]
      
      table1 <- matrix(as.numeric(table1))
      number <- subset(table1, table1[,1] > min & table1[,1] < max)  
      n <- length(number[,1])
      
      distWedge <- matrix(data=NA, ncol=1, nrow=0)
      matrixWedge <- matrix(data=NA, ncol=2, nrow=0)
      #look for original point coordinates to make distance matrix. 
      if(n>1){
        distWedge <- get_20furthest_points_in_wedge(n, tableA, time1, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, tableA, time1, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      if(n == 1){
        distWedge <- get_20furthest_points_in_wedge(n, tableA, time1, raster_to_points)
        matrixWedge <- get_20furthest_coordinates_in_wedge(n, tableA, time1, raster_to_points)
        furthest <- mean_dist_20furthest_points(distWedge)
        
        mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
        
        meanBearing <- bearingRhumb(combined_centroid, mean)
        meanDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
      } 
      
      
      if(n<1){
        mean<- combined_centroid
        meanBearing<- (min+max)/2
        meanDist<- 0
      }
      BearingTime1 <- rbind(BearingTime1, meanBearing)
      DistTime1 <- rbind(DistTime1, meanDist)
      MeansTime1 <- rbind(MeansTime1, mean)  
    }
    
    
    #NOW COMPARE BOTH TIMES
    ################################################################################
    
    DiffBearings <- BearingTime2 - BearingTime1
    DiffBearings <- as.numeric(DiffBearings)
    DiffDist <- DistTime2 - DistTime1
    DiffDist <- as.numeric(DiffDist)
    NewDiffDist <- cbind(spname,seq(30, 360, 30),DiffDist)
    NewMeansTime1 <- cbind(spname,seq(30, 360, 30),MeansTime1)
    NewMeansTime2 <- cbind(spname,seq(30, 360, 30),MeansTime2)
    NewDistTime1 <- cbind(spname, seq(30, 360, 30),DistTime1)
    NewDistTime2 <- cbind(spname, seq(30, 360, 30),DistTime2)
    
    make_individual_shift_map()
    
    }#close if loop of future layer contains 0 cells only
  
  ################################################################################
  #Lets summarize some of these results
  ################################################################################
  
  #This gets all the values from the loop into a single matrix FOR ALL SPECIES
  colnames(NewDiffDist) <- c("species", "rotation", "distance")
  TotalDiffDistTable<- rbind(TotalDiffDistTable,NewDiffDist)
  TotalMeansTime1<- rbind(TotalMeansTime1,NewMeansTime1)
  TotalMeansTime2 <- rbind(TotalMeansTime2, NewMeansTime2)
  TotalDistTime1 <- rbind(TotalDistTime1, NewDistTime1)
  TotalDistTime2 <- rbind(TotalDistTime2, NewDistTime2)
  colnames(TotalDistTime1) <- c("species", "rotation", "distance")
  colnames(TotalDistTime2) <- c("species", "rotation", "distance")
  row.names(TotalDistTime2) <- NULL
  row.names(TotalDistTime1) <- NULL
  TotalDistTime2 <- data.frame(TotalDistTime2)
  TotalDistTime1 <- data.frame(TotalDistTime1)
  ################################################################################
  #Calculate distance to coast
  
  DistToCoast <- matrix(data=NA, ncol=2, nrow=0)
  
  CoastPoints <- rasterToPoints(time1_raster, fun=function(time1_raster){time1_raster==1 | time1_raster==0})
  CoastRecs<-cbind(CoastPoints[,1],CoastPoints[,2])
  c <- nrow(CoastRecs)
  
  if (c>1){
    Coastcentroid <- geomean(CoastRecs)
  }
  
  if (c==1){
    Coastcentroid <- CoastRecs
  }
  
  Coast <- SpatialPoints(CoastRecs, proj4string=LATLON)
  nCoast <- length(Coast)
  
  table2 <- matrix(data=NA, nrow=0, ncol=1)
  
  for (d in 1:nCoast){
    dptCoast <- CoastPoints[d,1:2]
    bearingCoast <- bearingRhumb(combined_centroid, dptCoast)
    table2 <- rbind(table2, bearingCoast)
  }
  #create Final tables with no rows
  CoastMeans <- matrix(data=NA, ncol=2, nrow=0)
  CoastBearings <- matrix(data=NA, ncol=1, nrow=0)
  CoastDist <- matrix(data=NA, ncol=1, nrow=0)
  
  #lets do the loop
  a <- seq(0,360,by=30)
  na <- length(a)-1
  
   for (b in 1:na){
    min <- a[b]
    max <- a[b+1]
    
    table2 <- matrix(as.numeric(table2))
    number <- subset(table2, table2[,1] > min & table2[,1] < max)  
    n <- length(number[,1])
    
    distWedge <- matrix(data=NA, ncol=1, nrow=0)
    matrixWedge <- matrix(data=NA, ncol=2, nrow=0)
    #look for original point coordinates to make distance matrix. 
    if (n>1){
######
      distWedge <- get_20furthest_points_in_wedge(n, table2, time1, raster_to_points_coast)
      matrixWedge <- get_20furthest_coordinates_in_wedge(n, table2, time1, raster_to_points_coast)
      furthest <- mean_dist_20furthest_points(distWedge)
      
      mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
      
      CBearing <- bearingRhumb(combined_centroid, mean)
      CDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)
    }
    
    if (n==1){
      distWedge <- get_20furthest_points_in_wedge(n, table2, TimeDiff)
      matrixWedge <- get_20furthest_coordinates_in_wedge(n, table2, TimeDiff)
      furthest <- mean_dist_20furthest_points(distWedge)
      
      mean <- mean_points_table_for_wedge(furthest, distWedge, matrixWedge)
      
      CBearing <- bearingRhumb(combined_centroid, points)
      CDist <- distMeeus(combined_centroid, points, a=6378137, f=1/298.257223563)
    }

    if (n==0)
      
    {
      mean <- cbind(combined_centroid[,1], combined_centroid[,2])
      CBearing <- bearingRhumb(combined_centroid, mean)
      CDist <- distMeeus(combined_centroid, mean, a=6378137, f=1/298.257223563)    		
    }
    
    CoastBearings <- rbind(CoastBearings, CBearing)
    CoastDist <- rbind(CoastDist, CDist)
    CoastMeans <- rbind(CoastMeans, mean)  
  }
  
  DistToCoastPres <- CoastDist - DistTime1
  DistToCoastFut <- CoastDist - DistTime2
  
  #create table for Distance to Coast for All Species and All Wedges
  DistToCoastDiff <- DistToCoastPres - DistToCoastFut
  DistToCoastDiff <- cbind(spname, seq(30, 360, 30),DistToCoastDiff)
  DistCoastDiffTable <- rbind (DistCoastDiffTable, DistToCoastDiff)
  
  #Format DistTime1 and DistTime2 tables
  DistTime1 <- cbind(spname, seq(30, 360, 30), DistTime1)
  colnames(CoastDistTime1)<- c("species", "rotation", "distance")
  colnames(DistTime1)<- c("species", "rotation", "distance")
  DistTime2 <- cbind(spname, seq(30, 360, 30), DistTime2)
  colnames(CoastDistTime2)<- c("species", "rotation", "distance")
  colnames(DistTime2)<- c("species", "rotation", "distance")
  rownames(CoastDistTime1)<- NULL
  rownames(CoastDistTime2)<- NULL
  DistTime1 <- data.frame(DistTime1)
  DistTime2 <- data.frame(DistTime2)
  #add to All Species tables
  CoastDistTime1 <- rbind(CoastDistTime1, DistTime1)
  CoastDistTime2 <- rbind (CoastDistTime2, DistTime2)
  
  distance <- paste(CoastDistTime1[,3])
  distance <- as.numeric(distance)
  rotation <- paste(CoastDistTime1[,2])
  rotation <- as.numeric(rotation)

}  #close loop going through all species in folder


spname

###Now ditance from centroid for all species (mean) for time1


### create polar plot of distance from centroid for all species
name_time1 <- paste("change/MeanDistanceCentroidTo", time1, ".pdf", sep="")
pdf(name_time1)
plot1 <- make_polar_plot(TotalDistTime1, time1)
dev.off()

###Now ditance from centroid for all species (mean) for time2

name_time2 <- paste("change/MeanDistanceCentroidTo", time2, ".pdf", sep="")
pdf(name_time2)
plot3 <- make_polar_plot(TotalMeansTime2, time2)
dev.off()


####Now on with the difference in DISTANCE FROM CENTROID time2-time1  for all species (mean)
colnames(TotalDiffDistTable)<- c("species", "rotation", "MovDist")
rownames(TotalDiffDistTable) <- NULL
TotalDiffDistTable <- as.data.frame(TotalDiffDistTable)


pdf("change/ShiftOfEdge.pdf")
plot2 <- make_polar_plot(TotalDiffDistTable, TimeDiff)
dev.off()



write.table(DistCoastDiffTable, file="change/ShiftOfCoastDistance.csv", sep=",")
write.table(TotalMeansTime1, file="change/SpeciesMeanTime1.csv", sep=",")
write.table(TotalDiffDistTable, file="change/ShiftOfCentroid.csv", sep=",")
write.table(TotalMeansTime2, file="change/SpeciesMeanTime2.csv")



