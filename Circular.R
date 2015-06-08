#setwd("~/Desktop/models")
setwd("/scratch/lfs/cgermain/MAXENT/YEARLYBIOCLIM")
library(raster) 
library(sp)
library(rgeos)
library(geosphere)
library(fields)
library(plotrix)
library(maptools)
library(plyr)
library(reshape)
library(rgdal)
library(doBy)
library (eeptools)


ptm <- proc.time()
#list all present files
list <- list.files("PresentOut/BINARY_MAPS/", pattern="*threshold.asc")
nspecies <- length(list)
TotalDiffDist <- matrix(data=NA, nrow=0, ncol=3)
TotalAllMeans<- matrix(data=NA, nrow=0, ncol=4)
LATLON<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
flcoast <- readShapeLines("shapefiles/FLstate.shp",  proj4string=LATLON)

#list
nspecies

AllSpeciesDist <- matrix(data=NA, nrow=0, ncol=3)
FutureSpeciesDist <- matrix(data=NA, nrow=0, ncol=3)
DistCoastDiffTable <- matrix(data=NA, nrow=0, ncol=3)


#make loop for present-future niche shift

for (s in 1:nspecies){# for each species in the folder
#for (s in 465:470){
#s <- 465

presentFILE <- list[s]
parts <- strsplit(presentFILE, "_")
genus <- unlist(lapply(parts,"[",1))
species <- unlist(lapply(parts,"[",2))
spname <- paste(genus, species, sep="_")

spname

presentF <- paste("PresentOut/BINARY_MAPS/", presentFILE, sep="")
present <-raster(presentF)
PresPoints <- rasterToPoints(present, fun=function(present){present==1})
PresRecs<-cbind(PresPoints[,1],PresPoints[,2])
z <- nrow(PresRecs)


if(z>1){
	centroid <- geomean(PresRecs)
	}
		
if (z == 1)
   {
	centroid <- PresRecs
	}


futureFILE <- paste("Future70Out/Con7085/",spname, "_con.asc", sep="")
futureLayer <-raster(futureFILE)
futureLayer2 <- reclassify(futureLayer, c(-Inf,0.66,0,0.66,Inf,1))
futurePoints <- rasterToPoints(futureLayer2, fun=function(futureLayer2){futureLayer2==1})
nPoints <- nrow(futurePoints)
nPoints

################################################################################
#First calculate the distance to the mean 20 furthest points along each wedge. 
################################################################################


################################################################################
#If there are some present cells in the FUTURE layer 
################################################################################
if (nPoints>1){
futureRecs<-cbind(futurePoints[,1],futurePoints[,2])
futurecentroid <- geomean(futureRecs)

PresFutCentroid<-cbind((centroid[1]+futurecentroid[1])/2,(centroid[2]+futurecentroid[2])/2)



Future <- SpatialPoints(cbind(futurePoints[,1], futurePoints[,2]), proj4string=CRS("+proj=longlat +datum=WGS84"))
nFuture <- length(Future)



table1 <- matrix(data=NA, nrow=0, ncol=1)

for (d in 1:nFuture){
  dpt <- futurePoints[d,1:2]
#  bearing <- bearingRhumb(PresPastCentroid, dpt)
bearing <- bearingRhumb(PresFutCentroid, dpt)
  table1 <- rbind(table1, bearing)
					}

#create Final tables with no rows
futureMeans <- matrix(data=NA, ncol=2, nrow=0)
futureBearings <- matrix(data=NA, ncol=1, nrow=0)
futureDist <- matrix(data=NA, ncol=1, nrow=0)


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
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table1[,1]==x)
      coord <- futurePoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 		dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
      				} #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
  		}
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
	nP <- length(points)
#Here we add a condition in case n=1. In this case, geomean chokes and stops the script. If there is only 1 point in that wedge, it is its own mean. 
		 if(nP>1){
		    mean<-geomean(points)
			    }
    	if(nP<2) {
    		mean <- points
    		}
    		
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
    meanBearing <- bearingRhumb(PresFutCentroid, mean)
    meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)

  } 

################################################################################
#CALCULATIONS FOR wedge IF only 1 cell with 1 present

  if(n == 1){
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table1[,1]==x)
      coord <-  futurePoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 	  dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    				} #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    				}
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-points
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
 	meanBearing <- bearingRhumb(PresFutCentroid, mean)
 	meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  			} 

################################################################################
#CALCULATIONS FOR wedge IF 0 only

  if(n<1){
    mean<- 0
    meanBearing<- 0
    meanDist<- 0
  		}
futureBearings <- rbind(futureBearings, meanBearing)
futureDist <- rbind(futureDist, meanDist)
futureMeans <- rbind(futureMeans, mean)  
	}



#SAME THING FOR PRESENT LAYER
################################################################################


Pres <- SpatialPoints(cbind(PresPoints[,1], PresPoints[,2]), proj4string=CRS("+proj=longlat +datum=WGS84"))
nPres <- length(Pres)

tableA <- matrix(data=NA, nrow=0, ncol=1)

for (d in 1:nPres){
  dpt <- PresPoints[d,1:2]
#  bearing <- bearingRhumb(PresPastCentroid, dpt)
  bearing <- bearingRhumb(PresFutCentroid, dpt)
  tableA <- rbind(tableA, bearing)
					}
#create Final tables with no rows
AllMeans <- matrix(data=NA, ncol=2, nrow=0)
AllBearings <- matrix(data=NA, ncol=1, nrow=0)
AllDist <- matrix(data=NA, ncol=1, nrow=0)

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
    for (p in 1:n){
      x <- number[p,1]
      c <- which(tableA[,1]==x)
      coord <-  PresPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 	  dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    } #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    }
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-geomean(points)
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
 	meanBearing <- bearingRhumb(PresFutCentroid, mean)
 	meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  } 


################################################################################
#CALCULATIONS FOR PRESENT IF only 1 cell with value of 1 present
 
  if(n == 1){
    for (p in 1:n){
      x <- number[p,1]
      c <- which(tableA[,1]==x)
      coord <-  PresPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 	  dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    } #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    }
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-points
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
 	meanBearing <- bearingRhumb(PresFutCentroid, mean)
 	meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  } 


################################################################################
#CALCULATIONS FOR PRESENT IF 0 only
 
 if(n < 1){
    mean<- 0
    meanBearing<- 0
    meanDist<- 0
  }
  AllBearings <- rbind(AllBearings, meanBearing)
  AllDist <- rbind(AllDist, meanDist)
  AllMeans <- rbind(AllMeans, mean)  
}


#NOW COMPARE BOTH TIMES
################################################################################

DiffBearings <- AllBearings - futureBearings
DiffDist <- AllDist - futureDist
NewDiffDist <- cbind(spname,seq(30, 360, 30),DiffDist)
NewMeans <- cbind(spname,seq(30, 360, 30),AllMeans)

namePDF <- paste("Fut7085-Pres/", spname, "_plot.pdf", sep="")
#pdf("C6Plot.pdf")
pdf(namePDF)
plot(flcoast)
#map("state", "florida", myborder = 0)   
points(futurePoints, col=rgb(0, 0, 1, 0.15), pch="+")
points(PresPoints, col=" dim grey", pch=".")
points(AllMeans, col="black", pch="+", cex=2)
points(futureMeans, col="dark blue", pch="o", cex=2)
points(centroid, col="black", pch=16, cex=1)
points(futurecentroid, col="blue", pch=16, cex=1)
points(PresFutCentroid, col="red", pch=17, cex=2)
dev.off()

}#close if loop of future layer contains 1+ cells with value 1



################################################################################
#If there is 1 cell in the FUTURE layer 
################################################################################
if (nPoints == 1){
futureRecs<-cbind(futurePoints[,1],futurePoints[,2])
futurecentroid <- futureRecs

PresFutCentroid<-cbind((centroid[1]+futurecentroid[1])/2,(centroid[2]+futurecentroid[2])/2)



Future <- SpatialPoints(cbind(futurePoints[,1], futurePoints[,2]), proj4string=CRS("+proj=longlat +datum=WGS84"))
nFuture <- length(Future)



table1 <- matrix(data=NA, nrow=0, ncol=1)

for (d in 1:nFuture){
  dpt <- futurePoints[d,1:2]
#  bearing <- bearingRhumb(PresPastCentroid, dpt)
bearing <- bearingRhumb(PresFutCentroid, dpt)
  table1 <- rbind(table1, bearing)
				}

#create Final tables with no rows
futureMeans <- matrix(data=NA, ncol=2, nrow=0)
futureBearings <- matrix(data=NA, ncol=1, nrow=0)
futureDist <- matrix(data=NA, ncol=1, nrow=0)


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
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table1[,1]==x)
      coord <- futurePoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
       	      dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
					} #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
      		     }
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    nP <- length(points)
#Here we add a condition in case n=1. In this case, geomean chokes and stops the script. If there is only 1 point in that wedge, it is its own mean. 
      	  if(nP>1){
		    mean<-geomean(points)
				    }
    				    else {
    				    	 mean <- points
    					      }
						
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
    meanBearing <- bearingRhumb(PresFutCentroid, mean)
    meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)

  } 

################################################################################
#CALCULATIONS FOR wedge IF only 1 cell with 1 present

  if(n == 1){
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table1[,1]==x)
      coord <-  PresPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
         dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
      		   		    } #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
				}
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-points
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
     meanBearing <- bearingRhumb(PresFutCentroid, mean)
     meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
     	      	 } 

################################################################################
#CALCULATIONS FOR wedge IF 0 only

  if(n<1){
    mean<- 0
    meanBearing<- 0
    meanDist<- 0
    	       }
futureBearings <- rbind(futureBearings, meanBearing)
futureDist <- rbind(futureDist, meanDist)
futureMeans <- rbind(futureMeans, mean)  
	    }



#SAME THING FOR PRESENT LAYER
################################################################################


Pres <- SpatialPoints(cbind(PresPoints[,1], PresPoints[,2]), proj4string=CRS("+proj=longlat +datum=WGS84"))
nPres <- length(Pres)

tableA <- matrix(data=NA, nrow=0, ncol=1)

for (d in 1:nPres){
  dpt <- PresPoints[d,1:2]
#  bearing <- bearingRhumb(PresPastCentroid, dpt)
  bearing <- bearingRhumb(PresFutCentroid, dpt)
  tableA <- rbind(tableA, bearing)
				}
#create Final tables with no rows
AllMeans <- matrix(data=NA, ncol=2, nrow=0)
AllBearings <- matrix(data=NA, ncol=1, nrow=0)
AllDist <- matrix(data=NA, ncol=1, nrow=0)

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
    for (p in 1:n){
      x <- number[p,1]
      c <- which(tableA[,1]==x)
      coord <-  PresPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
         dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    } #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    }
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-geomean(points)
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
     meanBearing <- bearingRhumb(PresFutCentroid, mean)
     meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  } 


################################################################################
#CALCULATIONS FOR PRESENT IF only 1 cell with value of 1 present
 
  if(n == 1){
    for (p in 1:n){
      x <- number[p,1]
      c <- which(tableA[,1]==x)
      coord <-  PresPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
         dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    } #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    }
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-points
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
     meanBearing <- bearingRhumb(PresFutCentroid, mean)
     meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  } 


################################################################################
#CALCULATIONS FOR PRESENT IF 0 only
 
 if(n < 1){
    mean<- 0
    meanBearing<- 0
    meanDist<- 0
  }
  AllBearings <- rbind(AllBearings, meanBearing)
  AllDist <- rbind(AllDist, meanDist)
  AllMeans <- rbind(AllMeans, mean)  
}


#NOW COMPARE BOTH TIMES
################################################################################

DiffBearings <- AllBearings - futureBearings
DiffDist <- AllDist - futureDist
NewDiffDist <- cbind(spname,seq(30, 360, 30),DiffDist)
NewMeans <- cbind(spname,seq(30, 360, 30),AllMeans)

namePDF <- paste("Fut7085-Pres/",spname, "_plot.pdf", sep="")
#pdf("C6Plot.pdf")
pdf(namePDF)
plot(flcoast)
#map("state", "florida", myborder = 0)   
points(futurePoints, col=rgb(0, 0, 1, 0.15), pch="+")
points(PresPoints, col=" dim grey", pch=".")
points(AllMeans, col="black", pch="+", cex=2)
points(futureMeans, col="dark blue", pch="o", cex=2)
points(centroid, col="black", pch=16, cex=1)
points(futurecentroid, col="blue", pch=16, cex=1)
points(PresFutCentroid, col="red", pch=17, cex=2)
dev.off()

}#close if loop of future layer contains 1 cell with value 1






################################################################################
#CALCULATIONS FOR FUTURE IF only 0 
################################################################################


if (nPoints<1){

# Get results of 0 for future layer as it is empty
################################################################################

futureMeans <- matrix(data=0, ncol=2, nrow=12)
futureBearings <- matrix(data=0, ncol=1, nrow=12)
futureDist <- matrix(data=0, ncol=1, nrow=12)

PresFutCentroid <- centroid

#SAME THING FOR PRESENT LAYER
################################################################################

	Pres <- SpatialPoints(cbind(PresPoints[,1], PresPoints[,2]), proj4string=CRS("+proj=longlat +datum=WGS84"))
	nPres <- length(Pres)

	table1 <- matrix(data=NA, nrow=0, ncol=1)

	for (d in 1:nPres){
  		dpt <- PresPoints[d,1:2]
#  bearing <- bearingRhumb(PresPastCentroid, dpt)
  		bearing <- bearingRhumb(PresFutCentroid, dpt)
  		table1 <- rbind(table1, bearing)
					}
#create Final tables with no rows
AllMeans <- matrix(data=NA, ncol=2, nrow=0)
AllBearings <- matrix(data=NA, ncol=1, nrow=0)
AllDist <- matrix(data=NA, ncol=1, nrow=0)

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
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table1[,1]==x)
      coord <-  PresPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 	  dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    				} #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    				}
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-geomean(points)
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
 	meanBearing <- bearingRhumb(PresFutCentroid, mean)
 	meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  			} 

  if(n == 1){
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table1[,1]==x)
      coord <-  PresPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 	  dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    				} #close loop to gather all the points in distances in that wedge
    
    #now select the 20th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:20]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:20){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    				}
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-points
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
 	meanBearing <- bearingRhumb(PresFutCentroid, mean)
 	meanDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  			} 


  	if(n<1){
    		mean<- PresFutCentroid
    		meanBearing<- (min+max)/2
    		meanDist<- 0
  			}
  	AllBearings <- rbind(AllBearings, meanBearing)
  	AllDist <- rbind(AllDist, meanDist)
  	AllMeans <- rbind(AllMeans, mean)  
	}


#NOW COMPARE BOTH TIMES
################################################################################

DiffBearings <- AllBearings - futureBearings
DiffBearings <- as.numeric(DiffBearings)
DiffDist <- AllDist - futureDist
DiffDist <- as.numeric(DiffDist)
NewDiffDist <- cbind(spname,seq(30, 360, 30),DiffDist)
NewMeans <- cbind(spname,seq(30, 360, 30),AllMeans)


namePDF <- paste("Fut7085-Pres/", spname, "_plot.pdf", sep="")
#pdf("C6Plot.pdf")
pdf(namePDF)
plot(flcoast)
#map("state", "florida", myborder = 0)   
points(futurePoints, col=rgb(0, 0, 1, 0.15), pch="+")
points(PresPoints, col=" dim grey", pch=".")
points(AllMeans, col="black", pch="+", cex=2)
points(futureMeans, col="dark blue", pch="o", cex=2)
points(centroid, col="black", pch=16, cex=1)
points(futurecentroid, col="blue", pch=16, cex=1)
points(PresFutCentroid, col="red", pch=17, cex=2)
dev.off()
}#close if loop of future layer contains 0 cells only


################################################################################
#Lets summarize some of these results
################################################################################

#This gets all the values from the loop into a single matrix FOR ALL SPECIES

TotalDiffDist<- rbind(TotalDiffDist,NewDiffDist)
TotalAllMeans<- rbind(TotalAllMeans,NewMeans)


################################################################################
#Calculate distance to coast

DistToCoast <- matrix(data=NA, ncol=2, nrow=0)
 
CoastPoints <- rasterToPoints(present, fun=function(present){present==1 | present==0})
CoastRecs<-cbind(PresPoints[,1],PresPoints[,2])

c <- nrow(CoastRecs)

if (c>1){
Coastcentroid <- geomean(CoastRecs)
	      }

if (c==1){
Coastcentroid <- CoastRecs
	      }

Coast <- SpatialPoints(cbind(CoastPoints[,1], CoastPoints[,2]), proj4string=CRS("+proj=longlat +datum=WGS84"))
nCoast <- length(Coast)

table2 <- matrix(data=NA, nrow=0, ncol=1)

for (d in 1:nCoast){
  dptCoast <- CoastPoints[d,1:2]
#  bearing <- bearingRhumb(PresPastCentroid, dpt)
  bearingCoast <- bearingRhumb(PresFutCentroid, dptCoast)
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
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table2[,1]==x)
      coord <-  CoastPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 	  dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    			} #close loop to gather all the points in distances in that wedge
    
    #now select the 10th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:10]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:10){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    }
        
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean<-geomean(points)
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
 	CBearing <- bearingRhumb(PresFutCentroid, mean)
 	CDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)
  }
	
	 if (n==1){
    for (p in 1:n){
      x <- number[p,1]
      c <- which(table2[,1]==x)
      coord <-  CoastPoints[c,1:2]
      matrixWedge <- rbind(matrixWedge, coord)
      
#      dist<- distMeeus(PresPastCentroid, coord, a=6378137, f=1/298.257223563)
 	  dist<- distMeeus(PresFutCentroid, coord, a=6378137, f=1/298.257223563)
      distWedge <- rbind(distWedge, dist)
    			} #close loop to gather all the points in distances in that wedge
    
    #now select the 10th furthest points from the centroid and gather their coordinates to find centroid of those. 
    
    order <- distWedge[order(distWedge, decreasing=TRUE)]
    furthest <- order[1:10]
    
    furthestMatrix <- matrix(data=NA, ncol=2, nrow=0)
    for (x in 1:10){
      y <- which(distWedge==furthest[x])
      z <- matrixWedge[y,1:2]
      furthestMatrix <- rbind(furthestMatrix,z)
    }
        
    #meanmin is the centroid for the min-max wedge
    #name <- paste("mean", min, sep="")
    points <- cbind((furthestMatrix[,1]),(furthestMatrix[,2]))
    mean <- points
#    meanBearing <- bearingRhumb(PresPastCentroid, mean)
#    meanDist <- distMeeus(PresPastCentroid, mean, a=6378137, f=1/298.257223563)
 	CBearing <- bearingRhumb(PresFutCentroid, points)
 	CDist <- distMeeus(PresFutCentroid, points, a=6378137, f=1/298.257223563)
  }

	
	
	   if (n==0)
    
    	{
    	mean <- cbind(PresFutCentroid[,1], PresFutCentroid[,2])
    	 	CBearing <- bearingRhumb(PresFutCentroid, mean)
 			CDist <- distMeeus(PresFutCentroid, mean, a=6378137, f=1/298.257223563)    		
    	}
 
  CoastBearings <- rbind(CoastBearings, CBearing)
  CoastDist <- rbind(CoastDist, CDist)
  CoastMeans <- rbind(CoastMeans, mean)  

}



DistToCoastPres <- CoastDist - AllDist
DistToCoastFut <- CoastDist - futureDist

#create table for Distance to Coast for All Species and All Wedges
DistToCoastDiff <- DistToCoastPres - DistToCoastFut
DistToCoastDiff <- cbind(spname, seq(30, 360, 30),DistToCoastDiff)
DistCoastDiffTable <- rbind (DistCoastDiffTable, DistToCoastDiff)

#Format AllDist and futureDist tables
AllDist <- cbind(spname, seq(30, 360, 30), AllDist)
colnames(AllSpeciesDist)<- c("species", "rotation", "distance")
colnames(AllDist)<- c("species", "rotation", "distance")
futureDist <- cbind(spname, seq(30, 360, 30), futureDist)
colnames(FutureSpeciesDist)<- c("species", "rotation", "distance")
colnames(futureDist)<- c("species", "rotation", "distance")
rownames(AllSpeciesDist)<- NULL
rownames(FutureSpeciesDist)<- NULL
AllDist <- data.frame(AllDist)
futureDist <- data.frame(futureDist)
#add to All Species tables
AllSpeciesDist <- rbind(AllSpeciesDist, AllDist)
FutureSpeciesDist <- rbind (FutureSpeciesDist, futureDist)


distance <- paste(AllSpeciesDist[,3])
distance <- as.numeric(distance)
rotation <- paste(AllSpeciesDist[,2])
rotation <- as.numeric(rotation)



}  #close loop going through all species in folder

proc.time() - ptm  

spname

###Now ditance from centroid for all species (mean) FOR THE PRESENT
AllSpeciesDist$distance <- as.numeric(as.character(AllSpeciesDist$distance))
AllSpeciesDist$distance <- round(AllSpeciesDist$distance)
#MeanWedgeDist <- ddply(AllSpeciesDist, .(rotation), summarize, count=count(distance!=0), tot=sum(distance), mean=mean(distance), sd=sd(distance))
#Here we create a table for the mean distance per wedge for all species

test <- summaryBy(distance~rotation, data=AllSpeciesDist, FUN=function(x){c(m=mean(x), s=sd(x), c=length(x[which(x>0)]))})
test$rotation <- as.numeric(as.character(test$rotation))
test <- test[order(test$rotation),] 

### create polar plot of distance from centroid for all species
pdf("Fut7085-Pres/MeanDistanceFromCentroidForThePresent.pdf")
plot1 <- polar.plot(as.numeric((test$distance.m/1000)),as.numeric(as.character(test$rotation)),main="Mean Distance from centroid for the present", clockwise=TRUE, lwd=3,line.col=4, start=90, rp.type="p", label.pos=seq(0,350, 30), labels=seq(0,350, 30))
dev.off()

####Now on with the difference in DISTANCE FROM CENTROID FUTURE - PRESENT  for all species (mean)
colnames(TotalDiffDist)<- c("species", "rotation", "MovDist")
rownames(TotalDiffDist) <- NULL
TotalDiffDist <- as.data.frame(TotalDiffDist)
TotalDiffDist$MovDist<- as.numeric(as.character(TotalDiffDist$MovDist))
TotalDiffDist$MovDist<- round(TotalDiffDist$MovDist)


test2 <- summaryBy(MovDist~rotation, data=TotalDiffDist, FUN=function(x){c(m=mean(x), s=sd(x), c=length(x[which(x>0)]))})
test2$rotation <- as.numeric(as.character(test2$rotation))
test2 <- test2[order(test2$rotation),]

pdf("Fut7085-Pres/ShiftOfDistributionEdgeFuture-Present.pdf")
plot2 <- polar.plot(as.numeric((test2$MovDist.m/1000)),as.numeric(as.character(test2$rotation)),main="Shift of distribution edge: Future-Present", clockwise=TRUE, lwd=3,line.col=4, start=90, rp.type="p", label.pos=seq(0,350, 30), labels=seq(0,350, 30))
dev.off()

write.table(DistCoastDiffTable, file="Fut7085-Pres/Pres-FutCoastalDist.csv", sep=",")
write.table(TotalAllMeans, file="Fut7085-Pres/AllSpeciesMeans.csv", sep=",")
write.table(TotalDiffDist, file="Fut7085-Pres/Pres-FutCentroidDist.csv", sep=",")




