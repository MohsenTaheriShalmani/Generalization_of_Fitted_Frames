library(shapes)
library(rgl)
library(matlib)
library(RiemBase)
library(data.table)
library(Directional)
library(pracma)
library(numDeriv)
library(ggplot2)

#clear the environment
remove(list=ls())

#set working directory to file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load functions
if(TRUE){
  source("subFunctions/MathFunctions.R")
  source("subFunctions/meanFrames.R")
  source("subFunctions/normalsOfSkeletalSheetByTriangles.R")
  source("subFunctions/frameGenerator.R")
  source("subFunctions/rotateFrameForwardAndBackward.R")
}

#####################################################################################################
#####################################################################################################
# Read SPHARM-PDM meshes and ds-rep data 

#group labels
G1<-"CG"
G2<-"PD"

# read sreps
srepsDataG1<- read.csv(file=paste("../files/refined_17",G1,".csv",sep = ""), header=TRUE, sep=",")
srepsDataG2<- read.csv(file=paste("../files/refined_17",G2,".csv",sep = ""), header=TRUE, sep=",")
#View(srepsDataG1)
#View(srepsDataG2)

# # read sreps
# srepsDataG1<- read.csv(file=paste("../files/schizoData_CG.csv",sep = ""), header=TRUE, sep=",")
# srepsDataG2<- read.csv(file=paste("../files/schizoData_TG.csv",sep = ""), header=TRUE, sep=",")
# #View(srepsDataG1)
# #View(srepsDataG2)

#Read object meshes
meshPoints_G1 <- read.csv(file = paste("../files/meshPoints_17",G1,".csv",sep = ""),check.names = FALSE, header=TRUE, sep=",")
meshPoints_G2 <- read.csv(file = paste("../files/meshPoints_17",G2,".csv",sep = ""),check.names = FALSE, header=TRUE, sep=",")
#View(meshPoints_G1)
#View(meshPoints_G2)

# connections of triangular mesh
PolygonsCsv <- read.csv("../files/PolygonsCsv.csv")
polyMatrix<-cbind((PolygonsCsv$point1+1),(PolygonsCsv$point2+1),(PolygonsCsv$point3+1))

#Number of samples in each group
# nSamplesG1<-dim(meshPoints_G1)[2]
# nSamplesG1
# nSamplesG2<-dim(meshPoints_G2)[2]
# nSamplesG2
# or
nSamplesG1<-max(srepsDataG1$srepNumber)
nSamplesG1
nSamplesG2<-max(srepsDataG2$srepNumber)
nSamplesG2

#####################################################################################################
#####################################################################################################
# Extract number of srep GOPs

#No. samples and spokes
upSpoeksNumber<-max(srepsDataG1$SpokesNumber[which(srepsDataG1$srepNumber==1 & srepsDataG1$Spoke=='up')])
downSpoeksNumber<-max(srepsDataG1$SpokesNumber[which(srepsDataG1$srepNumber==1 & srepsDataG1$Spoke=='down')])
crestSpoksNumber<-max(srepsDataG1$SpokesNumber[which(srepsDataG1$srepNumber==1 & srepsDataG1$Spoke=='crest')])
nTotalRadii <- upSpoeksNumber + downSpoeksNumber + crestSpoksNumber
skelPointNo <- nTotalRadii-downSpoeksNumber
skelRange<-c(1:downSpoeksNumber,(2*downSpoeksNumber+1):nTotalRadii)

#####################################################################################################
#####################################################################################################
# Extract initial Skeletal and Boundary PDM of G1 and G2

source("subFunctions/readSrepsData.R")

tempG1<-readSrepsData(srepsData = srepsDataG1)
SkeletalPDMG1<-tempG1$SkeletalPDM
BoundaryPDMG1<-tempG1$BoundaryPDM
boundaryPlusSkeletal_G1<-tempG1$boundaryPlusSkeletal

tempG2<-readSrepsData(srepsData = srepsDataG2)
SkeletalPDMG2<-tempG2$SkeletalPDM
BoundaryPDMG2<-tempG2$BoundaryPDM
boundaryPlusSkeletal_G2<-tempG1$boundaryPlusSkeletal

#####################################################################################################
#####################################################################################################
# test shape correspondance between meshes and s-reps 
k<-1
k2<-0
for (i in srepsDataG1$file[which(srepsDataG1$SpokesNumber==1
                                 & srepsDataG1$Spoke=='crest')]) {
  if(i==names(meshPoints_G1)[k]){
    k2<-k2+1
  }else{
    cat("please find the correct list location for sample",i)
  }
  k<-k+1
}
if(k2==nSamplesG1){print("G1 group is OK. Please proceed!")}else{print("Please rearrenge the data! Do not proceed!!!")}

k<-1
k2<-0
for (i in srepsDataG2$file[which(srepsDataG2$SpokesNumber==1
                                 & srepsDataG2$Spoke=='crest')]) {
  if(i==names(meshPoints_G2)[k]){
    k2<-k2+1
  }else{
    cat("please find the correct list location for sample",i)
  }
  k<-k+1
}
if(k2==nSamplesG2){print("G2 group is OK. Please proceed!")}else{print("Please rearrenge the data! Do not proceed!!!")}


#####################################################################################################
#####################################################################################################
# Extract points of SPHARM-PDM meshes

numberOfPoints<-length(meshPoints_G1[,1])/3

spharmPDM_G1<-array(NA, dim = c(numberOfPoints,3,nSamplesG1))
k<-1
for (i in names(meshPoints_G1)){
  spharmPDM_G1[,,k]<-matrix(meshPoints_G1[[i]], ncol = 3, byrow = TRUE)
  k<-k+1
}

spharmPDM_G2<-array(NA, dim = c(numberOfPoints,3,nSamplesG2))
k<-1
for (i in names(meshPoints_G2)){
  spharmPDM_G2[,,k]<-matrix(meshPoints_G2[[i]], ncol = 3, byrow = TRUE)
  k<-k+1
}

#####################################################################################################
#####################################################################################################
# During model fitting some s-reps have correspondance issue that we fix it here.

# Fix correspondance issue among s-reps by fixFilipedSreps function
source("subFunctions/fixFlipIssue.R")

referenceSrep<-srepsDataG1[srepsDataG1$srepNumber==1,]

centerPoint<-16 #intrinsic centroid
p2<-13
p3<-17
v0<-BoundaryPDMG1[centerPoint,,1]-SkeletalPDMG1[centerPoint,,1]
reference_u0<-v0/sqrt(sum(v0^2))
v1<-SkeletalPDMG1[centerPoint,,1]-SkeletalPDMG1[p2,,1]
reference_u1<-v1/sqrt(sum(v1^2))
v2<-SkeletalPDMG1[centerPoint,,1]-SkeletalPDMG1[p3,,1]
reference_u2<-v2/sqrt(sum(v2^2))

# Initial grid size of SlicerSalt
gridSize = "5*7"

sorted_srepDataAdjustedG1<-fixFilipedSreps(srepsData = srepsDataG1,
                                           reference_u0 = reference_u0,
                                           reference_u1 = reference_u1,
                                           reference_u2 = reference_u2,
                                           centerPoint = centerPoint, p2 = p2, p3 = p3,
                                           gridSize = gridSize)
#View(sorted_srepDataAdjustedG1)
sorted_srepDataAdjustedG2<-fixFilipedSreps(srepsData = srepsDataG2,
                                           reference_u0 = reference_u0,
                                           reference_u1 = reference_u1,
                                           reference_u2 = reference_u2,
                                           centerPoint = centerPoint, p2 = p2, p3 = p3,
                                           gridSize = gridSize)
#View(sorted_srepDataAdjustedG2)


#####################################################################################################
#####################################################################################################
# Read adjusted Skeletal and Boundary PDM of G1 and G2

source("subFunctions/readSrepsData.R")

tempG1<-readSrepsData(srepsData = sorted_srepDataAdjustedG1)
SkeletalPDMG1<-tempG1$SkeletalPDM
BoundaryPDMG1<-tempG1$BoundaryPDM
boundaryPlusSkeletal_G1<-tempG1$boundaryPlusSkeletal

tempG2<-readSrepsData(srepsData = sorted_srepDataAdjustedG2)
SkeletalPDMG2<-tempG2$SkeletalPDM
BoundaryPDMG2<-tempG2$BoundaryPDM
boundaryPlusSkeletal_G2<-tempG2$boundaryPlusSkeletal

#####################################################################################################
#####################################################################################################
# Find spokes' directions and lengths (radii) after the correspondance adjustment

source("subFunctions/MathFunctions.R")

radii_G1<-array(NA, dim=c(nTotalRadii,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:nTotalRadii) {
    radii_G1[i,k]<-norm(boundaryPlusSkeletal_G1[i,,k]-
                          boundaryPlusSkeletal_G1[(i+nTotalRadii),,k],type = "2")
  }
}
radii_G2<-array(NA, dim=c(nTotalRadii,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:nTotalRadii) {
    radii_G2[i,k]<-norm(boundaryPlusSkeletal_G2[i,,k]-
                          boundaryPlusSkeletal_G2[(i+nTotalRadii),,k],type = "2")
  }
}

spokeDirections_G1<-array(NA, dim=c(nTotalRadii,3,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:nTotalRadii) {
    spokeDirections_G1[i,,k]<-convertVec2unitVec(BoundaryPDMG1[i,,k]-SkeletalPDMG1[i,,k])
  }
}
spokeDirections_G2<-array(NA, dim=c(nTotalRadii,3,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:nTotalRadii) {
    spokeDirections_G2[i,,k]<-convertVec2unitVec(BoundaryPDMG2[i,,k]-SkeletalPDMG2[i,,k])
  }
}

#####################################################################################################
#####################################################################################################
# Plot a ds-rep and its mesh from group1

sampleNo<-1 #choose sampleNo between 1 to nSamplesG1=108 to see other ds-reps
#plot
if(TRUE){
  open3d()
  srep1<-rbind(SkeletalPDMG1[,,sampleNo],BoundaryPDMG1[,,sampleNo])
  plot3d(SkeletalPDMG1[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  # plot mesh and normal vectors of a sample
  verts <- rbind(t(as.matrix(spharmPDM_G1[,,sampleNo])),1)
  trgls <- as.matrix(t(polyMatrix))
  tmesh <- tmesh3d(verts, trgls)
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
}
#####################################################################################################
#####################################################################################################
# Extra refinement
# Cut and stretch spokes so that their tips locate exactly at the boundary 
# We can skip this step because it does not significantly influence the result because
# it only affects the length of the spokes

# This step takes a few minutes!

source("subFunctions/cutAndStrechSpokes.R")
source("subFunctions/MathFunctions.R")

extraRefinement<-FALSE # set it to TRUE for extra refinemnt
# extraRefinement<-TRUE
if(extraRefinement){
  
  print("Start of extra refinement of Group 1.")
  tipOfCuttedSpokesG1<-array(NA,dim = dim(BoundaryPDMG1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    tipOfCuttedSpokesG1[,,i]<-cutAndStretchSpokes(BoundaryPDM = BoundaryPDMG1[,,i],
                                                   SkeletalPDM = SkeletalPDMG1[,,i],
                                                   meshPolygonData = PolygonsCsv,
                                                   spharmPDM = spharmPDM_G1[,,i])
    setTxtProgressBar(pb, i/nSamplesG1)
  }
  close(pb)
  print("Group 1 is done!")
  
  print("Start of extra refinement of Group 2.")
  tipOfCuttedSpokesG2<-array(NA,dim = dim(BoundaryPDMG2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    tipOfCuttedSpokesG2[,,i]<-cutAndStretchSpokes(BoundaryPDM = BoundaryPDMG2[,,i],
                                                   SkeletalPDM = SkeletalPDMG2[,,i],
                                                   meshPolygonData = PolygonsCsv,
                                                   spharmPDM = spharmPDM_G2[,,i])
    setTxtProgressBar(pb, i/nSamplesG2)
  }
  close(pb)
  print("Group 2 is done!")
  
  # Skip this part if you prefer s-reps without extra refinement
  BoundaryPDMG1<-tipOfCuttedSpokesG1
  BoundaryPDMG2<-tipOfCuttedSpokesG2
  
}

# Plot a extra refined ds-reps
sampleNo<-1 #choose sampleNo between 1 to nSamplesG1=108 to see other ds-reps
#plot
if(TRUE){
  open3d()
  srep1<-rbind(SkeletalPDMG1[,,sampleNo],BoundaryPDMG1[,,sampleNo])
  plot3d(SkeletalPDMG1[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 2,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  # plot mesh and normal vectors of a sample
  verts <- rbind(t(as.matrix(spharmPDM_G1[,,sampleNo])),1)
  trgls <- as.matrix(t(polyMatrix))
  tmesh <- tmesh3d(verts, trgls)
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
}

#####################################################################################################
#####################################################################################################
# Plot to see spoke correspondance
# We see spokes' tip and tails in seperated colors

#plot
if(TRUE){
  open3d()
  for (k in 1:nSamplesG2) {
    plot3d(SkeletalPDMG2[1:upSpoeksNumber,,k],type="p",col = "yellow",expand = 10,box=FALSE,add = TRUE)
    plot3d(BoundaryPDMG2[1:upSpoeksNumber,,k],type="p",col = "orange",expand = 10,box=FALSE,add = TRUE)
    plot3d(BoundaryPDMG2[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],type="p",col = "green",expand = 10,box=FALSE,add = TRUE)
    plot3d(SkeletalPDMG2[(2*upSpoeksNumber+1):nTotalRadii,,k],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
}

#####################################################################################################
#####################################################################################################
# Define lable of the frames for 5*7 grid of the skeletal sheet

# NB frame 16 is its own parent
framesCenters   <-c(16,13,10,7 ,4 ,1 ,2 ,3 ,19,22,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,29,30,26,27,23,24,20,21,17,18,14,15,11,12,8 ,9 ,5 ,6 ,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
framesParents   <-c(16,16,13,10,7 ,4 ,1 ,2 ,16,19,22,25,28,31,32,28,34,25,36,22,38,19,40,16,42,13,44,10,46,7 ,48,4 ,50,28,29,25,26,22,23,19,20,16,17,13,14,10,11,7 ,8 ,4 ,5 ,3 ,6 ,9 ,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51)
framesBackPoints<-c(13,16,13,10,7 ,4 ,1 ,2 ,16,19,22,25,28,31,32,28,34,25,36,22,38,19,40,16,42,13,44,10,46,7 ,48,4 ,50,28,29,25,26,22,23,19,20,16,17,13,14,10,11,7 ,8 ,4 ,5 ,3 ,6 ,9 ,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51)
framesFronts    <-c(19,10,7 ,4 ,1 ,2 ,3 ,52,22,25,28,31,32,33,62,35,63,37,64,39,65,41,66,43,67,45,68,47,69,49,70,51,71,30,61,27,60,24,59,21,58,18,57,15,56,12,55,9 ,54,6 ,53,rep(Inf,20)) #NB! crest frames don't have front point

sheetConnections<-c(2,5,8,11,14,17,20,23,26,29,32,34,36,38,40,42,44,46,48,50,2,3,6,9,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51,3,52:71,52)

# number of frames
numberOfFrames<-length(framesCenters)

#####################################################################################################
#####################################################################################################

# end of spokes at top and bottom parts
boundaryTop_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
for (k in 1:nSamplesG1) {
  boundaryTop_G1[,,k]<-rbind(BoundaryPDMG1[1:upSpoeksNumber,,k],BoundaryPDMG1[(2*upSpoeksNumber+1):nTotalRadii,,k])
}
boundaryTop_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
for (k in 1:nSamplesG2) {
  boundaryTop_G2[,,k]<-rbind(BoundaryPDMG2[1:upSpoeksNumber,,k],BoundaryPDMG2[(2*upSpoeksNumber+1):nTotalRadii,,k])
}
boundaryBottom_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
for (k in 1:nSamplesG1) {
  boundaryBottom_G1[,,k]<-rbind(BoundaryPDMG1[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],
                                BoundaryPDMG1[(2*upSpoeksNumber+1):nTotalRadii,,k])
}
boundaryBottom_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
for (k in 1:nSamplesG2) {
  boundaryBottom_G2[,,k]<-rbind(BoundaryPDMG2[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],
                                BoundaryPDMG2[(2*upSpoeksNumber+1):nTotalRadii,,k])
}

open3d()
for (i in 1:nSamplesG1) {
  plot3d(boundaryTop_G1[,,i],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryBottom_G1[,,i],type="p",col = "green",expand = 10,box=FALSE,add = TRUE)
}


#####################################################################################################
#####################################################################################################
# Calculate normal vectors of the skeletal sheet

# A general solution to calculate normals could be fitting a spline surface to the skeletal points.
# But spline fitting needs pre-alignment of skeletals to XY plane. Also it add variation to the normals.
# Here we generate normals based on the quadrilateral structure of the skeletal sheet.

#choose a method to find normals "sheetTriangles" or "splineFitting"
method2findNormals<-"sheetTriangles"
# method2findNormals<-"splineFitting"

if(method2findNormals=="sheetTriangles"){
  
  source("subFunctions/normalsOfSkeletalSheetByTriangles.R")
  
  skeletalSheet_G1<-SkeletalPDMG1[skelRange,,]
  skeletalSheet_G2<-SkeletalPDMG2[skelRange,,]
  medialNormals_G1<-array(NA,dim = dim(skeletalSheet_G1))
  normalsBoundaryTop_G1<-array(NA,dim = dim(skeletalSheet_G1))
  normalsBoundaryBottom_G1<-array(NA,dim = dim(skeletalSheet_G1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    medialNormals_G1[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = skeletalSheet_G1[,,i])
    normalsBoundaryTop_G1[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = boundaryTop_G1[,,i])
    normalsBoundaryBottom_G1[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = boundaryBottom_G1[,,i])
  }
  print("Group 1 is done!")
  medialNormals_G2<-array(NA,dim = dim(skeletalSheet_G2))
  normalsBoundaryTop_G2<-array(NA,dim = dim(skeletalSheet_G2))
  normalsBoundaryBottom_G2<-array(NA,dim = dim(skeletalSheet_G2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    medialNormals_G2[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = skeletalSheet_G2[,,i])
    normalsBoundaryTop_G2[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = boundaryTop_G2[,,i])
    normalsBoundaryBottom_G2[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = boundaryBottom_G2[,,i])
  }
  print("Group 2 is done!")
  
}else if(method2findNormals<-"splineFitting"){

  source("subFunctions/normalsOfSkeletalSheetBySpline.R")
  
  # Calculate normals by spline fitting. Takes a minute.
  # To fit spline, we use PCA + GPA to have skeletal sheets along the XY axis
  # Note that the GPA is not necessary, but we use it to make sure all the 
  # normals correspond to the northern of the ellipsoid
  sampleNo<-1 
  srep1<-rbind(SkeletalPDMG1[,,sampleNo],BoundaryPDMG1[,,sampleNo])
  centeredRefSrep<-prcomp(srep1,center = T)$x
  alignedByRefSrep_G1<-array(NA,dim = dim(boundaryPlusSkeletal_G1))
  for (i in 1:nSamplesG1) {
    alignedByRefSrep_G1[,,i]<-procOPA(centeredRefSrep,boundaryPlusSkeletal_G1[,,i],scale = F)$Bhat
  }
  alignedByRefSrep_G2<-array(NA,dim = dim(boundaryPlusSkeletal_G2))
  for (i in 1:nSamplesG2) {
    alignedByRefSrep_G2[,,i]<-procOPA(centeredRefSrep,boundaryPlusSkeletal_G2[,,i],scale = F)$Bhat
  }
  skeletalSheet_G1<-alignedByRefSrep_G1[skelRange,,]
  skeletalSheet_G2<-alignedByRefSrep_G2[skelRange,,]
  
  # Update skeletal PDMs
  SkeletalPDMG1<-alignedByRefSrep_G1[1:nTotalRadii,,]
  SkeletalPDMG2<-alignedByRefSrep_G2[1:nTotalRadii,,]
  
  # Update spokes directions
  # Since we applied PCA, we need to update spokes' directions in global coordinate system
  spokeDirections_G1<-array(NA, dim=c(nTotalRadii,3,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      spokeDirections_G1[i,,k]<-convertVec2unitVec(alignedByRefSrep_G1[i+nTotalRadii,,k]-alignedByRefSrep_G1[i,,k])
    }
  }
  spokeDirections_G2<-array(NA, dim=c(nTotalRadii,3,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      spokeDirections_G2[i,,k]<-convertVec2unitVec(alignedByRefSrep_G2[i+nTotalRadii,,k]-alignedByRefSrep_G2[i,,k])
    }
  }
  
  #calculate normals
  medialNormals_G1<-array(NA,dim = dim(skeletalSheet_G1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    temp<-normalsOfSkeletalSheetBySpline(centeredSkel = skeletalSheet_G1[,,i])
    medialNormals_G1[,,i]<-temp$medialNormals
  }
  close(pb)
  print("Group 1 is done!")
  medialNormals_G2<-array(NA,dim = dim(skeletalSheet_G2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    temp<-normalsOfSkeletalSheetBySpline(centeredSkel = skeletalSheet_G2[,,i])
    medialNormals_G2[,,i]<-temp$medialNormals
  }
  close(pb)
  print("Group 2 is done!")
  
}else{
  stop("Please specify method2findNormals!")
}

  
#####################################################################################################
#####################################################################################################
# Calculate frames vectors 

source("subFunctions/frameGenerator.R")

# frames in global coordinate system
if(TRUE){
  framesFirstVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesSecondVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesThirdVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesFirstVectorsTop_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesSecondVectorsTop_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesThirdVectorsTop_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesFirstVectorsBottom_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesSecondVectorsBottom_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesThirdVectorsBottom_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    temp<-frameGenerator(centeredSkel = skeletalSheet_G1[,,i],medialNormals = medialNormals_G1[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectors_G1[,,i]<-temp$framesFirstVec
    framesSecondVectors_G1[,,i]<-temp$framesSecondVec
    framesThirdVectors_G1[,,i]<-temp$framesThirdVec
    
    temp1<-frameGenerator(centeredSkel = boundaryTop_G1[,,i],medialNormals = normalsBoundaryTop_G1[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectorsTop_G1[,,i]<-(-temp1$framesFirstVec)
    framesSecondVectorsTop_G1[,,i]<-temp1$framesSecondVec
    framesThirdVectorsTop_G1[,,i]<-temp1$framesThirdVec
    
    temp2<-frameGenerator(centeredSkel = boundaryBottom_G1[,,i],medialNormals = normalsBoundaryBottom_G1[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectorsBottom_G1[,,i]<-temp2$framesFirstVec
    framesSecondVectorsBottom_G1[,,i]<-temp2$framesSecondVec
    framesThirdVectorsBottom_G1[,,i]<-temp2$framesThirdVec
    
  }
  close(pb)
  print("Group 1 is done!")
  framesFirstVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesSecondVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesThirdVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesFirstVectorsTop_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesSecondVectorsTop_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesThirdVectorsTop_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesFirstVectorsBottom_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesSecondVectorsBottom_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesThirdVectorsBottom_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    temp<-frameGenerator(centeredSkel = skeletalSheet_G2[,,i],medialNormals = medialNormals_G2[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectors_G2[,,i]<-temp$framesFirstVec
    framesSecondVectors_G2[,,i]<-temp$framesSecondVec
    framesThirdVectors_G2[,,i]<-temp$framesThirdVec
    
    temp1<-frameGenerator(centeredSkel = boundaryTop_G2[,,i],medialNormals = normalsBoundaryTop_G2[,,i],
                          framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectorsTop_G2[,,i]<-(-temp1$framesFirstVec)
    framesSecondVectorsTop_G2[,,i]<-temp1$framesSecondVec
    framesThirdVectorsTop_G2[,,i]<-temp1$framesThirdVec
    
    temp2<-frameGenerator(centeredSkel = boundaryBottom_G2[,,i],medialNormals = normalsBoundaryBottom_G2[,,i],
                          framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectorsBottom_G2[,,i]<-temp2$framesFirstVec
    framesSecondVectorsBottom_G2[,,i]<-temp2$framesSecondVec
    framesThirdVectorsBottom_G2[,,i]<-temp2$framesThirdVec
  }
  close(pb)
  print("Group 2 is done!")
}

#####################################################################################################
#####################################################################################################
# Plot LP-ds-rep

sampleNo<-1 #choose sampleNo between 1 to nSamplesG1=108 to see other ds-reps

#plot
if(TRUE){
  open3d()
  for (i in 2:numberOfFrames) {
    plot3d(rbind(skeletalSheet_G1[framesCenters[i],,sampleNo],skeletalSheet_G1[framesParents[i],,sampleNo])
           ,type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(boundaryTop_G1[framesCenters[i],,sampleNo],boundaryTop_G1[framesParents[i],,sampleNo])
           ,type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(boundaryBottom_G1[framesCenters[i],,sampleNo],boundaryBottom_G1[framesParents[i],,sampleNo])
           ,type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    
  }
  plot3d(skeletalSheet_G1[sheetConnections,,sampleNo],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryTop_G1[sheetConnections,,sampleNo],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryBottom_G1[sheetConnections,,sampleNo],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in framesCenters) {
    vectors3d(skeletalSheet_G1[i,,sampleNo]+framesFirstVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(skeletalSheet_G1[i,,sampleNo]+framesSecondVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(skeletalSheet_G1[i,,sampleNo]+framesThirdVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  for (i in framesCenters) {
    vectors3d(boundaryTop_G1[i,,sampleNo]+framesFirstVectorsTop_G1[i,,sampleNo],origin = boundaryTop_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryTop_G1[i,,sampleNo]+framesSecondVectorsTop_G1[i,,sampleNo],origin = boundaryTop_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryTop_G1[i,,sampleNo]+framesThirdVectorsTop_G1[i,,sampleNo],origin = boundaryTop_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  for (i in framesCenters) {
    vectors3d(boundaryBottom_G1[i,,sampleNo]+framesFirstVectorsBottom_G1[i,,sampleNo],origin = boundaryBottom_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryBottom_G1[i,,sampleNo]+framesSecondVectorsBottom_G1[i,,sampleNo],origin = boundaryBottom_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryBottom_G1[i,,sampleNo]+framesThirdVectorsBottom_G1[i,,sampleNo],origin = boundaryBottom_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  for (i in 1:nTotalRadii) {
    vectors3d(SkeletalPDMG1[i,,sampleNo]+(spokeDirections_G1[i,,sampleNo]*radii_G1[i,sampleNo]),origin = SkeletalPDMG1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="red", lwd=1)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
}

#####################################################################################################
#####################################################################################################
# Combine frames vectors to make SO(3) frames in global coordinate system

if(TRUE){
  frames_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  framesTop_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  framesBottom_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in framesCenters) {
      frames_G1[,,i,k]<- rbind(framesFirstVectors_G1[i,,k],
                               framesSecondVectors_G1[i,,k],
                               framesThirdVectors_G1[i,,k])
      framesTop_G1[,,i,k]<- rbind(framesFirstVectorsTop_G1[i,,k],
                               framesSecondVectorsTop_G1[i,,k],
                               framesThirdVectorsTop_G1[i,,k])
      framesBottom_G1[,,i,k]<- rbind(framesFirstVectorsBottom_G1[i,,k],
                               framesSecondVectorsBottom_G1[i,,k],
                               framesThirdVectorsBottom_G1[i,,k])
    }
  }
  frames_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  framesTop_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  framesBottom_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in framesCenters) {
      frames_G2[,,i,k]<- rbind(framesFirstVectors_G2[i,,k],
                               framesSecondVectors_G2[i,,k],
                               framesThirdVectors_G2[i,,k])
      framesTop_G2[,,i,k]<- rbind(framesFirstVectorsTop_G2[i,,k],
                                  framesSecondVectorsTop_G2[i,,k],
                                  framesThirdVectorsTop_G2[i,,k])
      framesBottom_G2[,,i,k]<- rbind(framesFirstVectorsBottom_G2[i,,k],
                                     framesSecondVectorsBottom_G2[i,,k],
                                     framesThirdVectorsBottom_G2[i,,k])
    }
  }
}

#####################################################################################################
#####################################################################################################
# Calculate children frames coordinates based on their parents frames

# redefine parents relations based on Prof. Pizer's article
# NB frame 16 is its own parent
framesCenters <-c(16,13,10,7 ,4 ,1 ,2 ,3 ,19,22,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,29,30,26,27,23,24,20,21,17,18,14,15,11,12,8 ,9 ,5 ,6 ,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
framesParents2 <-c(16,16,16,16,16,4 ,4 ,4 ,16,16,16,16,28,28,28,28,28,25,25,22,22,19,19,16,16,13,13,10,10,7 ,7 ,4 ,4 ,28,28,25,25,22,22,19,19,16,16,13,13,10,10,7 ,7 ,4 ,4 ,4 ,4 ,7 ,10,13,16,19,22,25,28,28,28,25,22,19,16,13,10,7 ,4 )
sampleNo<-1 #choose sampleNo between 1 to nSamplesG1=108 to see other ds-reps


source("subFunctions/rotateFrameForwardAndBackward.R")

if(TRUE){
  framesBasedOnParents_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  framesBasedOnParentsTop_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  framesBasedOnParentsBottom_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents2[i]
      framesBasedOnParents_G1[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                             vectors2rotate = frames_G1[,,k1,k])
      
      framesBasedOnParentsTop_G1[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k1,k],
                                                             vectors2rotate = framesTop_G1[,,k1,k])
      framesBasedOnParentsBottom_G1[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k1,k],
                                                                vectors2rotate = framesBottom_G1[,,k1,k])
    } 
  }
  print("Group 1 is done!")
  framesBasedOnParents_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  framesBasedOnParentsTop_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  framesBasedOnParentsBottom_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents2[i]
      framesBasedOnParents_G2[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                             vectors2rotate = frames_G2[,,k1,k])
      
      framesBasedOnParentsTop_G2[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k1,k],
                                                                vectors2rotate = framesTop_G2[,,k1,k])
      framesBasedOnParentsBottom_G2[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k1,k],
                                                                   vectors2rotate = framesBottom_G2[,,k1,k])
    } 
  }
  print("Group 2 is done!")
}

#####################################################################################################
#####################################################################################################
# Calculate spokes directions based on their frames

if(TRUE){
  spokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesDirectionsBasedOnFrames_G1[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,frameOfSpokeNo,k],
                                                                    vectors2rotate = spokeDirections_G1[i,,k])
      
    }
  }
  print("Group 1 is done!")
  spokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesDirectionsBasedOnFrames_G2[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,frameOfSpokeNo,k],
                                                                    vectors2rotate = spokeDirections_G2[i,,k])
      
      
    }
  }
  print("Group 2 is done!")
}

#####################################################################################################
#####################################################################################################
# Calculate connection lengths and directions based on their frames

if(TRUE){
  
  connections_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  connectionsLengths_G1<-array(NA,dim=c(numberOfFrames,nSamplesG1))
  connectionsBasedOnParentFrames_G1<-array(NA,dim = dim(connections_G1))
  for (k in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents2[i]
      tempVec<-skeletalSheet_G1[k1,,k]- skeletalSheet_G1[k2,,k]
      
      connectionsLengths_G1[k1,k]<-norm(tempVec,type = "2")
      
      if(norm(tempVec,type = "2")==0){
        connections_G1[k1,,k]<-c(0,0,0)
      }else{
        connections_G1[k1,,k]<-convertVec2unitVec(tempVec)
      }
      
      connectionsBasedOnParentFrames_G1[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                                      vectors2rotate = connections_G1[k1,,k])
    }
  }
  print("Group 1 is done!")
  connections_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  connectionsLengths_G2<-array(NA,dim=c(numberOfFrames,nSamplesG2))
  connectionsBasedOnParentFrames_G2<-array(NA,dim = dim(connections_G2))
  for (k in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents2[i]
      tempVec<-skeletalSheet_G2[k1,,k]- skeletalSheet_G2[k2,,k]
      
      connectionsLengths_G2[k1,k]<-norm(tempVec,type = "2")
      
      if(norm(tempVec,type = "2")==0){
        connections_G2[k1,,k]<-c(0,0,0)
      }else{
        connections_G2[k1,,k]<-convertVec2unitVec(tempVec)
      }
      
      connectionsBasedOnParentFrames_G2[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                                      vectors2rotate = connections_G2[k1,,k])
    }
  }
  print("Group 2 is done!")
}


#####################################################################################################
#####################################################################################################
# calculate LP-sizes

# LP sizes
LP_sizes_G1<-rep(NA,nSamplesG1)
for (i in 1:nSamplesG1) {
  LP_sizes_G1[i]<-sum(connectionsLengths_G1[,i])+sum(radii_G1[,i])
}
LP_sizes_G2<-rep(NA,nSamplesG2)
for (i in 1:nSamplesG2) {
  LP_sizes_G2[i]<-sum(connectionsLengths_G2[,i])+sum(radii_G2[,i])
}

#####################################################################################################
#####################################################################################################
#Removing or preserving the scale by LP-size

#choose the type of study "sizeAndShapeAnalysis" or "shapeAnalysis"

# typeOfStudy<-"shapeAnalysis"              #removing scale
typeOfStudy<-"sizeAndShapeAnalysis"     #preserving scale

if(typeOfStudy=="sizeAndShapeAnalysis"){
  sizes_G1<-rep(1,nSamplesG1)
  sizes_G2<-rep(1,nSamplesG2)
  
  radiiScaled_G1<-radii_G1 #we don't have scaling in size-and-shape analysis
  radiiScaled_G2<-radii_G2
  
  connectionsLengthsScaled_G1<-connectionsLengths_G1
  connectionsLengthsScaled_G2<-connectionsLengths_G2
  
}else if(typeOfStudy=="shapeAnalysis"){
  sizes_G1<-rep(NA,nSamplesG1)
  for (i in 1:nSamplesG1) {
    sizes_G1[i]<-sum(connectionsLengths_G1[,i])+sum(radii_G1[,i])
  }
  sizes_G2<-rep(NA,nSamplesG2)
  for (i in 1:nSamplesG2) {
    sizes_G2[i]<-sum(connectionsLengths_G2[,i])+sum(radii_G2[,i])
  }
  
  radiiScaled_G1<-sweep(radii_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  radiiScaled_G2<-sweep(radii_G2, 2, sizes_G2, "/") 
  
  connectionsLengthsScaled_G1<-sweep(connectionsLengths_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  connectionsLengthsScaled_G2<-sweep(connectionsLengths_G2, 2, sizes_G2, "/") 
  
}

#####################################################################################################
#####################################################################################################
# Calculate mean LP-ds-rep

# choose the method to calculate mean directions 
# PNS takes a few mintes!
typeOfMeanDirection<-"Frechet"
# typeOfMeanDirection<-"PNS"

# IF Boost==TRUE we start with initial frame based on mean n and b otherwise we start from the aligned centroids
Boost<-TRUE
plot4GradientDescent<-FALSE  # turn off plotting by plot4GradientDescent<-FALSE


#calculate mean frames in local and global coordinate systems
if(TRUE){
  source("subFunctions/meanFrames.R")
  
  meanFramesBasedOnParents_G1<-array(NA, dim = c(3,3,numberOfFrames))
  meanFramesBasedOnParentsTop_G1<-array(NA, dim = c(3,3,numberOfFrames))
  meanFramesBasedOnParentsBottom_G1<-array(NA, dim = c(3,3,numberOfFrames))
  for (k in framesCenters) {
    # readline(prompt="Press [enter] to continue")
    cat("Frame",k,"is done! \n")
    if(k==16){
      meanFramesBasedOnParents_G1[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
    }else{
      meanFramesBasedOnParents_G1[,,k]<-gradientDescent4meanFrame(framesBasedOnParents_G1[,,k,],
                                                                  Boost=TRUE,
                                                                  stepSize=0.001,
                                                                  threshold=1e-01,
                                                                  whileLimit=5000,
                                                                  method=typeOfMeanDirection,
                                                                  plotting=plot4GradientDescent)
    }
    meanFramesBasedOnParentsTop_G1[,,k]<-gradientDescent4meanFrame(framesBasedOnParentsTop_G1[,,k,],
                                                                   Boost=TRUE,
                                                                   stepSize=0.001,
                                                                   threshold=1e-01,
                                                                   whileLimit=5000,
                                                                   method=typeOfMeanDirection,
                                                                   plotting=plot4GradientDescent)
    meanFramesBasedOnParentsBottom_G1[,,k]<-gradientDescent4meanFrame(framesBasedOnParentsBottom_G1[,,k,],
                                                                      Boost=TRUE,
                                                                      stepSize=0.001,
                                                                      threshold=1e-01,
                                                                      whileLimit=5000,
                                                                      method=typeOfMeanDirection,
                                                                      plotting=plot4GradientDescent)
  }
  if(plot4GradientDescent==TRUE){
    open3d()
  }
  meanFramesBasedOnParents_G2<-array(NA, dim = c(3,3,numberOfFrames))
  meanFramesBasedOnParentsTop_G2<-array(NA, dim = c(3,3,numberOfFrames))
  meanFramesBasedOnParentsBottom_G2<-array(NA, dim = c(3,3,numberOfFrames))
  for (k in framesCenters) {
    if(k==16){
      meanFramesBasedOnParents_G2[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
    }else{
      meanFramesBasedOnParents_G2[,,k]<-gradientDescent4meanFrame(framesBasedOnParents_G2[,,k,],
                                                                  Boost=TRUE,
                                                                  stepSize=0.001,
                                                                  threshold=1e-02,
                                                                  whileLimit=5000,
                                                                  method=typeOfMeanDirection,
                                                                  plotting=plot4GradientDescent)
    }
    meanFramesBasedOnParentsTop_G2[,,k]<-gradientDescent4meanFrame(framesBasedOnParentsTop_G2[,,k,],
                                                                   Boost=TRUE,
                                                                   stepSize=0.001,
                                                                   threshold=1e-01,
                                                                   whileLimit=5000,
                                                                   method=typeOfMeanDirection,
                                                                   plotting=plot4GradientDescent)
    meanFramesBasedOnParentsBottom_G2[,,k]<-gradientDescent4meanFrame(framesBasedOnParentsBottom_G2[,,k,],
                                                                      Boost=TRUE,
                                                                      stepSize=0.001,
                                                                      threshold=1e-01,
                                                                      whileLimit=5000,
                                                                      method=typeOfMeanDirection,
                                                                      plotting=plot4GradientDescent)
  }
  
  meanFramesGlobalCoordinate_G1<-array(NA,dim = dim(meanFramesBasedOnParents_G1))
  meanFramesGlobalCoordinate_G1[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  meanFramesTopGlobalCoordinate_G1<-array(NA,dim = dim(meanFramesBasedOnParentsTop_G1))
  meanFramesBottomGlobalCoordinate_G1<-array(NA,dim = dim(meanFramesBasedOnParentsBottom_G1))
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents2[k]
    child_Index<-framesCenters[k]
    updatedParent<-meanFramesGlobalCoordinate_G1[,,parent_Index]
    meanFramesGlobalCoordinate_G1[,,child_Index]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = meanFramesBasedOnParents_G1[,,child_Index])
  }
  for (k in 1:numberOfFrames) {
    meanFramesTopGlobalCoordinate_G1[,,k]<-rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,k],
                                                                              vectorsInMainAxes = meanFramesBasedOnParentsTop_G1[,,k])
    meanFramesBottomGlobalCoordinate_G1[,,k]<-rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,k],
                                                                                 vectorsInMainAxes = meanFramesBasedOnParentsBottom_G1[,,k])
  }
  meanFramesGlobalCoordinate_G2<-array(NA,dim = dim(meanFramesBasedOnParents_G2))
  meanFramesGlobalCoordinate_G2[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  meanFramesTopGlobalCoordinate_G2<-array(NA,dim = dim(meanFramesBasedOnParentsTop_G2))
  meanFramesBottomGlobalCoordinate_G2<-array(NA,dim = dim(meanFramesBasedOnParentsBottom_G2))
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents2[k]
    child_Index<-framesCenters[k]
    updatedParent<-meanFramesGlobalCoordinate_G2[,,parent_Index]
    meanFramesGlobalCoordinate_G2[,,child_Index]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = meanFramesBasedOnParents_G2[,,child_Index])
  }
  for (k in 1:numberOfFrames) {
    meanFramesTopGlobalCoordinate_G2[,,k]<-rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,k],
                                                                              vectorsInMainAxes = meanFramesBasedOnParentsTop_G2[,,k])
    meanFramesBottomGlobalCoordinate_G2[,,k]<-rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,k],
                                                                                 vectorsInMainAxes = meanFramesBasedOnParentsBottom_G2[,,k])
  }
  print("Done!")
}

# Calculate mean spokes' directions based on frames
if(TRUE){
  meanSpokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G1[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G1[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G1[i,,]) 
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G1[i,,])
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-pns(spokesDirectionsBasedOnFrames_G1[i,,],sphere.type = sphereType)$PNS$mean 
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet")
    }
  }
  meanSpokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G2[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G2[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G2[i,,])
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G2[i,,])
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-pns(spokesDirectionsBasedOnFrames_G2[i,,],sphere.type = sphereType)$PNS$mean 
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet")
    }
  }
  print("Done!")
}

# Calculate mean spokes' directions based on global coordinate (using mean frames in global coordinate)
if(TRUE){
  meanSpokesDirectionsGlobalCoordinate_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    spokeNo<-i # 1<=spokeNo<=nTotalRadii
    frameOfSpokeNo<-NA
    if(spokeNo<=upSpoeksNumber){
      frameOfSpokeNo<-spokeNo
    }else if(spokeNo<=2*upSpoeksNumber){
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }
    meanSpokesDirectionsGlobalCoordinate_G1[i,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,frameOfSpokeNo],
                                         vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G1[i,])
    
  }
  meanSpokesDirectionsGlobalCoordinate_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    spokeNo<-i # 1<=spokeNo<=nTotalRadii
    frameOfSpokeNo<-NA
    if(spokeNo<=upSpoeksNumber){
      frameOfSpokeNo<-spokeNo
    }else if(spokeNo<=2*upSpoeksNumber){
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }
    meanSpokesDirectionsGlobalCoordinate_G2[i,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,frameOfSpokeNo],
                                         vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G2[i,])
    
  }
}

# Calculate geometric mean of spokes' lengths
if(TRUE){
  radiiMean_G1<-exp(rowMeans(log(radiiScaled_G1)))
  radiiMean_G2<-exp(rowMeans(log(radiiScaled_G2)))
}

# Calculate geometric mean of connections' lengths
if(TRUE){
  meanConnectionsLengths_G1<-exp(rowMeans(log(connectionsLengthsScaled_G1)))
  meanConnectionsLengths_G2<-exp(rowMeans(log(connectionsLengthsScaled_G2)))
}

# Calculate mean connection directions based on frames
if(TRUE){
  meanConnectionsBasedOnParentFrames_G1<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    if(i==framesCenters[1]){
      meanConnectionsBasedOnParentFrames_G1[i,]<-c(0,0,0)
    }else{
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G1[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        meanConnectionsBasedOnParentFrames_G1[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G1[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        meanConnectionsBasedOnParentFrames_G1[i,]<-frechetMean(connectionsBasedOnParentFrames_G1[i,,]) 
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G1[i,,])
        meanConnectionsBasedOnParentFrames_G1[i,]<-pns(connectionsBasedOnParentFrames_G1[i,,],sphere.type = sphereType)$PNS$mean  
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
  }
  meanConnectionsBasedOnParentFrames_G2<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    if(i==framesCenters[1]){
      meanConnectionsBasedOnParentFrames_G2[i,]<-c(0,0,0)
    }else{
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G2[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        meanConnectionsBasedOnParentFrames_G2[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G2[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        meanConnectionsBasedOnParentFrames_G2[i,]<-frechetMean(connectionsBasedOnParentFrames_G2[i,,]) 
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G2[i,,])
        meanConnectionsBasedOnParentFrames_G2[i,]<-pns(connectionsBasedOnParentFrames_G2[i,,],sphere.type = sphereType)$PNS$mean  
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
  }
  print("Done!")
}

# Calculate mean connection based on global coordinate (using mean frames in global coordinate)
if(TRUE){
  meanConnectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents2[i]
    meanConnectionsGlobalCoordinate_G1[k1,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,k2],
                                         vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G1[k1,])
  }
  meanConnectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents2[i]
    meanConnectionsGlobalCoordinate_G2[k1,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,k2],
                                         vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G2[k1,])
  }
}

#####################################################################################################
#####################################################################################################
# Convert mean LP-ds-rep to a GP-ds-rep

if(TRUE){
  meanPositions_G1<-array(NA,dim = c(numberOfFrames,3))
  meanPositions_G1[16,]<-c(0,0,0)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents2[i]
    meanPositions_G1[k1,]<-meanPositions_G1[k2,]+
      meanConnectionsLengths_G1[k1]*meanConnectionsGlobalCoordinate_G1[k1,]
    
  }
  meanPositions_G2<-array(NA,dim = c(numberOfFrames,3))
  meanPositions_G2[16,]<-c(0,0,0)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents2[i]
    meanPositions_G2[k1,]<-meanPositions_G2[k2,]+
      meanConnectionsLengths_G2[k1]*meanConnectionsGlobalCoordinate_G2[k1,]
    
  }
  meanSpokesTails_G1<-array(NA,dim = c(nTotalRadii,3))
  meanSpokesTips_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    frameOfSpokeNo<-NA
    if(i<=upSpoeksNumber){
      frameOfSpokeNo<-i
    }else if(i<=2*upSpoeksNumber){
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }
    
    meanSpokesTails_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]
    meanSpokesTips_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G1[i,]*radiiMean_G1[i]
  }
  meanSpokesTails_G2<-array(NA,dim = c(nTotalRadii,3))
  meanSpokesTips_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    frameOfSpokeNo<-NA
    if(i<=upSpoeksNumber){
      frameOfSpokeNo<-i
    }else if(i<=2*upSpoeksNumber){
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }
    
    meanSpokesTails_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]
    meanSpokesTips_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G2[i,]*radiiMean_G2[i]
  }
  print("Done!")
}
#####################################################################################################
#####################################################################################################
# Plot overlaid LP-ds-rep means of PD and CG

if(TRUE){
  open3d()
  srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  srep2<-rbind(meanSpokesTails_G2,meanSpokesTips_G2)*mean(sizes_G2)
  for (i in 1:nTotalRadii) {
    plot3d(srep2[c(i,(i+nTotalRadii)),],type="l",lwd = 1,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  # legend3d("topright", legend = paste(c('Mean G1', 'Mean G2')), pch = 16, col = c("blue","red"), cex=1, inset=c(0.02))
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  
}


#plot
if(TRUE){
  boundaryMeanTop_G1<-meanSpokesTips_G1[c(1:upSpoeksNumber,(upSpoeksNumber+downSpoeksNumber+1):nTotalRadii),]
  boundaryMeanBottom_G1<-meanSpokesTips_G1[c((upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),(upSpoeksNumber+downSpoeksNumber+1):nTotalRadii),]
  open3d()
  for (i in 2:numberOfFrames) {
    plot3d(rbind(meanPositions_G1[framesCenters[i],],meanPositions_G1[framesParents[i],])
           ,type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(boundaryMeanTop_G1[framesCenters[i],],boundaryMeanTop_G1[framesParents[i],])
           ,type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(rbind(boundaryMeanBottom_G1[framesCenters[i],],boundaryMeanBottom_G1[framesParents[i],])
           ,type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
    
  }
  plot3d(meanPositions_G1[sheetConnections,],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryMeanTop_G1[sheetConnections,],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  plot3d(boundaryMeanBottom_G1[sheetConnections,],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in framesCenters) {
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(meanPositions_G1[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  for (i in framesCenters) {
    vectors3d(boundaryMeanTop_G1[i,]+meanFramesTopGlobalCoordinate_G1[1,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryMeanTop_G1[i,]+meanFramesTopGlobalCoordinate_G1[2,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryMeanTop_G1[i,]+meanFramesTopGlobalCoordinate_G1[3,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  for (i in framesCenters) {
    vectors3d(boundaryMeanBottom_G1[i,]+meanFramesBottomGlobalCoordinate_G1[1,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryMeanBottom_G1[i,]+meanFramesBottomGlobalCoordinate_G1[2,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(boundaryMeanBottom_G1[i,]+meanFramesBottomGlobalCoordinate_G1[3,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  # for (i in 1:nTotalRadii) {
  #   plot3d(rbind(meanSpokesTails_G1[i,],
  #                meanSpokesTails_G1[i,]+(meanSpokesDirectionsGlobalCoordinate_G1[i,]*radiiMean_G1[i])),type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  # }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
}
#####################################################################################################
#####################################################################################################
# Hypothesis testing

# choose typeOfTest as "Parametric" or "Permutation"
# typeOfTest<-"Permutation" #nPerm default is 10000
typeOfTest<-"Parametric" # hotelling's test with normality assumption


# hypothesis test on LP size
pValues_LP_sizes<-meanDifferenceTest1D(LP_sizes_G1,LP_sizes_G2,type = typeOfTest) 
cat("pValue of LP sizes is:",pValues_LP_sizes,"\n")
boxplot(LP_sizes_G1, LP_sizes_G2, names = c("G1","G2"),main="LP size")
cat("sd LP size G1:",sd(LP_sizes_G1),"sd LP size G2:",sd(LP_sizes_G2),"\n")
cat("Mean LP size G1:",mean(LP_sizes_G1),"mean LP size G2:",mean(LP_sizes_G2),"\n")

# hypothesis test on spokes' lengths
pValues_TtestRadii<-rep(NA,nTotalRadii)
pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
for (i in 1:nTotalRadii) {
  setTxtProgressBar(pb, i) #create progress bar
  
  pValues_TtestRadii[i]<-meanDifferenceTest1D(log(radiiScaled_G1[i,]),
                                              log(radiiScaled_G2[i,]),
                                              type = typeOfTest)
}
# which(pValues_TtestRadii<=0.05)


# hypothesis test on connections' length
pValues_TtestConnectionsLengths<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for (i in 1:numberOfFrames) {
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValues_TtestConnectionsLengths[i]<-1
  }else{
    pValues_TtestConnectionsLengths[i]<-meanDifferenceTest1D(log(connectionsLengthsScaled_G1[i,]),
                                                             log(connectionsLengthsScaled_G2[i,]),
                                                             type = typeOfTest)
  }
}
# which(pValues_TtestConnectionsLengths<=0.05)


# choose euclideanization method typeOfStudy4directions as "PNS" or "tangent space"
# PNS takes a few minutes!
typeOfStudy4directions<-"tangent space"
# typeOfStudy4directions<-"PNS"

source("subFunctions/euclideanization.R")

# hypothesis test on spokes' directions based on local frames
pValspokesDirectionsBasedOnFrames<-rep(NA,nTotalRadii)
pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
for(i in 1:nTotalRadii){
  setTxtProgressBar(pb, i) #create progress bar
  #NB! euclideanization must contain two groups because it uses the pooled mean 
  euclideanizedTemp<-euclideanization(spokesDirectionsBasedOnFrames_G1[i,,],
                                      spokesDirectionsBasedOnFrames_G2[i,,],
                                      type = typeOfStudy4directions)
  
  pValspokesDirectionsBasedOnFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                       euclideanizedTemp$euclideanG2,
                                                                       type=typeOfTest)
  
}
# which(pValspokesDirectionsBasedOnFrames<=0.05)


# hypothesis test on connections' directions based on local frames
pValConnectionsBasedOnParentFrames<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValConnectionsBasedOnParentFrames[i]<-1
    next
  }
  
  euclideanizedTemp<-euclideanization(connectionsBasedOnParentFrames_G1[i,,],
                                      connectionsBasedOnParentFrames_G2[i,,],
                                      type = typeOfStudy4directions)
  
  pValConnectionsBasedOnParentFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                        euclideanizedTemp$euclideanG2,
                                                                        type=typeOfTest)

}
# which(pValConnectionsBasedOnParentFrames<=0.05)


# hypothesis test on frames
pValFrameZaxisBasedOnParent<-rep(NA,numberOfFrames)
pValFrameXaxisBasedOnParent<-rep(NA,numberOfFrames)
pValFrameYaxisBasedOnParent<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValFrameZaxisBasedOnParent[i]<-1
    pValFrameXaxisBasedOnParent[i]<-1
    pValFrameYaxisBasedOnParent[i]<-1
    next
  }
  euclideanizedTemp<-euclideanization(framesBasedOnParents_G1[1,,i,],
                                      framesBasedOnParents_G2[1,,i,],
                                      type = typeOfStudy4directions)
  
  pValFrameZaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type=typeOfTest)
  
  euclideanizedTemp<-euclideanization(framesBasedOnParents_G1[2,,i,],
                                      framesBasedOnParents_G2[2,,i,],
                                      type = typeOfStudy4directions)
  
  pValFrameXaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type=typeOfTest)
  
  euclideanizedTemp<-euclideanization(framesBasedOnParents_G1[3,,i,],
                                      framesBasedOnParents_G2[3,,i,],
                                      type = typeOfStudy4directions)
  
  pValFrameYaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type = typeOfTest)
}

# hypothesis test on frames Top
pValTopFrameZaxisBasedOnParent<-rep(NA,numberOfFrames)
pValTopFrameXaxisBasedOnParent<-rep(NA,numberOfFrames)
pValTopFrameYaxisBasedOnParent<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  euclideanizedTemp<-euclideanization(framesBasedOnParentsTop_G1[1,,i,],
                                      framesBasedOnParentsTop_G2[1,,i,],
                                      type = typeOfStudy4directions)
  
  pValTopFrameZaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type=typeOfTest)
  
  euclideanizedTemp<-euclideanization(framesBasedOnParentsTop_G1[2,,i,],
                                      framesBasedOnParentsTop_G2[2,,i,],
                                      type = typeOfStudy4directions)
  
  pValTopFrameXaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type=typeOfTest)
  
  euclideanizedTemp<-euclideanization(framesBasedOnParentsTop_G1[3,,i,],
                                      framesBasedOnParentsTop_G2[3,,i,],
                                      type = typeOfStudy4directions)
  
  pValTopFrameYaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type = typeOfTest)
}

# hypothesis test on frames Bottom
pValBottomFrameZaxisBasedOnParent<-rep(NA,numberOfFrames)
pValBottomFrameXaxisBasedOnParent<-rep(NA,numberOfFrames)
pValBottomFrameYaxisBasedOnParent<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  euclideanizedTemp<-euclideanization(framesBasedOnParentsBottom_G1[1,,i,],
                                      framesBasedOnParentsBottom_G2[1,,i,],
                                      type = typeOfStudy4directions)
  
  pValBottomFrameZaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type=typeOfTest)
  
  euclideanizedTemp<-euclideanization(framesBasedOnParentsBottom_G1[2,,i,],
                                      framesBasedOnParentsBottom_G2[2,,i,],
                                      type = typeOfStudy4directions)
  
  pValBottomFrameXaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type=typeOfTest)
  
  euclideanizedTemp<-euclideanization(framesBasedOnParentsBottom_G1[3,,i,],
                                      framesBasedOnParentsBottom_G2[3,,i,],
                                      type = typeOfStudy4directions)
  
  pValBottomFrameYaxisBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type = typeOfTest)
}

#####################################################################################################
#####################################################################################################
# plot significant GOPs

pValTopSpokesDirectionsBasedOnFrames<-pValspokesDirectionsBasedOnFrames[1:upSpoeksNumber] 
pValBottomSpokesDirectionsBasedOnFrames<-pValspokesDirectionsBasedOnFrames[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber)]   
pValCrestSpokesDirectionsBasedOnFrames<-pValspokesDirectionsBasedOnFrames[(upSpoeksNumber+downSpoeksNumber+1):nTotalRadii]   

pValTopSpokesPositions<-pValues_TtestRadii[1:upSpoeksNumber] 
pValBottomSpokesPositions<-pValues_TtestRadii[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber)]   
pValCrestSpokesPositions<-pValues_TtestRadii[(upSpoeksNumber+downSpoeksNumber+1):nTotalRadii]   


pvalues_LP_ds_rep <- c(pValTopSpokesPositions,             #n_f = 71-20
                       pValBottomSpokesPositions,          #n_f = 71-20
                       # pValCrestSpokesPositions,         
                       pValues_TtestConnectionsLengths[1:upSpoeksNumber],    #n_f = 71-20
                       pValTopSpokesDirectionsBasedOnFrames,   #n_f = 71-20
                       pValBottomSpokesDirectionsBasedOnFrames,  #n_f = 71-20  
                       # pValCrestSpokesDirectionsBasedOnFrames,   
                       pValConnectionsBasedOnParentFrames[1:upSpoeksNumber], #n_f = 71-20
                       pValFrameZaxisBasedOnParent[1:upSpoeksNumber],         #n_f = 71-20 
                       pValFrameXaxisBasedOnParent[1:upSpoeksNumber],         #n_f = 71-20
                       pValFrameYaxisBasedOnParent[1:upSpoeksNumber],         #n_f = 71-20
                       pValTopFrameZaxisBasedOnParent[1:upSpoeksNumber],      #n_f = 71-20
                       pValTopFrameXaxisBasedOnParent[1:upSpoeksNumber],      #n_f = 71-20
                       pValTopFrameYaxisBasedOnParent[1:upSpoeksNumber],     #n_f = 71-20
                       pValBottomFrameZaxisBasedOnParent[1:upSpoeksNumber],  #n_f = 71-20
                       pValBottomFrameXaxisBasedOnParent[1:upSpoeksNumber],  #n_f = 71-20
                       pValBottomFrameYaxisBasedOnParent[1:upSpoeksNumber])  #n_f = 71-20


# n_s<-nTotalRadii
# n_f<-numberOfFrames

n_f<-upSpoeksNumber
n_f

length(pvalues_LP_ds_rep)
alpha<-0.05
significantPvalues<-which(pvalues_LP_ds_rep<=alpha)
significantPvalues

#adjust p-values by Benjamini-Hochberg
FDR<-0.1
pvalues_LP_ds_rep_BH<-p.adjust(pvalues_LP_ds_rep,method = "BH")
significantPvalues_BH<-which(pvalues_LP_ds_rep_BH<=FDR)
significantPvalues_BH

#1
significantTopSpokesPositions<-significantPvalues[which(significantPvalues<=n_f)]
significantTopSpokesPositions
significantTopSpokesPositions_BH<-significantPvalues_BH[which(significantPvalues_BH<=n_f)]
significantTopSpokesPositions_BH
#2
significantBottomSpokesPositions<-significantPvalues[which(n_f+1<=significantPvalues
                                                          & significantPvalues<=(2*n_f))]-n_f
significantBottomSpokesPositions
significantBottomSpokesPositions_BH<-significantPvalues_BH[which(n_f+1<=significantPvalues_BH
                                                              & significantPvalues_BH<=(2*n_f))]-n_f
significantBottomSpokesPositions_BH
#3
significantConnectionsLengths<-significantPvalues[which(2*n_f+1<=significantPvalues
                                                        & significantPvalues<=(3*n_f))]-2*n_f
significantConnectionsLengths
significantConnectionsLengths_BH<-significantPvalues_BH[which((2*n_f+1)<=significantPvalues_BH
                                                              & significantPvalues_BH<=(3*n_f))]-2*n_f
significantConnectionsLengths_BH
#4
significantTopspokesDirections<-significantPvalues[which((3*n_f+1)<=significantPvalues
                                                      & significantPvalues<=(4*n_f))]-(3*n_f)
significantTopspokesDirections
significantTopspokesDirections_BH<-significantPvalues_BH[which((3*n_f+1)<=significantPvalues_BH
                                                            & significantPvalues_BH<=(4*n_f))]-(3*n_f)
significantTopspokesDirections_BH
#5
significantBottomspokesDirections<-significantPvalues[which((4*n_f+1)<=significantPvalues
                                                         & significantPvalues<=(5*n_f))]-(4*n_f)
significantBottomspokesDirections
significantBottomspokesDirections_BH<-significantPvalues_BH[which((4*n_f+1)<=significantPvalues_BH
                                                               & significantPvalues_BH<=(5*n_f))]-(4*n_f)
significantBottomspokesDirections_BH
#6
significantConnectionsDirections<-significantPvalues[which((5*n_f+1)<=significantPvalues
                                                           & significantPvalues<=(6*n_f))]-(5*n_f)
significantConnectionsDirections
significantConnectionsDirections_BH<-significantPvalues_BH[which((5*n_f+1)<=significantPvalues_BH
                                                              & significantPvalues_BH<=(6*n_f))]-(5*n_f)
significantConnectionsDirections_BH
#7
significantFrameZaxis<-significantPvalues[which((6*n_f+1)<=significantPvalues
                                                & significantPvalues<=(7*n_f))]-(6*n_f)
significantFrameZaxis
significantFrameZaxis_BH<-significantPvalues_BH[which((6*n_f+1)<=significantPvalues_BH
                                                   & significantPvalues_BH<=(7*n_f))]-(6*n_f)
significantFrameZaxis_BH
#8
significantFrameXaxis<-significantPvalues[which((7*n_f+1)<=significantPvalues
                                                & significantPvalues<=(8*n_f))]-(7*n_f)
significantFrameXaxis
significantFrameXaxis_BH<-significantPvalues_BH[which((7*n_f+1)<=significantPvalues_BH
                                                   & significantPvalues_BH<=(8*n_f))]-(7*n_f)
significantFrameXaxis_BH
#9
significantFrameYaxis<-significantPvalues[which((8*n_f+1)<=significantPvalues
                                                & significantPvalues<=(9*n_f))]-(8*n_f)
significantFrameYaxis
significantFrameYaxis_BH<-significantPvalues_BH[which((8*n_f+1)<=significantPvalues_BH
                                                   & significantPvalues_BH<=(9*n_f))]-(8*n_f)
significantFrameYaxis_BH
#10
significantTopFrameZaxis<-significantPvalues[which((9*n_f+1)<=significantPvalues
                                                   & significantPvalues<=(10*n_f))]-(9*n_f)
significantTopFrameZaxis
significantTopFrameZaxis_BH<-significantPvalues_BH[which((9*n_f+1)<=significantPvalues_BH
                                                      & significantPvalues_BH<=(10*n_f))]-(9*n_f)
significantTopFrameZaxis_BH
#11
significantTopFrameXaxis<-significantPvalues[which((10*n_f+1)<=significantPvalues
                                                   & significantPvalues<=(11*n_f))]-(10*n_f)
significantTopFrameXaxis
significantTopFrameXaxis_BH<-significantPvalues_BH[which((10*n_f+1)<=significantPvalues_BH
                                                      & significantPvalues_BH<=(11*n_f))]-(10*n_f)
significantTopFrameXaxis_BH
#12
significantTopFrameYaxis<-significantPvalues[which((11*n_f+1)<=significantPvalues
                                                   & significantPvalues<=(12*n_f))]-(11*n_f)
significantTopFrameYaxis
significantTopFrameYaxis_BH<-significantPvalues_BH[which((11*n_f+1)<=significantPvalues_BH
                                                      & significantPvalues_BH<=(12*n_f))]-(11*n_f)
significantTopFrameYaxis_BH
#13
significantBottomFrameZaxis<-significantPvalues[which((12*n_f+1)<=significantPvalues
                                                      & significantPvalues<=(13*n_f))]-(12*n_f)
significantBottomFrameZaxis
significantBottomFrameZaxis_BH<-significantPvalues_BH[which((12*n_f+1)<=significantPvalues_BH
                                                         & significantPvalues_BH<=(13*n_f))]-(12*n_f)
significantBottomFrameZaxis_BH
#14
significantBottomFrameXaxis<-significantPvalues[which((13*n_f+1)<=significantPvalues
                                                      & significantPvalues<=(14*n_f))]-(13*n_f)
significantBottomFrameXaxis
significantBottomFrameXaxis_BH<-significantPvalues_BH[which((13*n_f+1)<=significantPvalues_BH
                                                         & significantPvalues_BH<=(14*n_f))]-(13*n_f)
significantBottomFrameXaxis_BH
#15
significantBottomFrameYaxis<-significantPvalues[which((14*n_f+1)<=significantPvalues
                                                      & significantPvalues<=(15*n_f))]-(14*n_f)
significantBottomFrameYaxis
significantBottomFrameYaxis_BH<-significantPvalues_BH[which((14*n_f+1)<=significantPvalues_BH
                                                         & significantPvalues_BH<=(15*n_f))]-(14*n_f)
significantBottomFrameYaxis_BH


#######################
# final plot
open3d()
boundaryMeanTop_G1<-meanSpokesTips_G1[c(1:upSpoeksNumber,(upSpoeksNumber+downSpoeksNumber+1):nTotalRadii),]
boundaryMeanBottom_G1<-meanSpokesTips_G1[c((upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),(upSpoeksNumber+downSpoeksNumber+1):nTotalRadii),]
criticalPoints<-c(1,2,3,31,32,33,52:71)

#plot significants
#1
for (i in significantTopSpokesPositions[!significantTopSpokesPositions %in% criticalPoints]) {
  # c3d <- cube3d(color="magenta")
  # c3d2 <- c3d %>% scale3d(0.2,0.2,0.2) %>%
  #   translate3d(x = boundaryMeanTop_G1[i,1],y = boundaryMeanTop_G1[i,2],z = boundaryMeanTop_G1[i,3])
  # shade3d(c3d2)
  spheres3d(boundaryMeanTop_G1[i,], radius = 0.32,col = "red")
}
#2
for (i in significantBottomSpokesPositions[!significantBottomSpokesPositions %in% criticalPoints]) {
  # c3d <- cube3d(color="magenta")
  # c3d2 <- c3d %>% scale3d(0.2,0.2,0.2) %>%
  #   translate3d(x = boundaryMeanBottom_G1[i,1],y = boundaryMeanBottom_G1[i,2],z = boundaryMeanBottom_G1[i,3])
  # shade3d(c3d2)
  spheres3d(boundaryMeanBottom_G1[i,], radius = 0.32,col = "red")
}
#3
for (i in significantConnectionsLengths[!significantConnectionsLengths %in% criticalPoints]) {
  # c3d <- cube3d(color="magenta")
  # c3d2 <- c3d %>% scale3d(0.2,0.2,0.2) %>%
  #   translate3d(x = meanPositions_G1[i,1],y = meanPositions_G1[i,2],z = meanPositions_G1[i,3])
  # shade3d(c3d2)
  spheres3d(meanPositions_G1[i,], radius = 0.32,col = "red")
}
#4
for (i in significantTopspokesDirections[!significantTopspokesDirections %in% criticalPoints]) {
  spheres3d(boundaryMeanTop_G1[i,], radius = 0.32,col = "red")
}
#5
for (i in significantBottomspokesDirections[!significantBottomspokesDirections %in% criticalPoints]) {
  spheres3d(boundaryMeanBottom_G1[i,], radius = 0.32,col = "red")
}
#6
for (i in significantConnectionsDirections[!significantConnectionsDirections %in% criticalPoints]) {
  spheres3d(meanPositions_G1[i,], radius = 0.32,col = "red")
}
# significant frames
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  if((i %in% significantFrameZaxis & i %in% significantFrameXaxis) | 
     (i %in% significantFrameZaxis & i %in% significantFrameYaxis) |
     (i %in% significantFrameXaxis & i %in% significantFrameYaxis)){
    vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
  }
}
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  if((i %in% significantTopFrameZaxis & i %in% significantTopFrameXaxis) | 
     (i %in% significantTopFrameZaxis & i %in% significantTopFrameYaxis) |
     (i %in% significantTopFrameXaxis & i %in% significantTopFrameYaxis)){
    vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[1,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[2,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[3,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
  }
}
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  if((i %in% significantBottomFrameZaxis & i %in% significantBottomFrameXaxis) | 
     (i %in% significantBottomFrameZaxis & i %in% significantBottomFrameYaxis) |
     (i %in% significantBottomFrameXaxis & i %in% significantBottomFrameYaxis)){
    vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[1,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[2,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[3,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
  }
}
# #7
# open3d()
# for (i in significantFrameZaxis[!significantFrameZaxis %in% criticalPoints]) {
#   vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #8
# for (i in significantFrameXaxis[!significantFrameXaxis %in% criticalPoints]) {
#   vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #9
# for (i in significantFrameYaxis[!significantFrameYaxis %in% criticalPoints]) {
#   vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #10
# for (i in significantTopFrameZaxis[!significantTopFrameZaxis %in% criticalPoints]) {
#   vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[1,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #11
# for (i in significantTopFrameXaxis[!significantTopFrameXaxis %in% criticalPoints]) {
#   vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[2,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #12
# for (i in significantTopFrameYaxis[!significantTopFrameYaxis %in% criticalPoints]) {
#   vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[3,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #13
# for (i in significantBottomFrameZaxis[!significantBottomFrameZaxis %in% criticalPoints]) {
#   vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[1,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #14
# for (i in significantBottomFrameXaxis[!significantBottomFrameXaxis %in% criticalPoints]) {
#   vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[2,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #15
# for (i in significantBottomFrameYaxis[!significantBottomFrameYaxis %in% criticalPoints]) {
#   vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[3,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }

#plot all the rest
#plot the center point
plot3d(rbind(meanPositions_G1[16,],meanPositions_G1[16,])
       ,type="s",radius = 0.45,col = "darkblue",expand = 10,box=FALSE,add = TRUE) 
#plot spine
spineIndex <-c(33,32,31,28,25,22,19,16,13,10,7,4,1,2,3)
# spineIndex <-c(28,25,22,19,16,13,10,7,4)
for (i in spineIndex) {
  parent_Index<-which(framesCenters==i)
  plot3d(rbind(meanPositions_G1[i,],meanPositions_G1[framesParents[parent_Index],])
         ,type="l",lwd = 4,col = "black",expand = 10,box=FALSE,add = TRUE) 
}
for (i in 2:numberOfFrames) {
  plot3d(rbind(meanPositions_G1[framesCenters[i],],meanPositions_G1[framesParents[i],])
         ,type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
  plot3d(rbind(boundaryMeanTop_G1[framesCenters[i],],boundaryMeanTop_G1[framesParents[i],])
         ,type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
  plot3d(rbind(boundaryMeanBottom_G1[framesCenters[i],],boundaryMeanBottom_G1[framesParents[i],])
         ,type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
  
}
plot3d(meanPositions_G1[sheetConnections,],type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
plot3d(boundaryMeanTop_G1[sheetConnections,],type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
plot3d(boundaryMeanBottom_G1[sheetConnections,],type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)

for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[1,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[2,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[3,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[1,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[2,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[3,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
}
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  spheres3d(meanPositions_G1[i,], radius = 0.32,col = "lightblue")
  spheres3d(boundaryMeanTop_G1[i,], radius = 0.32,col = "lightblue")
  spheres3d(boundaryMeanBottom_G1[i,], radius = 0.32,col = "lightblue")
}
# decorate3d(xlab = "x", ylab = "y", zlab = "z",
#            box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
# pp <- par3d(no.readonly=TRUE)
# dput(pp, file="plotView.R", control = "all")
# pp <- dget("plotView.R")
par3d(userMatrix=pp$userMatrix)

################
# BH adjustment
#plot significants
open3d()
#1
for (i in significantTopSpokesPositions_BH[!significantTopSpokesPositions_BH %in% criticalPoints]) {
  # c3d <- cube3d(color="magenta")
  # c3d2 <- c3d %>% scale3d(0.2,0.2,0.2) %>%
  #   translate3d(x = boundaryMeanTop_G1[i,1],y = boundaryMeanTop_G1[i,2],z = boundaryMeanTop_G1[i,3])
  # shade3d(c3d2)
  spheres3d(boundaryMeanTop_G1[i,], radius = 0.32,col = "red")
}
#2
for (i in significantBottomSpokesPositions_BH[!significantBottomSpokesPositions_BH %in% criticalPoints]) {
  # c3d <- cube3d(color="magenta")
  # c3d2 <- c3d %>% scale3d(0.2,0.2,0.2) %>%
  #   translate3d(x = boundaryMeanBottom_G1[i,1],y = boundaryMeanBottom_G1[i,2],z = boundaryMeanBottom_G1[i,3])
  # shade3d(c3d2)
  spheres3d(boundaryMeanBottom_G1[i,], radius = 0.32,col = "red")
}
#3
for (i in significantConnectionsLengths_BH[!significantConnectionsLengths_BH %in% criticalPoints]) {
  # c3d <- cube3d(color="magenta")
  # c3d2 <- c3d %>% scale3d(0.2,0.2,0.2) %>%
  #   translate3d(x = meanPositions_G1[i,1],y = meanPositions_G1[i,2],z = meanPositions_G1[i,3])
  # shade3d(c3d2)
  spheres3d(meanPositions_G1[i,], radius = 0.32,col = "red")
}
#4
for (i in significantTopspokesDirections_BH[!significantTopspokesDirections_BH %in% criticalPoints]) {
  spheres3d(boundaryMeanTop_G1[i,], radius = 0.32,col = "red")
}
#5
for (i in significantBottomspokesDirections_BH[!significantBottomspokesDirections_BH %in% criticalPoints]) {
  spheres3d(boundaryMeanBottom_G1[i,], radius = 0.32,col = "red")
}
#6
for (i in significantConnectionsDirections_BH[!significantConnectionsDirections_BH %in% criticalPoints]) {
  spheres3d(meanPositions_G1[i,], radius = 0.32,col = "red")
}
# significant frames
k<-0
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  if((i %in% significantFrameZaxis_BH & i %in% significantFrameXaxis_BH) | 
     (i %in% significantFrameZaxis_BH & i %in% significantFrameYaxis_BH) |
     (i %in% significantFrameXaxis_BH & i %in% significantFrameYaxis_BH)){
    vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    k<-k+1
  }
}
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  if((i %in% significantTopFrameZaxis_BH & i %in% significantTopFrameXaxis_BH) | 
     (i %in% significantTopFrameZaxis_BH & i %in% significantTopFrameYaxis_BH) |
     (i %in% significantTopFrameXaxis_BH & i %in% significantTopFrameYaxis_BH)){
    vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[1,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[2,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[3,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    k<-k+1
  }
}
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  if((i %in% significantBottomFrameZaxis_BH & i %in% significantBottomFrameXaxis_BH) | 
     (i %in% significantBottomFrameZaxis_BH & i %in% significantBottomFrameYaxis_BH) |
     (i %in% significantBottomFrameXaxis_BH & i %in% significantBottomFrameYaxis_BH)){
    vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[1,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[2,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[3,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
    k<-k+1
  }
}
print(k)

# #7
# for (i in significantFrameZaxis_BH[!significantFrameZaxis_BH %in% criticalPoints]) {
#   vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #8
# for (i in significantFrameXaxis_BH[!significantFrameXaxis_BH %in% criticalPoints]) {
#   vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #9
# for (i in significantFrameYaxis_BH[!significantFrameYaxis_BH %in% criticalPoints]) {
#   vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #10
# for (i in significantTopFrameZaxis_BH[!significantTopFrameZaxis_BH %in% criticalPoints]) {
#   vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[1,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #11
# for (i in significantTopFrameXaxis_BH[!significantTopFrameXaxis_BH %in% criticalPoints]) {
#   vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[2,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #12
# for (i in significantTopFrameYaxis_BH[!significantTopFrameYaxis_BH %in% criticalPoints]) {
#   vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[3,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #13
# for (i in significantBottomFrameZaxis_BH[!significantBottomFrameZaxis_BH %in% criticalPoints]) {
#   vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[1,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #14
# for (i in significantBottomFrameXaxis_BH[!significantBottomFrameXaxis_BH %in% criticalPoints]) {
#   vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[2,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }
# #15
# for (i in significantBottomFrameYaxis_BH[!significantBottomFrameYaxis_BH %in% criticalPoints]) {
#   vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[3,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="red", lwd=4)
# }

#plot all the rest
#plot the center point
plot3d(rbind(meanPositions_G1[16,],meanPositions_G1[16,])
       ,type="s",radius = 0.45,col = "darkblue",expand = 10,box=FALSE,add = TRUE) 
#plot spine
spineIndex <-c(33,32,31,28,25,22,19,16,13,10,7,4,1,2,3)
# spineIndex <-c(28,25,22,19,16,13,10,7,4)
for (i in spineIndex) {
  parent_Index<-which(framesCenters==i)
  plot3d(rbind(meanPositions_G1[i,],meanPositions_G1[framesParents[parent_Index],])
         ,type="l",lwd = 4,col = "black",expand = 10,box=FALSE,add = TRUE) 
}
for (i in 2:numberOfFrames) {
  plot3d(rbind(meanPositions_G1[framesCenters[i],],meanPositions_G1[framesParents[i],])
         ,type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
  plot3d(rbind(boundaryMeanTop_G1[framesCenters[i],],boundaryMeanTop_G1[framesParents[i],])
         ,type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
  plot3d(rbind(boundaryMeanBottom_G1[framesCenters[i],],boundaryMeanBottom_G1[framesParents[i],])
         ,type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
  
}
plot3d(meanPositions_G1[sheetConnections,],type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
plot3d(boundaryMeanTop_G1[sheetConnections,],type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)
plot3d(boundaryMeanBottom_G1[sheetConnections,],type="l",lwd = 0.5,col = "grey",expand = 10,box=FALSE,add = TRUE)

for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[1,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[2,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(meanPositions_G1[i,]+1.5*meanFramesGlobalCoordinate_G1[3,,i],origin = meanPositions_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[1,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[2,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanTop_G1[i,]+1.5*meanFramesTopGlobalCoordinate_G1[3,,i],origin = boundaryMeanTop_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[1,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[2,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
  vectors3d(boundaryMeanBottom_G1[i,]+1.5*meanFramesBottomGlobalCoordinate_G1[3,,i],origin = boundaryMeanBottom_G1[i,],headlength = 0.2,radius = 2/10, col="lightblue", lwd=4)
}
for (i in framesCenters[!framesCenters %in% criticalPoints]) {
  spheres3d(meanPositions_G1[i,], radius = 0.32,col = "lightblue")
  spheres3d(boundaryMeanTop_G1[i,], radius = 0.32,col = "lightblue")
  spheres3d(boundaryMeanBottom_G1[i,], radius = 0.32,col = "lightblue")
}
# decorate3d(xlab = "x", ylab = "y", zlab = "z",
#            box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
pp <- par3d(no.readonly=TRUE)
dput(pp, file="plotView.R", control = "all")
pp <- dget("plotView.R")
par3d(userMatrix=pp$userMatrix)
