source("subFunctions/ray_triangle_intersection.R")

cutAndStretchSpokes <- function(BoundaryPDM,SkeletalPDM,spharmPDM,meshPolygonData) {
  numberOfPoints<-dim(spharmPDM)[1]
  tipOfCuttedSpokes<-array(NA, dim = dim(BoundaryPDM))
  for (i in 1:nTotalRadii) {
    spokeDirectionGlobal<-convertVec2unitVec(BoundaryPDM[i,]-SkeletalPDM[i,])
    spokeTail<-SkeletalPDM[i,]
    intersections<-array(NA,dim = c(dim(meshPolygonData)[1],3))
    for (j in 1:dim(meshPolygonData)[1]) {
      p1<-polyMatrix[j,1]
      p2<-polyMatrix[j,2]
      p3<-polyMatrix[j,3]
      point1<-spharmPDM[p1,]
      point2<-spharmPDM[p2,]
      point3<-spharmPDM[p3,]
      tempIntersection<-rayTriangleIntersection(rayOrigin = spokeTail,rayDirection = spokeDirectionGlobal,
                                                triangleVertex1 = point1,triangleVertex2 = point2,triangleVertex3 = point3)
      intersections[j,]<-tempIntersection
    }
    distances<-rep(Inf,dim(meshPolygonData)[1])
    for (k in 1:dim(meshPolygonData)[1]) {
      if(!is.na(intersections[k,1])){
        distances[k]<-norm(intersections[k,]-spokeTail,type = "2")
      }
    }
    
    tipOfCuttedSpokes[i,]<-intersections[which.min(distances),]
  }
  return(tipOfCuttedSpokes) #new boundary
}
