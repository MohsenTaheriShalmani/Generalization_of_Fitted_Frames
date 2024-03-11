
secondFrameVectorFunction <- function(p1,vertex,p2,normalVec) {
  
  u1<-convertVec2unitVec(p2-vertex)
  u2<-convertVec2unitVec(p1-vertex)
  
  v<-convertVec2unitVec(u1-u2)

  #projection of tempPoint on tangent space
  #sum(v*normalVec) is dot product that is the point distance from the plane
  projected_point <- v-(sum(v*normalVec)/norm(normalVec,type = "2"))*normalVec
  
  secondVec<-convertVec2unitVec(projected_point)
  
  return(secondVec)
}

# generate the frames
frameGenerator <- function(centeredSkel,medialNormals,
                           framesCenters,framesBackPoints,framesFronts) {
  
  numberOfFrames<-length(framesCenters)
  
  # NB!!! number of frames is equal to upSpoeksNumber
  framesFirstVec<-array(NA,dim = c(numberOfFrames,3))
  framesSecondVec<-array(NA,dim = c(numberOfFrames,3))
  framesThirdVec<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    
    a<-framesBackPoints[i]
    b<-framesCenters[i]
    c<-framesFronts[i]
    
    framesFirstVec[b,]<-medialNormals[b,]
    
    if(c==Inf){
      framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[a,],vertex = centeredSkel[b,],
                                                     p2 =convertVec2unitVec(centeredSkel[b,]-centeredSkel[a,])+centeredSkel[b,],
                                                     normalVec = medialNormals[b,])
    }else{
      framesSecondVec[b,]<-secondFrameVectorFunction(p1 = centeredSkel[a,],vertex = centeredSkel[b,],
                                                             p2 = centeredSkel[c,],normalVec = medialNormals[b,])
    }
    
    framesThirdVec[b,]<-convertVec2unitVec(myCrossProduct(framesFirstVec[b,],framesSecondVec[b,]))
  }
  
  result<-list(framesFirstVec=framesFirstVec,
               framesSecondVec=framesSecondVec,
               framesThirdVec=framesThirdVec)
}


