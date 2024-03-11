library(shapes)
library(RiemBase)
source("subFunctions/MathFunctions.R")

euclideanization <- function(directions1, directions2, type="tangent space") {
  
  nSamplesG1<-dim(directions1)[2]
  nSamplesG2<-dim(directions2)[2]
  
  pooledDirection<-cbind(directions1,directions2)
  
  # For extremely concentrated data we use Mardia mean direction
  pcaData1<-prcomp(t(pooledDirection))
  if(pcaData1$sdev[1]<1e-02 | pcaData1$sdev[2]<1e-02){
    mu_g<-convertVec2unitVec(colMeans(t(pooledDirection)))
    
    R <- rotMat(mu_g,c(0,0,1))
    shiftedG1<-R%*%directions1
    shiftedG2<-R%*%directions2
    
    #log transfer to the tangent space
    logG1<-t(LogNPd(shiftedG1))
    logG2<-t(LogNPd(shiftedG2))
    
    result<-list(euclideanG1=logG1, euclideanG2=logG2)
    
  }else if(type=="tangent space"){
    
    # mean by Fre'chet mean
    allDirTemp<-t(pooledDirection)
    data1 <- list()
    for (j in 1:dim(allDirTemp)[1]){
      data1[[j]] <-allDirTemp[j,]
    }
    data2 <- riemfactory(data1, name="sphere")
    # Fre'chet Mean
    out1<- rbase.mean(data2)
    mu_g<-as.vector(out1$x)
    
    #rotate data to the north pole 
    R <- rotMat(mu_g,c(0,0,1))
    shiftedG1<-R%*%directions1
    shiftedG2<-R%*%directions2
    
    #log transfer to the tangent space
    logG1<-t(LogNPd(shiftedG1))
    logG2<-t(LogNPd(shiftedG2))
    
    # theta1<-acos(shiftedG1[3,])
    # logPointG1x<-shiftedG1[1,]*theta1/sin(theta1)
    # logPointG1y<-shiftedG1[2,]*theta1/sin(theta1)
    # logPointG1z<-rep(1,nSamplesG1)
    # logG1<-cbind(logPointG1x,logPointG1y,logPointG1z)
    # theta2<-acos(shiftedG2[3,])
    # logPointG2x<-shiftedG2[1,]*theta2/sin(theta2)
    # logPointG2y<-shiftedG2[2,]*theta2/sin(theta2)
    # logPointG2z<-rep(1,nSamplesG2)
    # logG2<-cbind(logPointG2x,logPointG2y,logPointG2z)
    # result<-list(euclideanG1=logG1[,1:2], euclideanG2=logG2[,1:2])
    
    result<-list(euclideanG1=logG1, euclideanG2=logG2)
    
  }else if(type=="PNS"){
    
    # use pns instead of Fre'chet mean
    typeOfSphere<-kurtosisTestFunction(pooledDirection)
    pnsDirection<-pns(pooledDirection,sphere.type = typeOfSphere) #pooled directions
    res_G1<-t(pnsDirection$resmat[,1:nSamplesG1])
    res_G2<-t(pnsDirection$resmat[,(nSamplesG1+1):(nSamplesG1+nSamplesG2)])
    
    result<-list(euclideanG1=res_G1, euclideanG2=res_G2)
    
  }else{
    stop("Please choose the type of analysis i.e., 'PNS' or 'tangent space'!")
  }
  
  return(result)
}

