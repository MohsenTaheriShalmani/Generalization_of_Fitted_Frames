library(RiemBase)
library(pracma)
library(matlib)

# Math functions
# function to calculate Mahalanobis distance (Hotelling Metric)
MahalanobisDistance<-function(X,Y){
  k<-dim(X)[2]
  nx<-dim(X)[1]
  ny<-dim(Y)[1]
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  # T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  T2<-t(meanX-meanY)%*%Ginv(S*(1/nx+1/ny))%*%(meanX-meanY) # Ginv is the Moore-Penrose generalized inverse
  return(as.double(T2))
}

MahalanobisDistance_1D<-function(X,Y){
  nx<-length(X)
  ny<-length(Y)
  Sx<-var(X)
  Sy<-var(Y)
  meanX<-mean(X)
  meanY<-mean(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled variance
  T2<-abs(meanX-meanY)*(1/(S*sqrt(1/nx+1/ny)))
  return(T2)
}

# Hotelling T2 test
HotellingT2<-function(X,Y){
  if(dim(X)[2]!=dim(Y)[2]){
    stop("Dimention Error!\n")
  }
  k<-dim(X)[2]
  nx<-dim(X)[1]
  ny<-dim(Y)[1]
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  # T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  T2<-t(meanX-meanY)%*%Ginv(S*(1/nx+1/ny))%*%(meanX-meanY) # Ginv is the Moore-Penrose generalized inverse
  F_value<-((n-k)/(k*(n-1)))*T2
  df1<-k
  df2<-n-k
  p_value<-1-pf(F_value,df1,df2)
  return(p_value)
}

permutation1D <- function(A,B,nPerm=10000) {
  if(!is.null(dim(A)) | !is.null(dim(B))){
    stop("Dimention of A & B must be equal 1!")
  }
  T0<-MahalanobisDistance_1D(A,B) #observed test statistics
  nSamplesG1<-length(A)
  nSamplesG2<-length(B)
  testStatistics<-rep(NA,nPerm)
  pooledGroup<-c(A,B)
  for (j in 1:nPerm) {
    permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
    groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
    groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
    
    tempX<-pooledGroup[groupA_elementNumbers]
    tempY<-pooledGroup[groupB_elementNumbers]
    
    testStatistics[j]<-MahalanobisDistance_1D(tempX,tempY)
  }
  pVal<-sum(testStatistics>T0)/nPerm
  return(pVal)
}

permutationMultivariate <- function(A,B,nPerm=10000) {
  if(dim(A)[2]!=dim(B)[2] | is.null(dim(A)) | is.null(dim(B))){
    stop("Dimention of A & B must be equal and greater than 1!")
  }
  T0<-MahalanobisDistance(A,B) #observed test statistics
  nSamplesG1<-dim(A)[1]
  nSamplesG2<-dim(B)[1]
  testStatistics<-rep(NA,nPerm)
  pooledGroup<-rbind(A,B)
  for (j in 1:nPerm) {
    permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
    groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
    groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
    
    tempX<-pooledGroup[groupA_elementNumbers,]
    tempY<-pooledGroup[groupB_elementNumbers,]
    
    testStatistics[j]<-MahalanobisDistance(tempX,tempY)
  }
  pVal<-sum(testStatistics>T0)/nPerm
  return(pVal)
}


meanDifferenceTest1D <- function(A,B,type="Parametric",nPerm=10000) {
  if(!is.null(dim(A)) | !is.null(dim(B))){
    stop("Dimention of A & B must be equal 1!")
  }
  if(type=="Parametric"){
    pVal<-t.test(A,B)$p.value
    return(pVal)
  }else if(type=="Permutation"){
    
    T0<-MahalanobisDistance_1D(A,B) #observed test statistics
    nSamplesG1<-length(A)
    nSamplesG2<-length(B)
    testStatistics<-rep(NA,nPerm)
    pooledGroup<-c(A,B)
    for (j in 1:nPerm) {

      permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
      groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
      groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
      
      tempX<-pooledGroup[groupA_elementNumbers]
      tempY<-pooledGroup[groupB_elementNumbers]
      
      testStatistics[j]<-MahalanobisDistance_1D(tempX,tempY)
    }
    pVal<-(0.5+sum(testStatistics>=T0))/(nPerm+1)
    return(pVal)
  }else{
    stop("Please verify the type as Parametric or Permutation!")
  }
}

meanDifferenceTestMultivariate <- function(A,B,type="Parametric",nPerm=10000) {
  if(dim(A)[2]!=dim(B)[2] | is.null(dim(A)) | is.null(dim(B))){
    stop("Dimention of A & B must be equal and greater than 1!")
  }
  if(type=="Parametric"){
    pVal<-HotellingT2(A,B)
    return(pVal)
  }else if(type=="Permutation"){
    
    T0<-MahalanobisDistance(A,B) #observed test statistics
    nSamplesG1<-dim(A)[1]
    nSamplesG2<-dim(B)[1]
    testStatistics<-rep(NA,nPerm)
    pooledGroup<-rbind(A,B)
    for (j in 1:nPerm) {
      
      permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
      groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
      groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
      
      tempX<-pooledGroup[groupA_elementNumbers,]
      tempY<-pooledGroup[groupB_elementNumbers,]
      
      testStatistics[j]<-MahalanobisDistance(tempX,tempY)
    }
    pVal<-(0.5+sum(testStatistics>=T0))/(nPerm+1)
    return(pVal)
  }else{
    stop("Please verify the type as Parametric or Permutation!")
  }
}



# library(BisRNA)
# The fisher.method function is also available in BisRNA library
fisher.method<-function (pvalues)
{
  df <- 2 * length(pvalues)
  global_pValue<-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = FALSE)
  # global_pValue<-1-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = TRUE)
  return(global_pValue)
}

# convert vectors to unit vectors
convertVec2unitVec <- function(vec) {
  if(norm(vec,type = "2")==0){
    stop("vector is zero!")
  }
  return(vec/norm(vec,type = "2"))
}

# cross product of 2 vectors
myCrossProduct <- function(v,u) {
  return(c(v[2]*u[3]-v[3]*u[2],v[3]*u[1]-v[1]*u[3],v[1]*u[2]-v[2]*u[1]))
}

# frechet mean

frechetMean <- function(directions) {
  
  allDirTemp<-t(directions)
  data1 <- list()
  for (j in 1:dim(allDirTemp)[1]){
    data1[[j]] <-allDirTemp[j,]
  }
  data2 <- riemfactory(data1, name="sphere")
  ### Compute Fre'chet Mean
  out1<- rbase.mean(data2)
  meanFrechet<-as.vector(out1$x)
  
  return(meanFrechet)
  
}


# calculate unit normal vector of triangle mesh
unitNormalOfTriangle <- function(point1,point2,point3) {
  a<-point2-point1
  b<-point3-point1
  
  normalVec<-c((a[2]*b[3]-a[3]*b[2]),-(a[1]*b[3]-a[3]*b[1]),(a[1]*b[2]-a[2]*b[1]))
  triangleArea<-sqrt(sum(normalVec^2))/2
  unitNormal<-normalVec/sqrt(sum(normalVec^2))
  
  return(unitNormal)
  
}

#this function is exactly like rgl::rotate3d(vec, angle, x, y, z)
rotationAboutXYZaxis <- function(vector,angle,axis=1) {
  if(axis==1){
    rotationMatrixTemp<-matrix(c(1,0,0,0,cos(angle),-sin(angle),0,sin(angle),cos(angle)),nrow = 3,byrow =T)
  }else if(axis==2){
    rotationMatrixTemp<-matrix(c(cos(angle),0,sin(angle),0,1,0,-sin(angle),0,cos(angle)),nrow = 3,byrow =T)  
  }else if(axis==3){
    rotationMatrixTemp<-matrix(c(cos(angle),-sin(angle),0,sin(angle),cos(angle),0,0,0,1),nrow = 3,byrow =T)  
  }else{
    stop("Axis is not acceptable! It must be 1 for X 2 for Y and 3 for Z")
  }
  return(as.vector(vector%*%rotationMatrixTemp))
}

# rotate3d(c(2,2,2), pi/4, 1, 0, 0)
# rotationAboutXYZaxis(c(2,2,2),angle = pi/4,axis = 1)


# generate random von Mises distribution on circle in radian
# converted code of Sungkyu Jung 2013, and Byungwon Kim 2017, MATLAB randvonMises.m
# mean is in[0,2pi] and kappa>0
randVonMises <- function(mean, kappa, n) {
  tau<-1+sqrt(1+4*kappa^2)
  rho<-(tau-sqrt(2*tau))/(2*kappa)
  r<-(1+rho^2)/(2*rho)
  
  u1<-runif(n)
  z<-cos(pi*u1)
  f<-(1+r*z)/(r+z)
  c<-kappa*(r-f)
  u2<-runif(n)
  acceptanceid<-(c*(2-c)-u2>0) | (log(c/u2)+1-c>=0)
  u3<-runif(sum(acceptanceid))
  theta<-sign(u3-0.5)*acos(f[acceptanceid])
  nnow<-length(theta)
  
  while (n>nnow) {
    n_more<-ceiling(n/nnow*(n-nnow))
    u1<-runif(n_more)
    z<-cos(pi*u1)
    f<-(1+r*z)/(r+z)
    c<-kappa*(r-f)
    u2<-runif(n_more)
    acceptanceid<-(c*(2-c)-u2>0) | (log(c/u2)+1-c>=0)
    u3<-runif(sum(acceptanceid))
    thetamore<-sign(u3-0.5)*acos(f[acceptanceid])
    
    theta<-c(theta, thetamore)
    nnow<-length(theta)
  }
  
  theta<-theta[1:n] + mean
  
  return(theta)
}

# randVonMisesSamples<-randVonMises(mean = pi/4,kappa = 10,n = 2000)
# hist(randVonMisesSamples)
# plotshapes(cbind(cos(randVonMisesSamples),sin(randVonMisesSamples)))


# converted code of Sungkyu Jung MATLAB randS2.m
# generate random sample of small sphere on S2 of the second kind 
# mu0, mu1 are directions and kappa0>1, kappa1>1
randS2 <- function(mu0,mu1,kappa0,kappa1,n) {
  
  mu0<-mu0/norm(mu0,type = "2")
  mu1<-mu1/norm(mu1,type = "2")
  nu<-sum(mu0*mu1)
  
  #generate Bingham-Mardia random vectors by the north pole
  
  x<-rnorm(n = n,mean = nu,sd = 1/sqrt(2*kappa0))
  x<-x[x<1 & x>-1]
  nnow<-length(x)
  while(n>nnow){
    n_more<-ceiling(n/nnow*(n-nnow))
    xmore<-rnorm(n = n_more,mean = nu,sd = 1/sqrt(2*kappa0))
    xmore<-xmore[xmore<1 & xmore>-1]
    x<-c(x,xmore)
    nnow<-length(x)
  }
  z<-x[1:n]
  
  #generate von Mises for longitude that c=mu0-nu*mu1 is parallel to x-axis
  theta<-randVonMises(mean = 0, kappa = kappa1, n = n)
  X_axis_northpole<-cbind(sqrt(1-z^2)*cos(theta),sqrt(1-z^2)*sin(theta), z)
  
  cx<-(mu1-nu*mu0)/sqrt(1-nu^2)
  cy<-cross(mu0,cx)
  cz<-mu0
  
  #rotate
  X<-X_axis_northpole%*%rbind(cx,cy,cz)
  
  return(X)
}

# draw circle on unit sphere S2 by center of small circle and r
# converted code of Sungkyu Jung MATLAB drawCircleS2.m
drawCircleS2 <- function(center,theta) {
  # NB!!! theta is the angle from center
  if(theta==pi/2){
    t<-cbind(cos(seq(0,2*pi,length.out = 50)),sin(seq(0,2*pi,length.out = 50)),rep(0,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }else{
    t<-cbind(sin(theta)*cos(seq(0,2*pi,length.out = 50)),sin(theta)*sin(seq(0,2*pi,length.out = 50)),cos(theta)*rep(1,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }
  spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
  plot3d(sCirc,type="l",col = "black",expand = 10,box=TRUE,add = TRUE)
}


# copyright belongs to Sungkyu Jung
# converted code from MATLAB
# output is the type of fitted circle 
# likelihood ratio test from shapes package
LRTpval2 <- function(resGreat,resSmall,n) {
  chi2 <- max(n*log(sum(resGreat^2)/sum(resSmall^2)))
  pval <- 1-pchisq(q = chi2, df = 1, lower.tail = T) # likelihood test p-value Also you can use chi2cdf(chi2,1) from library(PEIP) like matlab
}

kurtosisTestFunction <- function(sphericalData, alpha=0.1) {
  ndata<-dim(sphericalData)[2]
  
  subsphereSmall<-getSubSphere(sphericalData,geodesic = "small")
  subsphereGreat<-getSubSphere(sphericalData,geodesic = "great")
  
  currentSphere<-sphericalData
  
  rSmall<-subsphereSmall$r                     # rSmall is rs in matlab
  centerSmall<-subsphereSmall$center           # NB! center is the centerSmall is pnsSmall$PNS$orthaxis[[1]]
  # and centers in matlab
  resSmall <- acos(t(centerSmall)%*%currentSphere)-rSmall  # NB!!! resSmall==(pnsSmall$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  rGreat<-subsphereGreat$r                     # rGreat is rg in matlab
  centerGreat<-subsphereGreat$center           # centerGreat is centers in matlab
  resGreat <- acos(t(centerGreat)%*%currentSphere)-rGreat  # NB!!! resGreat==(pnsGreat$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  # LRTpval is the likelihood ratio test from 'shapes' package
  # Chi-squared statistic for a likelihood test
  pval1 <- LRTpval(resGreat,resSmall,n = ndata)
  pval1
  
  if(pval1>alpha){
    print('great by likelihood ratio test')
    return('great')
    break
  }
  
  # # equivalently we can find pval by pns function
  # pnsTest2<-pns(sphericalData)
  # pnsTest2$PNS$pvalues
  # sum(pnsTest2$resmat[2,]==resSmall)
  
  # kurtosis test routine
  X <- LogNPd(rotMat(centerSmall) %*% currentSphere)
  
  # plot3d(t(sphericalData),type="p",expand=10, add=TRUE)
  # plot3d(t(rbind(X,rep(1,dim(X)[2]))),type="p",col = "blue",expand=10, add=TRUE)
  
  # Note that the tangential point is the center of the small circle
  d<-dim(X)[1]
  n<-dim(X)[2]
  normX2 <- colSums(X^2)
  kurtosis <- sum( normX2^2 ) / n / ( sum( normX2 ) / (d * (n-1)) )^2
  M_kurt <- d * (d+2)^2 / (d+4)
  V_kurt <- (1/n) * (128*d*(d+2)^4) / ((d+4)^3*(d+6)*(d+8))
  pval2 <- pnorm((kurtosis - M_kurt) / sqrt(V_kurt))
  
  if(pval2>alpha){
    return('great')
  }else{
    # drawCircleS2(normalVec = centerSmall,radius = rSmall)
    return('small')
  }
}

