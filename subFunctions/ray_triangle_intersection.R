
# return the intersection of a ray with a triangle if the intersection exist
rayTriangleIntersection <- function(rayOrigin,rayDirection,triangleVertex1,triangleVertex2,triangleVertex3) {
  
  O<-rayOrigin #origin of the ray
  D<-rayDirection #direction of the ray 
  A<-triangleVertex1#triangle vertices
  B<-triangleVertex2
  C<-triangleVertex3
  
  E1<-B-A
  E2<-C-A
  N<-myCrossProduct(E1,E2)
  
  det<-(-sum(D*N))
  invdet <- 1/det
  AO<-O-A
  DAO<-myCrossProduct(AO,D)
  u<-sum(E2*DAO)*invdet
  v<-(-sum(E1*DAO)*invdet)
  t<-sum(AO*N)*invdet
  if (abs(det) >= 1e-6 & t >= 0 & u >= 0 & v >= 0 & (u+v) <= 1){
    intersection<-O + t * D
  }else{
    intersection<-c(NA,NA,NA)
  }
  return(intersection)
}



