
library(OpenImageR)
Im=OpenImageR::readImage("Training/1AT.jpg")
dim(Im)
red=Im[,,1]
green=Im[,,2]
blue=Im[,,3]
imageShow(red)
blue
# concatenate them in a long vector
as.vector(red)
youhavetodo=c(as.vector(red),as.vector(green), as.vector(blue))

# size of var covar matrix= 3x3, for principal components
m=cbind(as.vector(red),as.vector(green),as.vector(blue));m
dim(m)
out=princomp(m)
pc1=matrix(out$scores[,1],nrow = nrow(red),ncol = ncol(green))
imageShow(pc1)
pc2=matrix(out$scores[,2],nrow = nrow(red),ncol = ncol(green))
imageShow(pc2)
pc3=matrix(out$scores[,3],nrow = nrow(red),ncol = ncol(green))
imageShow(pc3)

out$sdev/sum(out$sdev)


###ACTUAL THING
library(dplyr)
red=c()
green=c()
blue=c()
Files = list.files(path="Training/");Files

for(i in seq_along(Files)){
  Im = readImage(paste0("Training/",Files[i]))
  ri=as.vector(Im[,,1])
  gi=as.vector(Im[,,2])
  bi=as.vector(Im[,,3])
  red=c(red,ri)
  green=c(green,gi)
  blue=c(blue,bi)
}


m=cbind(red,green,blue)
n = 36000
num<-c(10:19,1,20:25,2:9)
label=c()
for (i in num){
  temp = rep(i,n*6)
  label=c(label,temp)
}
m=cbind(m,label)
m

minkowski_distance <- function(x, y, p) {
  # Ensure both vectors have the same length
  if(length(x) != length(y)) {
    stop("Vectors x and y must have the same length.")
  }
  
  # Calculate Minkowski distance
  dist <- sum(abs(x - y)^p)^(1/p)
  return(dist)
}

