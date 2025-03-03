
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
Files = list.files(path="Training/");Files
#m = vector("list", length(Files))
m = matrix(NA, nrow=150, ncol = 108000)

for(i in seq_along(Files)){
  Im = readImage(paste0("Training/",Files[i]))
  ri=as.vector(Im[,,1])
  gi=as.vector(Im[,,2])
  bi=as.vector(Im[,,3])
  m[i,] = cbind(ri,gi,bi)
}
dim(m) ## m == X in formula
m[1,1]
obs_mean = colMeans(m)
m2 = m-colMeans(m) #subtracting every mean from every row
n=nrow(m2)
m_ =t(m2)
Sigma_s= (m2%*%m_)/(n-1)
dim(Sigma_s)
eig=eigen(Sigma_s)
e_vec_s=eig$vectors
e_vec_l = m_%*%e_vec_s #== P in formula
dim(e_vec_l)
e_vals=diag(eig$values)
round(eig$values/sum(eig$values),3)

PCAs=e_vec_l[1:10,]%*%m2
PCAs=m2%*%e_vec_l

dim(PCAs)
# n = 36000
num<-c(10:19,1,20:25,2:9)
label=c()
for (i in num){
  temp = rep(i,6)
  label=c(label,temp)
}
M=cbind(as.factor(label),m2)






###FISHER SHIT
# in summary, you get Sb, Sw and then multiply them to get w, then with w you can get the fisher value
# n == 150, p == 108000
# HERE'S THE TRICK, we compare each one of our classes to the OVERALL means and such
# or rather, to the ALTERNATIVE MEAN, so all the classes EXCEPT the one we're looking at
# we have classes c1-c25
num<-c(10:19,1,20:25,2:9)

cmeans = matrix(NA, nrow=25, ncol = 108000)
# get class means
mcopy = m
for (i in num){
  cmeans[i,] = colMeans(mcopy[1:6,])
  mcopy = mcopy[-c(1:6),]
}
cmeans[1,]

# get external means
extmeans = matrix(NA, nrow=25, ncol = 108000)

for (i in num){
  mcopy = m
  start = 6*(i-1)#step size
  mcopy = mcopy[-c(seq(start+1,start+6,1)),] #delete targeted class
  extmeans[i,] = colMeans(mcopy)
}
extmeans[1,]

# Get Sb
for (i in num){
  diff[i,j] = cmeans[j,] - cmeans[i,]
}




















