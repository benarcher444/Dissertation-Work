contourplot = function(target,n = 50,minx=-2,maxx=2,miny=-2,maxy=2){
  n = 31
  dens = matrix(0,n,n)
  l1 = seq(minx,maxx,length = n)
  l2 = seq(miny,maxy,length = n)
  for (i in 1:n){
    for (j in 1:n){
      dens[i,j] = target(l1[i],l2[j])
    }
  }
  return (list(dens,l1,l2))
}

MVND = function(M,S,SS,x,y){
  X = c(x,y)
  return (as.numeric(1/(2*pi*sqrt(det(S)))*exp(-1/2*(X-M)%*%SS%*%(X-M))))
}

M1 = c(-0.9,-0.8)
S1 = matrix(c(0.1,0,0,0.1),nrow = 2)
SS1 = solve(S1)
W1 = 3
M2 = c(-0.8,0.8)
S2 = matrix(c(0.2,0,0,0.2),nrow = 2)
SS2 = solve(S2)
W2 = 3
M3 = c(1,0.5)
S3 = matrix(c(0.3,0,0,0.3),nrow = 2)
SS3 = solve(S3)
W3 = 3

target = function(x,y){
  return (1*MVND(M1,S1,SS1,x,y)/W1+1*MVND(M2,S2,SS2,x,y)/W2+1*MVND(M3,S3,SS3,x,y)/W3)
}

contourplotout = contourplot(target, n = 1000, minx = -2, maxx = 2, miny = -2, maxy = 2)

tdens = contourplotout[[1]]
l1 = contourplotout[[2]]
l2 = contourplotout[[3]]

par(mfrow = c(1,2))
contour(l1,l2,tdens, main = 'Contour Plot of the Target Distribution', xlab = 'x', ylab = 'y')
persp(l1,l2,tdens, phi = 10, theta  = 20, xlab = 'x', ylab = 'y', main = '3D Perspective of the Target Distribution', zlab = 'Density')