library(plot3D)

proposal = function(x,y){
  X = c(x,y)
  return (as.numeric(1/(2*pi*sqrt(det(S)))*exp(-1/2*(X-M)%*%solve(S)%*%(X-M))))
}

set.seed(42)
mx = 0
my = 0
sx = 0.9
sy = 0.9
p = 0

S = matrix(c(sx**2,p*sx*sy,p*sx*sy,sy**2),ncol = 2)
M = c(mx,my)

N = 1000000

x = rep(0,N)
y = rep(0,N)
a = 0

x[1] = 0
y[1] = 0

contourplotout = contourplot(proposal)

pdens = contourplotout[[1]]
l1 = contourplotout[[2]]
l2 = contourplotout[[3]]

for (i in 2:N){
  ptm = proc.time()
  current_x = x[i-1]
  current_y = y[i-1]
  z = rnorm(2,0,1)
  proposed_x = mx + sx*z[1]
  proposed_y = my + sy*(p*z[1] + sqrt(1-p**2)*z[2])
  A = target(proposed_x,proposed_y)*proposal(current_x,current_y)/(target(current_x,current_y)*proposal(proposed_x,proposed_y))
  if (runif(1) < A) {
    x[i] = proposed_x
    y[i] = proposed_y
    a = a + 1
  }
  else{
    x[i] = current_x
    y[i] = current_y
  }
  if (i%%10000 == 0){
    print(i/N)
  }
}

n = 50
z=table(cut(x,n),cut(y,n))
print(a/N)

par(mfrow = c(1,3))
contour(l1,l2,tdens, main = 'Contour Plot of Target Distribution', xlab = 'x', ylab = 'y')
contour(l1,l2,pdens, main = 'Contour Plot of Proposal Distribution with Sampled Points', xlab = 'x', ylab = 'y')
points(x[seq(0,N,length = N/100)],y[seq(0,N,length = N/100)])
image2D(z=z,border = 'black',x = seq(-2,2,length = n), y = seq(-2,2,length = n),xlab = 'x', ylab = 'y', main = 'Heat Map of Samples')

phi = 10
theta = 20
persp(l1,l2,tdens, phi = phi, theta = theta, main = '3D Contour Plot of Target Distribution',xlab = 'x', ylab = 'y', zlab = 'Density')
persp(l1,l2,pdens, phi = phi, theta = theta, main = '3D Contour Plot of Proposal Distribution',xlab = 'x', ylab = 'y', zlab = 'Density')
hist3D(z=z, border = 'black', phi = phi, theta = theta, main = '3D histogram of Sampled Points',xlab = 'x', ylab = 'y', zlab = 'Density')

par(mfrow = c(1,2))
hist(x, probability = TRUE, breaks = 100,xlab = 'x', ylab = 'Density')
lines(xseq,Xtotal, type = 'l', col = 'red')
hist(y, probability = TRUE, breaks = 100,xlab = 'x', ylab = 'Density')
lines(xseq,ytotal, type = 'l', col = 'red')

par(mfrow = c(1,4))
plot((N-1000):N,x[(N-1000):N],xlab = 'Iteration', ylab = 'X Value',main ='X Value Sampled At Each Iteration',type ='l')
plot((N-1000):N,y[(N-1000):N],xlab = 'Iteration', ylab = 'Y Value',main ='Y Value Sampled At Each Iteration',type ='l')

acf(x, main = 'X Samples Autocorrelation')
acf(y, main = 'Y Samples Autocorrelation')
