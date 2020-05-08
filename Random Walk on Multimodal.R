N = 1000000

x = rep(0,N)
y = rep(0,N)

set.seed(42)
x[1] = 0
y[1] = 0
a = 0
s = 1.5

for(i in 2:N){
  current_x = x[i-1]
  proposed_x = current_x + rnorm(1,0,s)
  current_y = y[i-1]
  proposed_y = current_y + rnorm(1,0,s)
  A = target(proposed_x,proposed_y)/target(current_x,current_y)
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

n = 100
z=table(cut(x,n),cut(y,n))
print(a/N)

par(mfrow = c(1,3))
contour(l1,l2,tdens, main = 'Contour Plot of Target Distribution', xlab = 'x', ylab = 'y')
plot(x[seq(0,N,length = N/100)],y[seq(0,N,length = N/100)], xlab = 'x', ylab = 'y', main = 'Points Sampled')
image2D(z=table(cut(x,n),cut(y,n)),border = 'black', x = seq(-2,2,length = n), y = seq(-2,2,length = n), xlab = 'x', ylab = 'y', main = 'Heat Map of Samples')

phi = 10
theta = 20
par(mfrow = c(1,2))
persp(l1,l2,tdens, phi = phi, theta = theta, main = '3D Contour Plot of Target Distribution')
hist3D(z=z, border = 'black', phi = phi, theta = theta, main = '3D Histogram of Sampled Points')

par(mfrow = c(1,2))
hist(x, probability = TRUE, breaks = 100,xlab = 'x', ylab = 'Density')
lines(xseq,Xtotal, type = 'l', col = 'red')
hist(y, probability = TRUE, breaks = 100,xlab = 'x', ylab = 'Density')
lines(xseq,ytotal, type = 'l', col = 'red')

par(mfrow = c(1,4))
plot((N-1000 across t):N,x[(N-1000):N],xlab = 'Iteration', ylab = 'X Value',main ='X Value Sampled At Each Iteration',type ='l')
plot((N-1000):N,y[(N-1000):N],xlab = 'Iteration', ylab = 'Y Value',main ='Y Value Sampled At Each Iteration',type ='l')

acf(x, main = 'X Samples Autocorrelation')
acf(y, main = 'Y Samples Autocorrelation')

acceptance = a/N
acceptance
