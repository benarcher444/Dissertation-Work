proposal = function(x,y,xp,yp){
  return (exp(-1/(4*tau)*sum((c(x,y)-c(xp,yp)-tau*gradlogtarget(x,y))^2)))
}

gradlogtarget = function(x,y,h = 0.000001){
  xgrad = (log(target(x+h,y))-log(target(x,y)))/h
  ygrad = (log(target(x,y+h))-log(target(x,y)))/h
  return(c(xgrad,ygrad))
}

logtarget = function(x,y){
  return (log(target(x,y)))
}

contourplotout = contourplot(logtarget, n = 1000, minx = -2, maxx = 2, miny = -2, maxy = 2)

ldens = contourplotout[[1]]
l1 = contourplotout[[2]]
l2 = contourplotout[[3]]

N = 1000000

x = rep(0,N)
y = rep(0,N)

X = cbind(x,y)

draw = FALSE
set.seed(42)
x[1] = 0
y[1] = 0
a = 0
s = 1
t = 0.08
par(mfrow = c(1,1))

if (draw){
  contour(l1,l2,tdens)
}

for(i in 2:N){
  current_x = x[i-1]
  current_y = y[i-1]
  grad = gradlogtarget(current_x,current_y)
  rn1 = rnorm(1,0,1)
  rn2 = rnorm(1,0,1)
  proposed_x = current_x + t*grad[1] + sqrt(2*t)*rn1
  proposed_y = current_y + t*grad[2] + sqrt(2*t)*rn2
  
  A = min(1,target(proposed_x,proposed_y)*proposal(current_x,current_y,proposed_x,proposed_y)/(target(current_x,current_y)*proposal(proposed_x,proposed_y,current_x,current_y)))
  #print(A)
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
  if (draw){
    lines(c(current_x,current_x + t*grad[1]),c(current_y,current_y + t*grad[2]), col = 'red')
    lines(c(current_x + t*grad[1],current_x + t*grad[1] + sqrt(2*t)*rn1),c(current_y + t*grad[2],current_y + t*grad[2] + sqrt(2*t)*rn2), col = 'blue')
    lines(c(x[i-1],x[i]),c(y[i-1],y[i]))
    points(x[i],y[i])
    Sys.sleep(0.2)
  }
}

n = 20
z=table(cut(x,n),cut(y,n))
print(a/N)

#par(mfrow = c(1,3))
#contour(l1,l2,tdens, main = 'Contour Plot of Target Distribution', xlab = 'x', ylab = 'y')
#plot(x[seq(0,N,length = N/100)],y[seq(0,N,length = N/100)], xlab = 'x', ylab = 'y', main = 'Points Sampled')
#image2D(z=table(cut(x,n),cut(y,n)),border = 'black', x = seq(-2,2,length = n), y = seq(-2,2,length = n), xlab = 'x', ylab = 'y', main = 'Heat Map of Samples')

phi = 10
theta = 20
par(mfrow = c(1,3))
persp(l1,l2,tdens, phi = phi, theta = theta, main = '3D Contour Plot of Target Distribution')
persp(l1,l2,ldens, phi = phi, theta = theta, main = '3D Contour Plot of Log-Target Distribution')
hist3D(z=z, border = 'black', phi = phi, theta = theta, main = '3D Histogram of Sampled Points')

n = 100
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

acceptance = a/N
acceptance
