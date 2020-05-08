f = function(x){
  return(sin(3*x))
}

noisey = function(x){
  return(f(x) + rnorm(length(x),0,noises))
}

noises = 0.4

x = runif(3000,-3,3)
xseq = seq(-4,4,length = 10000)

y = noisey(x)

x_scaled = (x-mean(x))/sd(x)

dataset = cbind(x_scaled,y)

nrows = length(x)
ncols = 2

proportion = 0.7

training = dataset[1:round(proportion*nrows),]    #Take the first 90% of rows for taraining
x_training = as.matrix(training[,1:(ncols-1)])
y_training = training[,ncols]

test = dataset[round(proportion*nrows):nrows,]
x_test = as.matrix(test[,1:(ncols-1)])
y_test = test[,ncols]

input_size <- ncols-1 #define the network architecture
hidden_size <- c(9,9,9) 
output_size <- 1
sizes = c(input_size,hidden_size,output_size)
width = length(sizes)

plot(x_training*sd(x)+mean(x),y_training, col = 'blue')
points(x_test*sd(x)+mean(x),y_test,col = 'red')
lines(xseq,f(xseq))