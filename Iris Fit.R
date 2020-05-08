library(mcmc)

data(iris)

dataset = iris #load in the data set

nrows = length(dataset[,1])
ncols = length(dataset[1,])

head(dataset)
rows = sample(nrows)
dataset = dataset[rows,]  #Shuffle the data set and order it x's then y's

proportion = 0.7

training = dataset[1:round(proportion*nrows),]    #Take the first 70% of rows for training
x_training = as.matrix(training[,1:(ncols-1)])
y_training = training[,ncols] 

test = dataset[round(proportion*nrows):nrows,]
x_test = as.matrix(test[,1:(ncols-1)])
y_test = test[,ncols] 

input_size <- ncols-1 #define the network architecture
hidden_size <- c(10,10,10) 
output_size <- 3
sizes = c(input_size,hidden_size,output_size)
width = length(sizes)

pairs(iris[,c(1,2,3,4,5)], col = iris[,5],pch = c(rep(1,72),3,rep(1,10),3,rep(1,66)))

