setwd("C:/Users/Ben Archer/Documents/Notes/Third Year Notes/Diss Work/Dissertation_Work/FMNIST")

trainingset = read.csv("fashion-mnist_train.csv")
testset = read.csv("fashion-mnist_test.csv")

y_training = trainingset$label
x_training = trainingset[,2:785]

y_test = testset$label
x_test = testset[,2:785]

#Standardisation
x_training = x_training/255
y_training = y_training + 1

x_test = x_test/255
y_test = y_test + 1

#x_training = apply(x_training, 2, as.numeric)
#x_test = apply(x_test, 2, as.numeric)

x_training = as.matrix(x_training)
x_test = as.matrix(x_test)

ncols = length(x_training[1,])+1
training_length = length(x_training[,1])
test_length = length(x_test[,1])

draw = function(image){ #To draw each image
  m = matrix(image*255, 28, 28) 
  par(mar = c(0,0,0,0))
  m = apply(m, 2, as.numeric)
  m = apply(m, 1, rev)
  m = apply(m, 1, rev)
  image(z = m, useRaster = TRUE, axes = FALSE, col = gray.colors(256))
  par(mfrow=c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
}

input_size <- ncols-1 #define the network architecture
hidden_size <- c(20,20,20) 
output_size <- max(y_test)
sizes = c(input_size,hidden_size,output_size)
width = length(sizes)
