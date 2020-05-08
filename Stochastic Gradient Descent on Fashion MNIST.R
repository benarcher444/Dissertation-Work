init_zeros = function(){ #Creates a matrix in the right format
  beta <- list()
  for (i in 1:(width-1)){
    beta[[i]] = matrix(rep(0,(sizes[i]+1)*sizes[i+1]), nrow = sizes[i]+1, ncol = sizes[i+1])
  }
  return (beta)
}

beta_shape = init_zeros()

betaunflatten = function(betalist){   ####Functions for converting beta between matrix and just vector form for ease of use
  beta_matrices = list()
  upto = 0
  for(i in 1:(width-1)) {
    beta_matrices[[i]] <- matrix(betalist[(upto+1):(upto+length(beta_shape[[i]]))],
                                 nrow = nrow(beta_shape[[i]]),
                                 ncol = ncol(beta_shape[[i]]))
    upto <- upto + length(beta_shape[[i]])
  }
  return (beta_matrices)
}

betaflatten = function(beta){
  return (c(beta, recursive = TRUE))
}

beta_length = length(betaflatten(beta_shape))

init_beta = function(){
  beta = list()
  for ( i in 1:(width-1)){
    beta[[i]] = runif((sizes[i]+1)*sizes[i+1],-sqrt(6/(sizes[i]+sizes[i+1])),sqrt(6/(sizes[i]+sizes[i+1])))
  }
  return (c(beta,recursive = TRUE))
}

act_fn <- function(x) {  #Define activation function
  ifelse(x > 0, x, x/1000)
}

act_fn_grad = function(x){
  ifelse(x > 0, 1, 1/1000)
}

out_fn <- function(x) { #Define Output Function
  exp(x)/sum(exp(x))
}

yi_fn = function(x,y){
  (exp(x[y])-sum(exp(x)))/sum(exp(x))
}

out_delta = function(x,y){
  delta = out_fn(x)
  delta[y] = yi_fn(x,y)
  return(delta)
}

acc_check = function(beta,x,y){
  beta = betaunflatten(beta)
  n = length(y)                  
  hidden_val = list()
  hidden_val[[1]] = act_fn(cbind(1,x) %*% beta[[1]])
  for (i in 2:(width-2)){
    hidden_val[[i]] = act_fn(cbind(1,hidden_val[[i-1]]) %*% beta[[i]])
  }
  output_val = cbind(1,hidden_val[[width - 2]]) %*% beta[[width - 1]]
  output = t(apply(output_val, 1, out_fn))
  correct = 0
  whichmax = c()
  ifcorrect = c()
  for (i in 1:dim(x)[1]){
    whichmax = c(whichmax,which.max(output[i,]))
    ifcorrect = c(ifcorrect,whichmax[i] == as.numeric(y[i]))
    if (ifcorrect[i]){
      correct = correct + 1
    }
  }
  return(cbind(output,whichmax,y,ifcorrect,correct/n))
}

tight_acc_check = function(beta,x,y){
  beta = betaunflatten(beta)
  n = length(y)                  
  hidden_val = list()
  hidden_val[[1]] = act_fn(cbind(1,x) %*% beta[[1]])
  for (i in 2:(width-2)){
    hidden_val[[i]] = act_fn(cbind(1,hidden_val[[i-1]]) %*% beta[[i]])
  }
  output_val = cbind(1,hidden_val[[width - 2]]) %*% beta[[width - 1]]
  output = t(apply(output_val, 1, out_fn))
  correct = 0
  whichmax = c()
  ifcorrect = c()
  for (i in 1:dim(x)[1]){
    whichmax = c(whichmax,which.max(output[i,]))
    ifcorrect = c(ifcorrect,whichmax[i] == as.numeric(y[i]))
    if (ifcorrect[i]){
      correct = correct + 1
    }
  }
  return(correct/n)
}

epsilon = function(n){ 
  return (1/((n*decay+addon)*timeson))    
}

ah_unflatten = function(a){
  newa = list()
  upto = 0
  for (i in 1:(width-1)){
    newa[[i]] = matrix(a[(upto+1):(upto+sizes[i+1])], ncol = sizes[i+1])
    upto = upto + sizes[i+1]
  }
  return (newa)
}

loss = function(beta, x, y){
  n = length(x[,1])
  beta = betaunflatten(beta)
  h_vals = list()
  h_vals[[1]] = act_fn(cbind(1,x) %*% beta[[1]])
  for (i in 2:(width-2)){
    h_vals[[i]] = act_fn(cbind(1,h_vals[[i-1]]) %*% beta[[i]])
  }
  h_vals[[width-1]] = cbind(1,h_vals[[width - 2]]) %*% beta[[width - 1]]
  output = t(apply(h_vals[[width-1]], 1, out_fn))
  llh = 0
  for (i in 1:n){ 
    llh = llh + log(as.numeric(output[i,y[i]]))
  }
  return (-llh/n)
}

forward = function(beta, x){
  a_vals = list()
  h_vals = list()
  a_vals[[1]] = c(1,x) %*% beta[[1]]
  h_vals[[1]] = act_fn(a_vals[[1]])
  for (i in 2:(width-2)){
    a_vals[[i]] = cbind(1,h_vals[[i-1]]) %*% beta[[i]]
    h_vals[[i]] = act_fn(a_vals[[i]])
  }
  a_vals[[width-1]] = cbind(1,h_vals[[width - 2]]) %*% beta[[width - 1]]
  h_vals[[width-1]] = out_fn(a_vals[[width-1]])
  return (c(a_vals,h_vals, recursive = TRUE))
}

backprop = function(beta, x, y){
  
  beta = betaunflatten(beta)
  ah_vals = forward(beta,x)
  
  a_vals_full_list = ah_vals[1:(length(ah_vals)/2)]
  h_vals_full_list = ah_vals[(length(ah_vals)/2+1):length(ah_vals)]
  
  a_vals = ah_unflatten(a_vals_full_list)
  h_vals = ah_unflatten(h_vals_full_list)
  
  deltas = list()
  grad = list()
  
  deltas[[width-1]] = matrix(out_delta(a_vals[[width-1]],y))
  grad[[width-1]] = as.matrix(c(1,h_vals[[width-2]]))%*%t(deltas[[width-1]])
  
  for (i in (width-2):2){
    deltas[[i]] = t(act_fn_grad(a_vals[[i]]))*beta[[i+1]][-1,]%*%deltas[[i+1]]
    grad[[i]] = as.matrix(c(1,h_vals[[i-1]]))%*%t(matrix(deltas[[i]]))
  }
  
  deltas[[1]] = t(act_fn_grad(a_vals[[1]]))*beta[[2]][-1,]%*%deltas[[2]]
  grad[[1]] = as.matrix(c(1,x))%*%t(matrix(deltas[[1]]))
  return (c(grad, recursive = TRUE))
}

confcalc = function(beta,x,y){
  beta = betaunflatten(beta)
  n = length(y)                  
  hidden_val = list()
  hidden_val[[1]] = act_fn(cbind(1,x) %*% beta[[1]])
  for (i in 2:(width-2)){
    hidden_val[[i]] = act_fn(cbind(1,hidden_val[[i-1]]) %*% beta[[i]])
  }
  output_val = cbind(1,hidden_val[[width - 2]]) %*% beta[[width - 1]]
  output = t(apply(output_val, 1, out_fn))
  correct = 0
  whichmax = c()
  ifcorrect = c()
  confmat = matrix(rep(0,100),nrow = 10)
  for (i in 1:dim(x)[1]){
    confmat[which.max(output[i,]),y[i]] = confmat[which.max(output[i,]),y[i]] + 1
    whichmax = c(whichmax,which.max(output[i,]))
    ifcorrect = c(ifcorrect,whichmax[i] == as.numeric(y[i]))
    if (ifcorrect[i]){
      correct = correct + 1
    }
  }
  return(list(confmat,correct/n))
}

stochasticgd = function(x, y, n, iter){
  t = 0
  beta = init_beta()
  N = length(x[,1])
  
  if (TRUE){
    training_loss = loss(beta,x,y)
    test_loss = loss(beta,x_test,y_test)
    #training_accuracy = tight_acc_check(beta,x,y)
    #test_accuracy = tight_acc_check(beta,x_test,y_test)
    print(c(t,test_loss))
    training_loss_list = c(training_loss)
    test_loss_list = c(test_loss)
    #training_accuracy_list = c(training_accuracy)
    #test_accuracy_list = c(test_accuracy)
  }
  
  for (i in 1:iter){
    t = t + 1
    row_indexs = sample(N,size = n, replace = TRUE)
    x_rows = x[row_indexs,]
    y_rows = y[row_indexs]
    total_grad = rep(0,beta_length)
    for (i in 1:n){
      total_grad = total_grad + backprop(beta, x_rows[i,], y_rows[i])
    }
    grad = total_grad/n
    beta = beta - epsilon(t)*grad
    
    if (t%%takes == 0){
      training_loss = loss(beta,x,y)
      test_loss = loss(beta,x_test,y_test)
      #training_accuracy = tight_acc_check(beta,x,y)
      #test_accuracy = tight_acc_check(beta,x_test,y_test)
      print(c(t,test_loss))
      training_loss_list = c(training_loss_list,training_loss)
      test_loss_list = c(test_loss_list,test_loss)
      #training_accuracy_list = c(training_accuracy_list,training_accuracy)
      #test_accuracy_list = c(test_accuracy_list,test_accuracy)
    }
  }
  return (list(beta,training_loss_list,test_loss_list))
}


takes = 1000
addon = 10000
timeson = 0.001
decay = 1
n = 30
iter = 10000

ptm = proc.time()
output = stochasticgd(x_training,y_training,n,iter)
print(proc.time()-ptm)

final_beta = output[[1]]
training_loss_list = output[[2]]
test_loss_list = output[[3]]

plot(seq(0,iter,takes),training_loss_list, col = 'blue', type = 'l', lwd = 2, main = 'Stochastic Gradient Descent on Fashion MNIST', ylab = 'Cross Entropy Loss', xlab = 'Iteration')
lines(seq(0,iter,takes),test_loss_list, col = 'red', lwd = 2)

acc_check(final_beta,x_training,y_training)[1:10,]

training_loss = loss(final_beta,x_training,y_training)
training_loss
