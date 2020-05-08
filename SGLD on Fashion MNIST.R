library(mcmc)
library(mcmcse)

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

out_delta = function(x,y){
  delta = out_fn(x)
  delta[y] = out_fn(x)[y]-1
  return(delta)
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

loglikelihood = function(beta, x, y){
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
  return (llh)
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

logprior = function(beta){
  return (sum(dnorm(beta, log = TRUE)))
}

proposal = function(betastar,beta,grad,tau){
  return (exp(-1/(4*tau)*sum((betastar-beta-tau*grad)^2)))
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

outputs = function(betalist, x, y){
  output_list = list()
  N = length(betalist[,1])
  for (j in 1:N){
    if (j%%round(N/100) == 0){
      print (j/N)
    }
    n = length(x[,1])
    beta = betaunflatten(betalist[j,])
    h_vals = list()
    h_vals[[1]] = act_fn(cbind(1,x) %*% beta[[1]])
    for (i in 2:(width-2)){
      h_vals[[i]] = act_fn(cbind(1,h_vals[[i-1]]) %*% beta[[i]])
    }
    h_vals[[width-1]] = cbind(1,h_vals[[width - 2]]) %*% beta[[width - 1]]
    output = t(apply(h_vals[[width-1]], 1, out_fn))
    output_list[[j]] = output
  }
  return (output_list)
}

conditional_prob = function(output, i){
  count = rep(0,10)
  N = length(output)
  for (j in 1:N){
    count = count + output[[j]][i,]/N
  }
  return(count)
}

accuracy = function(output,y){
  correct = matrix(rep(0,100),nrow = 10)
  total_correct = 0
  jlength = length(output[[1]][,1])
  for (j in 1:jlength){
    if (j%%round(jlength/10) == 0){
      print(j/jlength)
    }
    prob = conditional_prob(output, j)
    max = which.max(prob)
    if (max == as.numeric(y[j])){
      total_correct = total_correct+1
    } 
    correct[max,y[j]] = correct[max,y[j]] + 1
  }
  total_correct = total_correct/length(output[[1]][,1])
  return(list(correct,total_correct))
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

SGLD = function(beta, x, y, iter, n, t = 0, score = TRUE){
  N = length(x[,1])
  beta_list = matrix(beta,nrow = 1)
  if(score){
    test_loss = loss(beta,x_test,y_test)
    test_accuracy = tight_acc_check(beta,x_test,y_test)
    print(c(t,test_loss,test_accuracy))
    test_loss_list = c(test_loss)
    test_accuracy_list = c(test_accuracy)
  }
  
  for (i in 2:iter){
    t = t + 1
    row_indexs = sample(N,size = n, replace = FALSE)
    x_rows = x[row_indexs,]
    y_rows = y[row_indexs]
    total_grad = rep(0,beta_length)
    for (i in 1:n){
      total_grad = total_grad - backprop(beta, x_rows[i,], y_rows[i])
    }
    likelihood_grad = N*total_grad/n
    logprior_grad = -beta
    grad = logprior_grad + likelihood_grad
    beta = beta + epsilon(t)*grad + sqrt(2*epsilon(t))*rnorm(1,0,1)
    if (score == TRUE){
      if (t%%takes == 0){
        test_loss = loss(beta,x_test,y_test)
        test_accuracy = tight_acc_check(beta,x_test,y_test)
        print(c(t,test_loss,test_accuracy))
        test_loss_list = c(test_loss_list,test_loss)
        test_accuracy_list = c(test_accuracy_list,test_accuracy)
      }
      beta_list = rbind(beta_list,beta)
    }
  }
  return (beta_list)
}

iter = 1000
takes = 10000
divisions = 100
breaks = iter/divisions
burn = 4000
addon = 50000
timeson = 10000
decay = 1
batchsize = 30

init_t = 1

ptmbegin = proc.time()
output = SGLD(burn_beta,x_training,y_training,breaks,batchsize, t = init_t)
proc.time() - ptmbegin
final_beta = output[breaks,]

for (i in 2:divisions){
  #ptm = proc.time()
  new_output = SGLD(final_beta,x_training,y_training,breaks,batchsize, t = breaks*(i-1) + init_t)
  #print(proc.time() - ptm)
  final_beta = new_output[breaks,]
  output = rbind(output,new_output)
}
proc.time()-ptmbegin

print(dim(output))

#output_training = outputs(output,x_training,y_training)
#correct_training = accuracy(output_training, y_training)

output_test =  outputs(output,x_test,y_test)
correct_test = accuracy(output_test, y_test)

burn_beta = final_beta
