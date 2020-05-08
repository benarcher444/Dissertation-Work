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

out_fn <- function(x) { #Define Output Function
  exp(x)/sum(exp(x))
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
  return (llh/n)
}

logprior = function(beta){
  sum(dnorm(beta, log = TRUE))
}

logpost = function(beta){
  return (loglikelihood(beta,x_training,y_training)+logprior(beta))
}

outputs = function(betalist, x, y){
  output_list = list()
  N = length(betalist[,1])
  for (j in 1:N){
    if (j%%round(N/10) == 0){
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
  count = rep(0,3)
  for (j in 1:length(output)){
    count = count + output[[j]][i,]
  }
  prob = count/length(output)
  return(prob)
}

accuracy = function(output,y){
  correct = matrix(rep(0,9),nrow = 3)
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
    }else{
      print(j)
    }
    correct[max,y[j]] = correct[max,y[j]] + 1
  }
  total_correct = total_correct/length(output[[1]][,1])
  return(list(correct,total_correct))
}

random_walk = function(beta,N,sigma){
  a = 0
  beta_list = list()
  beta_list[[1]] = beta
  for (i in 2:N){
    proposed_beta = beta + rnorm(beta_length,0,sigma)
    A = exp(logpost(proposed_beta))/exp(logpost(beta))
    if (runif(1) < A) {
      beta = proposed_beta
      a = a + 1
    }
    beta_list[[i]] = beta
    if (i%%100 == 0){
      print(i/N)
    }
  }
  print(a/N)
  return (beta_list)
}

burn_in = 10000
iter = 100000
blen = 1
sigma = 0.1

burn = random_walk(init_beta(),burn_in,sigma/2)

final_beta = burn[[burn_in]]

beta_list = random_walk(final_beta,iter,sigma)