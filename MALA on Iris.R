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

old_conditional_prob = function(output, i){
  count = c(0,0,0)
  for (j in 1:length(output)){
    max = which.max(output[[j]][i,])
    count[max] = count[max] + 1
  }
  prob = count/length(output)
  return(prob)
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
    } else{
      print(j)
    }
    correct[max,y[j]] = correct[max,y[j]] + 1
  }
  total_correct = total_correct/length(output[[1]][,1])
  return(list(correct,total_correct))
}

MALA = function(beta, x, y, iter, tau){
  t = 0
  acc = 0
  beta_list = matrix(beta,nrow = 1)
  N = length(x[,1])
  llh_grad_old = 0
  for (i in 1:N){
    llh_grad_old = llh_grad_old + backprop(beta, x[i,], y[i])
  }
  prior_grad_old = -beta
  total_grad_old = llh_grad_old + prior_grad_old
  
  for (i in 2:iter){
    t = t + 1
    betastar = beta_list[i-1,] + tau*total_grad_old + sqrt(2*tau)*rnorm(beta_length,0,1)
    #print((tau*total_grad_old + sqrt(2*tau)*1)[1:5])
    llh_grad_prop = 0
    for (j in 1:N){
      llh_grad_prop = llh_grad_prop + backprop(betastar, x[j,], y[j])
    }
    prior_grad_prop = - betastar
    total_grad_prop = llh_grad_prop + prior_grad_prop
    
    targets = exp(loglikelihood(betastar,x,y)+logprior(betastar))/exp(loglikelihood(beta_list[i-1,],x,y)+logprior(beta_list[i-1,]))
    #print(targets)
    proposals = proposal(beta_list[i-1,],betastar,total_grad_prop,tau)/proposal(betastar,beta_list[i-1,],total_grad_old,tau)
    #print(proposals)
    
    A = min(1,targets*proposals)
    #print(A)
    u = runif(1,0,1)
    
    if (u < A){
      beta_list = rbind(beta_list,betastar)
      total_grad_old = total_grad_prop
      acc = acc+1
    }
    else{
      beta_list = rbind(beta_list,beta_list[i-1,])
    }
  }
  return (list(beta_list,acc/iter))
}

breaks = 100
iter = 10000
tau = 0.000008
burn = 1000

initptm = proc.time()
output = MALA(init_beta(),x_training,y_training,iter/breaks,tau)
proc.time() - initptm
batch = output[[1]]
final_beta = batch[iter/breaks,]
acceptance_rate = output[[2]]
print(1)
for (i in 2:breaks){
  ptm = proc.time()
  new_output = MALA(final_beta,x_training,y_training,iter/breaks,tau)
  print(proc.time() - ptm)
  final_beta = new_output[[1]][iter/breaks,]
  acceptance_rate = c(acceptance_rate,new_output[[2]])
  batch = rbind(batch,new_output[[1]])
  print(i)
}
print(proc.time()-initptm)

mean(acceptance_rate)

output_training = outputs(batch[burn:iter,],x_training,y_training)
correct_training = accuracy(output_training, y_training)

output_test =  outputs(batch[burn:iter,],x_test,y_test)
correct_test = accuracy(output_test, y_test)


