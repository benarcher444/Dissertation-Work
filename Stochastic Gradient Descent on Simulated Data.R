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
  return (rnorm(beta_length,0,stand))
}

act_fn <- function(x) {  #Define activation function
  ifelse(x > 0, x, x/1000)
}

act_fn_grad = function(x){
  ifelse(x > 0, 1, 1/1000)
}

out_fn <- function(x) { #Define Output Function
  x
}

out_fn_grad = function(x){
  ifelse(x>0,1,1)
}

epsilon = function(n){ 
  return (1/(n*decay+addon))    
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

plotter = function(beta, x){
  beta = betaunflatten(beta)
  h_vals = list()
  h_vals[[1]] = act_fn(cbind(1,x) %*% beta[[1]])
  for (i in 2:(width-2)){
    h_vals[[i]] = act_fn(cbind(1,h_vals[[i-1]]) %*% beta[[i]])
  }
  h_vals[[width-1]] = out_fn(cbind(1,h_vals[[width - 2]]) %*% beta[[width - 1]])
  return (h_vals[[width-1]])
}
  
loss = function(beta, x, y){
  n = length(x)
  beta = betaunflatten(beta)
  h_vals = list()
  h_vals[[1]] = act_fn(cbind(1,x) %*% beta[[1]])
  for (i in 2:(width-2)){
    h_vals[[i]] = act_fn(cbind(1,h_vals[[i-1]]) %*% beta[[i]])
  }
  h_vals[[width-1]] = out_fn(cbind(1,h_vals[[width - 2]]) %*% beta[[width - 1]])
  return (sum((h_vals[[width-1]]-y)^2)/n)
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
  
  deltas[[width-1]] = 2*(h_vals[[width-1]]-y)*out_fn_grad(a_vals[[width-1]])
  grad[[width-1]] = as.matrix(c(1,h_vals[[width-2]]))%*%matrix(deltas[[width-1]])
  
  for (i in (width-2):2){
    deltas[[i]] = t(act_fn_grad(a_vals[[i]]))*beta[[i+1]][-1,]%*%deltas[[i+1]]
    grad[[i]] = as.matrix(c(1,h_vals[[i-1]]))%*%t(matrix(deltas[[i]]))
  }
  
  deltas[[1]] = t(act_fn_grad(a_vals[[1]]))*beta[[2]][-1,]%*%deltas[[2]]
  grad[[1]] = as.matrix(c(1,x))%*%t(matrix(deltas[[1]]))
  
  return (c(grad, recursive = TRUE))
}

stochasticgd = function(x, y, n, iter, alpha){
  t = 0
  beta = init_beta()
  N = length(x)
  training_loss = loss(beta,x,y)
  test_loss = loss(beta,x_test,y_test)
  print(c(t,training_loss,test_loss))
  training_loss_list = c(training_loss)
  test_loss_list = c(test_loss)
  for (i in 1:iter){
    t = t + 1
    row_indexs = sample(N,size = n, replace = TRUE)
    x_rows = x[row_indexs]
    y_rows = y[row_indexs]
    total_grad = rep(0,beta_length)
    for (i in 1:n){
      total_grad = total_grad + backprop(beta, x_rows[i], y_rows[i])
    }
    #print(total_grad)
    grad = total_grad/n
    beta = beta - epsilon(t)*grad
    if (t%%100 == 0){
        training_loss = loss(beta,x,y)
        test_loss = loss(beta,x_test,y_test)
        print(c(t,training_loss,test_loss))
        training_loss_list = c(training_loss_list,training_loss)
        test_loss_list = c(test_loss_list,test_loss)

    }
  }
  return (list(beta,training_loss_list,test_loss_list))
}


stand = 1
addon = 300
decay = 0.8
n = 30
iter = 100000
alpha = 0.0001

ptm = proc.time()
burn = stochasticgd(x_training,y_training,n,1,alpha)
expected_time = (proc.time() - ptm)*iter
print(expected_time)

ptm = proc.time()
output = stochasticgd(x_training,y_training,n,iter,alpha)
print(proc.time()-ptm)

final_beta = output[[1]]
training_loss_list = output[[2]]
test_loss_list = output[[3]]

plot(seq(0,iter,100),training_loss_list, col = 'blue', ylim = c(0,1), type = 'l', lwd = 2)
lines(seq(0,iter,100),test_loss_list, col = 'red', lwd = 2)
abline(h = noises**2, lty = 3)

xseq = seq(-3,3,length = 1000)
scaledxseq = (xseq-mean(x))/sd(x)
plot(xseq,plotter(final_beta,scaledxseq),type = 'l', lwd = 2, col = 'red', ylim = c(-2,2.5), main = 'Medium Neural Network Fit on Unscaled Data Set',ylab = 'y',xlab = 'x')
points(x_training*sd(x) + mean(x),y_training)
points(x_test*sd(x) + mean(x),y_test,col = 'blue')
legend(-3,2.4,legend = c('Network Fit','sin(3x)'), col = c('red','black'), lty = c(1,1), lwd = rep(2,2))
lines(xseq,f(xseq), type = 'l', lwd = 3)
lines(xseq,plotter(final_beta,scaledxseq),type = 'l', lwd = 3, col = 'red')

takes  = 50
plot(seq(0,iter,100*takes),training_loss_list[seq(1,iter/100+1,takes)], col = 'blue', ylim = c(0.05,0.4), type = 'l', lwd = 2, main = 'Loss Graph', ylab = 'Mean Squared Error Loss Function Value', xlab = 'Iteration')
lines(seq(0,iter,100*takes),test_loss_list[seq(1,iter/100+1,takes)], col = 'red', lwd = 2)
abline(h = noises**2, lty = 3,lwd = 2)
legend(6.7e+4,0.402,legend = c('Training Loss','Test Loss','Sigma Squared'), col = c('blue','red','black'), lty = c(1,1,3), lwd = rep(2,2,2))
