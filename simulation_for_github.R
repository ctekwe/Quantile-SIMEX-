library(stats)
library(sn)
library(splines)
library(ald)
library(mvtnorm)
library(boot)
library(fda)
library(ggplot2)
library(lattice)
library(readxl)
library(plyr)
library(magrittr)
library(Matrix)
library(fda)
library(SparseM)
library(quantreg)
library(corpcor)
library(MASS)
library(reshape)
library(grid)
library(gridExtra)
library(knitr)
library(doParallel)
library(foreach)


### extrapolation function 
extra_quad = function(k,beta_hat_lambda,lambda,lambda_app){
  ## k is the number of basis functions for the functional variable
  ## beta_hat_lambda is the estimated beta(t) from SIMEX simulation step
  ## lambda is the squence of lambda 
  ## lambda_app = -1 for SIMEX estimators 
  gamma = rep(0,k)  
  lambda2 = lambda*lambda
  for(m in 1:k)
  {
    testr<-beta_hat_lambda[m,]
    lr<-lm(testr~lambda+lambda2)
    coeff<-lr$coefficients
    gamma[m]<- coeff[1]+lambda_app*coeff[2]+lambda_app*lambda_app*coeff[3]
  }
  
  return(gamma)
}

#### functions used to generate data 
met = c("f1","f2","f3","f4","f5", "f6", "f7")
Betafunc = function(t, met){
  switch(met,
         f1 = sin(2*pi*t),  ## f1 = sin(2*pi*t),
         f2 = 1/(1+exp(4*2*(t-.5))),
         f3 = sin(pi*(8*(t-.5))/2)/ (1 + (2*(8*(t-.5))^2)*(sign(8*(t-.5))+1)),
         f4 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(.5*sin(8*(t-.5))+1)),
         f5 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(sin(8*(t-.5))+1)),
         f6 = 1/(1+exp((t))),
         f7 = 1/(1+exp(4*2*(t-.5))) +1
  )
}


######## function for selecting number of basis ###
### based on BIC 
### the output of this function is a BIC value
L_BIC=function(k, Y, W, tau, EF=EF,time_interval){
  ### Y is response, W is observed matrix of W(t); 
  ### k is the number of basis and should be at least 4; 
  ##  tau is quantile level; 
  ## time_interval is interval of observed time
  cubic_bs = bs(time_interval, df=k, degree=3, intercept = TRUE)
  pred = t(t(W)-colMeans(W))%*%cubic_bs/length(time_interval)
  
  coef_quant = rq(Y ~ pred +EF, tau=tau)$coefficients
  n = length(Y)
  eps = Y-cbind(rep(1,n),pred, EF)%*%coef_quant
  loss = eps*ifelse(eps<0, tau-1, tau) 
  
  re=log(mean(loss))+(k+1)*log(n)/n
  return(re)
  
}



main.delta=function(n,c,sd_x,sd_w,p_sim,fi,fx,seeds,sim_liter){
  require(fda)
  ## @seeds is the the seeding basis
  ## @n is the size of simulated data
  ## @c is the tuning parameter for generating delta(t) (delta(t) = c*sin(2*pi*t) +1)
  ## @sd_x is the sd of X(t)
  ## @sd_w is the sd of W(t)
  ## @p_sim is the quantile number (0.25, 0.5, 0.75, 0.95)
  ## @fi is the function for generating Y
  ## @fx is the function for generating X(t)
  ## @sim_lier is the number of replications we do for one simulation
  
  set.seed(seeds)
  t = 100  ## Number of time points
  k_max = 7 ## The maximum number of basis functions to be considered
  a  = seq(0, 1, length.out=t) #create a sequence between 0,1 of length = t
  
  sd_m=0.25 ## standard deviation of instrumental variable 
  sigma_y = .1 ###standard deviation of Y
  sigma_ef = 0.5 ###standard deviation of EF covariate
  p = 0.6 ## prob for the EF binary predictor
 
  ############ create covariance matrix for error terms #################
  rho_x=0.5 ## correlation coefficient of X(t)
  rho_w=0.5 ## correlation coefficient of W(t)
  rho_m=0.5 ## correlation coefficient of M(t)
  
  ###compound symmetric covariance matrix ###
  Sigma_m=sd_m^2*(matrix(rho_m,nrow = t, ncol = t)+diag(rep(1-rho_m,t)))
  Sigma_x=sd_x^2*(matrix(rho_x,nrow = t, ncol = t)+diag(rep(1-rho_x,t)))
  Sigma_w=sd_w^2*(matrix(rho_w,nrow = t, ncol = t)+diag(rep(1-rho_w,t)))
  
  zero = rep(0,times=t)
  
  ############### create some matrices to store results #############
  
  #### functional variable #####
  beta_X = matrix(nrow=t,ncol=sim_liter) #results matrix from the quantile regression based on the true X
  beta_simex = matrix(nrow=t,ncol=sim_liter)  #SIMEX estimator when assuming delta(t) is 1
  beta_naive = matrix(nrow=t,ncol=sim_liter)  ## naive estimator
  beta_simex.delta = matrix(nrow=t,ncol=sim_liter)  ## SIMEX estimator when delta(t) is estimated as raw ratio
  beta_simex.delta.known = matrix(nrow=t,ncol=sim_liter)  ## SIMEX estimator when delta(t) is known
  beta_simex.deltas = matrix(nrow=t,ncol=sim_liter)  ## SIMEX estimator when delta(t) is estimated as smoothed ratio
 
  ###### continuous error-free covariate ####
  beta_X_EF = matrix(nrow=1,ncol=sim_liter) ## ## benchmark estimator 
  beta_naive_EF = matrix(nrow=1,ncol=sim_liter) # SIMEX estimator when assuming delta(t) is 1
  beta_W_EF = matrix(nrow=1,ncol=sim_liter)  ### EF estimator when ignore ME in functional varibale 
  beta_naive_EF.delta = matrix(nrow=1,ncol=sim_liter)
  beta_naive_EF.deltas = matrix(nrow=1,ncol=sim_liter)
  beta_naive_EF.delta.known = matrix(nrow=1,ncol=sim_liter)
  ###### binary error-free covariate ####
  beta_X_EF.b = matrix(nrow=1,ncol=sim_liter) ## benchmark estimator 
  beta_naive_EF.b = matrix(nrow=1,ncol=sim_liter) ## SIMEX estimator when assuming delta(t) is 1
  beta_W_EF.b = matrix(nrow=1,ncol=sim_liter)  ### EF estimator when ignore ME in functional variable 
  beta_naive_EF.delta.b = matrix(nrow=1,ncol=sim_liter)
  beta_naive_EF.deltas.b = matrix(nrow=1,ncol=sim_liter)
  beta_naive_EF.delta.known.b = matrix(nrow=1,ncol=sim_liter)
  
  
  #### number of basis functions k #### 
  selected_k_naive=c() ## k based on BIC value of naive model 
  selected_k_bench=c() ## k based on BIC value of benchmark model (with X(t))
  selected_k_n=c() ## k based on sample size 
  
  #######  estimated delta(t) #####
  delta_t.hat = matrix(nrow=t, ncol = sim_liter) ## raw delta(t)
  delta_t.hat.smooth = matrix(nrow=t, ncol = sim_liter) ## smoothed delta(t)
  
  
    
  for(iter in 1:sim_liter){
    
    #-------------------------------------------------------------#
    ################### Generate X, W, M, Y,and EF ################
    #-------------------------------------------------------------#
    
    ######### generate error terms #####
    err_x = mvrnorm(n, zero, Sigma_x) ### error term of X(t)
    err_w = mvrnorm(n, zero, Sigma_w)### error term of W(t)
    err_m = mvrnorm(n, zero, Sigma_m)### error term of M(t)
    err_ef = rnorm(n, 0, sigma_ef) ### error term of continuous error-free covariate
  
    
    ### generate delta(t)
    delta_t =c*Betafunc(a,"f1") +1 # c is the constant might change for different delta(t)
    
    X_t = matrix(rep(Betafunc(a,fx),n),nrow=n,byrow=TRUE) + err_x 
    M_t = t(apply(X_t, 1, function(s) delta_t*s)) + err_m 
    W_t = X_t + err_w
    
    EF = rnorm(n, 0, sigma_ef)  ## EF predictor that is NOT correlated with X(t)
    EF.b = rbinom(n, 1,p) ## binary predictor 
    
     #### generate Y ### 
    gamma.1 = 2 ## true association of Y and EF
    gamma.3 = 0.6 ## true association of Y and biarny EF
    
    Y_norm = crossprod(t(X_t),Betafunc(a,fi))/length(a)+ ## true association of Y and x(t)
      crossprod(t(EF), gamma.1)+ ## true association of Y and EF
      crossprod(t(EF.b), gamma.3)+
      rnorm(n,0,sd=sigma_y) 

    #-------------------------------------------------------------#
    ######################### estimate delta_t   #################
    #-------------------------------------------------------------#
    ### delta(t) is estimated as the ratio of M(t) and W(t) based on the assumption 
    ### that M(t) = delta(t)*W(t)
    delta_t.hat.smooth[,iter] = lowess(a, (colMeans(M_t))/(colMeans(W_t)), f=1/3)$y ## smoothed delta(t)
    delta_t.hat[,iter] = (colMeans(M_t))/(colMeans(W_t)) ## non-smoothed delta(t) 
    
    
    ## create new M_t by dividing the old one with estimated delta_t
    M_t.new = M_t/matrix(rep(delta_t.hat[,iter],n), nrow=n, ncol=t, byrow = T)  
    M_t.news = M_t/matrix(rep(delta_t.hat.smooth[,iter],n), nrow=n, ncol=t, byrow = T) 
    M_t.known = M_t/matrix(rep(delta_t,n), nrow=n, ncol=t, byrow = T) 
    
    #-------------------------------------------------------------#
    ############################# Select k  #######################
    #-------------------------------------------------------------#
    
    BIC_value_naive = sapply(5:k_max, L_BIC, Y=Y_norm, W = W_t, tau = p_sim, EF=EF, time_interval = a)
    k = (5:k_max)[which.min(BIC_value_naive)]
    selected_k_naive[iter]=k
    
    BIC_value_bench = sapply(5:k_max, L_BIC, Y=Y_norm, W = X_t, tau = p_sim,  EF=EF,time_interval = a)
    selected_k_bench[iter]=(5:k_max)[which.min(BIC_value_bench)]
    
    k = ceiling(n^(1/5))+2
    selected_k_n[iter] = k
    bs2 = bs(a, df = k, degree = 3, intercept = TRUE)
    
    #-------------------------------------------------------------#
    #################### basis expansion of M_t,W_t,X_t  ##########
    #-------------------------------------------------------------#
    
    W_i = t(t(W_t)-colMeans(W_t))%*%bs2/length(a)   #(n,k) matrix inner product
    X_i = t(t(X_t)-colMeans(X_t))%*%bs2/length(a)   #(n,k) matrix
    M_i = t(t(M_t)-colMeans(M_t))%*%bs2/length(a)   #(n,k) matrix
    M_i.new = t(t(M_t.new)-colMeans(M_t.new))%*%bs2/length(a)   #(n,k) matrix
    M_i.news = t(t(M_t.news)-colMeans(M_t.news))%*%bs2/length(a)   #(n,k) matrix
    M_i.known = t(t(M_t.known)-colMeans(M_t.known))%*%bs2/length(a)   #(n,k) matrix
     
    #-------------------------------------------------------------#
    #################### Benchmark regression model  ##############
    #-------------------------------------------------------------#
    
    model = rq(Y_norm ~ X_i +EF + EF.b, tau=p_sim)
    quantreg_X = model$coefficients[2:(k+1)]
    beta_X[,iter] = crossprod(t(bs2),quantreg_X)
    beta_X_EF[,iter] = model$coefficients[(k+2)]  ## estimate of EF when functional variable is true measurement
    beta_X_EF.b[,iter] = model$coefficients[(k+3)] 
    
    #-------------------------------------------------------------#
    #################### Naive regression model  ##############
    #-------------------------------------------------------------#
    
    model.W = rq(Y_norm ~ W_i +EF+ EF.b, tau=p_sim)
    quantreg_W = model.W$coefficients[2:(k+1)]
    beta_naive[,iter] = crossprod(t(bs2),quantreg_W)
    beta_W_EF[,iter] = model.W$coefficients[(k+2)] ## estimate of EF when ME in functional variable is ignored
    beta_W_EF.b[,iter] = model.W$coefficients[(k+3)]
    
    #-------------------------------------------------------------#
    ############ Get a reasonable estimate for Sigma_uu  #########
    #-------------------------------------------------------------#
    
    ##### estimated covariance matrix of X(t) ####
    #######by assumption cov(W(t),M(t))/delta= Sigma_xx
    Sigma_xx =try( cov(W_i, M_i),T)     
    Sigma_xx.delta = try( cov(W_i, M_i.new),T)  
    Sigma_xx.deltas = try( cov(W_i, M_i.news),T) 
    Sigma_xx.delta.known = try( cov(W_i, M_i.known),T) 
    
    ### estimated covariance matrix of W(t) ####
    Sigma_ww = var(W_i)  ########covariance matrix of observed surrogate var(W_ic) (centered)
    
    ###### estimate covariance matrix of measurement error U(t) ###
    Sigma_uu.delta = try(Sigma_ww - Sigma_xx.delta,T) ######covariance matrix of W|x from W=X+U
    Sigma_uu.delta = try(make.positive.definite(as.matrix(forceSymmetric(Sigma_uu.delta))) ,T)
    
    Sigma_uu = try(Sigma_ww - Sigma_xx,T)
    Sigma_uu = try(make.positive.definite(as.matrix(forceSymmetric(Sigma_uu))) ,T)
    
    Sigma_uu.deltas = try(Sigma_ww - Sigma_xx.deltas,T)     #######covariance matrix of W|x from W=X+U
    Sigma_uu.deltas = try(make.positive.definite(as.matrix(forceSymmetric(Sigma_uu.deltas))) ,T)
    
    Sigma_uu.delta.known = try(Sigma_ww - Sigma_xx.delta.known,T)     #######covariance matrix of W|x from W=X+U
    Sigma_uu.delta.known = try(make.positive.definite(as.matrix(forceSymmetric(Sigma_uu.delta.known))) ,T)
    
    
    #---------------------------------------------------------------#
    ######################### SIME simulation step  #################
    #---------------------------------------------------------------# 
    
    ##########Simulation step of the SIMEX procedure, see page 101 of Ray's book##############################
    B = 100                        ####number of repliates     
    lambda = seq(0.0001,2.0001,.05)              ###get a set of monotonically increasing small numbers
  
    #### when delta(t) is assumed to be 1 #####
      gamma.simex = lapply(seq(1:B), function(b){
        sapply(lambda, function(s) {
          set.seed(b+iter)
          U_b = try(mvrnorm(n, rep(0, ncol(W_i)), Sigma_uu, empirical = TRUE),T)
          W_lambda =try( W_i + (sqrt(s)*U_b) ,T)
          
          model =  try(rq(Y_norm ~ W_lambda+ EF+ EF.b, tau=p_sim),T)
          try(model$coefficients ,T)
        })
      } )

    #### when delta(t) is estimated as raw ratio #####
      gamma_simex.delta = lapply(seq(1:B), function(b){
        sapply(lambda, function(s) {
          set.seed(b+iter)
          U_b = try(mvrnorm(n, rep(0, ncol(W_i)), Sigma_uu.delta, empirical = TRUE),T)
          W_lambda =try( W_i + (sqrt(s)*U_b) ,T)

          model=  try(rq(Y_norm ~ W_lambda+ EF+ EF.b, tau=p_sim),T)
          try(model$coefficients ,T)
        })
      } )

    #### when delta(t) is estimated as smoothed ratio #####
      gamma_simex.deltas = lapply(seq(1:B), function(b){
        sapply(lambda, function(s) {
          set.seed(b+iter)
          U_b.deltas = try(mvrnorm(n, rep(0, ncol(W_i)), Sigma_uu.deltas, empirical = TRUE),T)
          W_lambda.deltas = try(W_i + (sqrt(s)*U_b.deltas) ,T)

          model =  try(rq(Y_norm ~ W_lambda.deltas+ EF+ EF.b, tau=p_sim),T)
          try(model$coefficients ,T)
        })
      } )

      #### when delta(t) is known #####
      gamma_simex.delta.known = lapply(seq(1:B), function(b){
        sapply(lambda, function(s) {
          set.seed(b+iter)
          U_b.delta.known = try(mvrnorm(n, rep(0, ncol(W_i)), Sigma_uu.delta.known, empirical = TRUE),T)
          W_lambda.delta.known = try(W_i + (sqrt(s)*U_b.delta.known) ,T)
          
          model =  try(rq(Y_norm ~ W_lambda.delta.known+ EF+ EF.b, tau=p_sim),T)
          try(model$coefficients ,T)
        })
      } )
      
    #-------------------------------------------------------------#
    ######################### extrapolation step  #################
    #-------------------------------------------------------------#
    #### get average accross B ###
    gamma_simex.ave = try(Reduce("+", gamma.simex)/B,T)
    gamma_simex.delta.ave = try(Reduce("+", gamma_simex.delta)/B,T)
    gamma_simex.deltas.ave = try(Reduce("+", gamma_simex.deltas)/B,T)
    gamma_simex.delta.known.ave = try(Reduce("+", gamma_simex.delta.known)/B,T)
    
    
    gamma.t.simex = (try(extra_quad((k+3), gamma_simex.ave,lambda,-1),T)) ## all regression coefficients
    gamma.t.simex.delta = (try(extra_quad((k+3), gamma_simex.delta.ave,lambda,-1),T))
    gamma.t.simex.deltas = (try(extra_quad((k+3), gamma_simex.deltas.ave,lambda,-1),T))
    gamma.t.simex.delta.known = (try(extra_quad((k+3), gamma_simex.delta.known.ave,lambda,-1),T))
    
   
    
    #-------------------------------------------------------------#
    ######################### assume delta is 1  #################
    #-------------------------------------------------------------#
    beta_simex[,iter] = try(crossprod(t(bs2),(gamma.t.simex[2:(k+1)])),T)
    beta_naive_EF[,iter]=try((gamma.t.simex[k+2]),T)
    beta_naive_EF.b[,iter]=try((gamma.t.simex[k+3]),T)

    
    #-------------------------------------------------------------#
    ######################### delta estimated  ###############
    #-------------------------------------------------------------#
    
    beta_simex.delta[,iter] = try(crossprod(t(bs2),(gamma.t.simex.delta[2:(k+1)])),T)
    beta_naive_EF.delta[,iter]=try((gamma.t.simex.delta[k+2]),T)
    beta_naive_EF.delta.b[,iter]=try((gamma.t.simex.delta[k+3]),T)

    
    #-------------------------------------------------------------#
    ######################### delta estimated smoothed ############
    #-------------------------------------------------------------#
    
    
    beta_simex.deltas[,iter] = try(crossprod(t(bs2),(gamma.t.simex.deltas[2:(k+1)])),T)
    beta_naive_EF.deltas[,iter]=try((gamma.t.simex.deltas[k+2]),T)
    beta_naive_EF.deltas.b[,iter]=try((gamma.t.simex.deltas[k+3]),T)

    
    #-------------------------------------------------------------#
    ######################### delta known #######################
    #-------------------------------------------------------------#
    beta_simex.delta.known[,iter] = try(crossprod(t(bs2),(gamma.t.simex.delta.known[2:(k+1)])),T)
    beta_naive_EF.delta.known[,iter]=try((gamma.t.simex.delta.known[k+2]),T)
    beta_naive_EF.delta.known.b[,iter]=try((gamma.t.simex.delta.known[k+3]),T)
  
  
  print(iter)
  
  } ### end of simulation iterations
  
  
    #-------------------------------------------------------------#
    ####################### output of the function ################
    #-------------------------------------------------------------#
  
  re=list(beta_X=beta_X,beta_simex=beta_simex,beta_naive=beta_naive,
          selected_k_naive=selected_k_naive,selected_k_bench=selected_k_bench,
          beta_simex.delta = beta_simex.delta, ## SIMEX estimator when delta is assumed to be 1
          beta_simex.delta.known = beta_simex.delta.known, ### SIMEX estimator of functional variable when delta is known
          beta_simex.deltas = beta_simex.deltas, ### SIMEX estimator of functional variable when delta is estimated 
          beta_X_EF=beta_X_EF,beta_W_EF=beta_W_EF,
          beta_X_EF.b=beta_X_EF.b,beta_W_EF.b=beta_W_EF.b,
          beta_naive_EF=beta_naive_EF, #### Estimator of EF covariate when ME in functional variable is corrected by assume delta is 1 
          beta_naive_EF.delta= beta_naive_EF.delta,
          beta_naive_EF.delta.known= beta_naive_EF.delta.known,
          beta_naive_EF.deltas= beta_naive_EF.deltas,
          beta_naive_EF.b=beta_naive_EF.b, #### Estimator of the binary EF covariate when ME in functional variable is corrected by assume delta is 1 
          beta_naive_EF.delta.b= beta_naive_EF.delta.b,
          beta_naive_EF.delta.known.b= beta_naive_EF.delta.known.b,
          beta_naive_EF.deltas.b= beta_naive_EF.deltas.b,
          delta_t= delta_t.hat,
          delta_t.smooth= delta_t.hat.smooth) #### Estimator of EF covariate when ME in functional variable is corrected by estimating delta 
  
  return(re)
  
} ### end of simulation function 

#-------------------------------------------------------------#
########################### Example ###########################
#-------------------------------------------------------------#


### run a 500-replications simulation 
res = main.delta(n=200,c= 0,sd_x=1.5,sd_w=0.75,p_sim=0.5,fi = "f1",fx = "f7",seeds=123,sim_liter=500)




