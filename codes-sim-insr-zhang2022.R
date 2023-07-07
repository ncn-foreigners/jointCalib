## international statistical review
## code based on
## Zhang, Shixiao, Peisong Han, and Changbao Wu. ‘Calibration Techniques Encompassing Survey Sampling, Missing Data Analysis and Causal Inference’. International Statistical Review, 11 August 2022, insr.12518. https://doi.org/10.1111/insr.12518.

### nonprobability sample
library(laeken) ## weightedQuantile
library(data.table)
source("functions.R")
set.seed(2023-07-03)
n <- 3000 # sample size inc
gamma <- c(1,-0.5,0.25,0.1)         
beta1 <- c(210,27.4,13.7,13.7,13.7) 
X1 <- rnorm(n)
X2 <- rnorm(n) 
X3 <- rnorm(n)
X4 <- rnorm(n)
mu <- beta1[1]+beta1[2]*X1+beta1[3]*X2+beta1[4]*X3+beta1[5]*X4
Z1 <- exp(X1/2)
Z2 <- X2/(1+exp(X1))+10
Z3 <- (X1*X3/25+0.6)^3
Z4 <- (X2+X4+20)^2
ps <- 1/(1+exp(gamma[1]*X1+gamma[2]*X2+gamma[3]*X3+gamma[4]*X4)) 
Y <- rnorm(n,mean=mu)
Y_q_true <- quantile(Y, seq(0.1, 0.9,0.1))


## Non-probability sample
results <- list()
R <- 5000
for (r in 1:R) {
    set.seed(r)
    if (r %% 500 == 0) print(r)
    Trt <- rbinom(n,size=1,prob=ps)     # use this for setting 1
    n1 <- sum(Trt)
    n0 <- n - n1
    pop <- data.frame(Y,X1,X2,X3,X4,Z1,Z2,Z3,Z4, Trt)
    sample <- subset(pop, Trt == 1)
    sample$w_b <- n/nrow(sample)
    
    w_res <- calib_quantiles(X_q = cbind(sample$X1,sample$X2,sample$X3,sample$X4),  
                             X = cbind(sample$X1, sample$X2,sample$X3,sample$X4), 
                             d = sample$w_b,  
                             N=n,  
                             totals = c(sum(pop$X1), sum(pop$X2), sum(pop$X3), sum(pop$X2)), 
                             totals_q = list(quantile(pop$X1, seq(0.1, 0.9, 0.1)),
                                             quantile(pop$X2, seq(0.1, 0.9, 0.1)),
                                             quantile(pop$X3, seq(0.1, 0.9, 0.1)),
                                             quantile(pop$X4, seq(0.1, 0.9, 0.1))), 
                             method = "raking",
                             backend = "sampling")
    
    # w_res_z <- calib_quantiles(X_q = cbind(sample$Z1,sample$Z2,sample$Z3,sample$Z4),  
    #                          X = cbind(sample$Z1, sample$Z2,sample$Z3,sample$Z4), 
    #                          d = sample$w_b,  
    #                          N=n,  
    #                          totals = c(sum(pop$Z1), sum(pop$Z2), sum(pop$Z3), sum(pop$Z4)), 
    #                          totals_q = list(quantile(pop$Z1, seq(0.1, 0.9, 0.1)),
    #                                          quantile(pop$Z2, seq(0.1, 0.9, 0.1)),
    #                                          quantile(pop$Z3, seq(0.1, 0.9, 0.1)),
    #                                          quantile(pop$Z4, seq(0.1, 0.9, 0.1))), 
    #                          method = "raking",
    #                          backend = "sampling")
    
    sample$w_cal <- w_res$w
    #sample$w_cal_z <- w_res_z$w
    
    ## standard calibration
    
    w_res_cal <- calibWeights(X = cbind(sample$X1, sample$X2,sample$X3,sample$X4),
                              d = sample$w_b,
                              totals = c(sum(pop$X1), sum(pop$X2), sum(pop$X3), sum(pop$X2)),
                              method = "raking")
    
    sample$w_cal_st <- sample$w_b*w_res_cal
    
    res <- cbind(
      with(sample, weightedQuantile(Y, w_cal, seq(0.1,0.9,0.1))),
      #with(sample, weightedQuantile(Y, w_cal_z, seq(0.1,0.9,0.1))),
      with(sample, weightedQuantile(Y, w_cal_st, seq(0.1,0.9,0.1))),
      seq(0.1,0.9,0.1)
    )
  
    results[[r]] <- as.data.frame(res)
}

results_df <- rbindlist(results, idcol = "r")

results_df[, .(q=mean(V1), s=mean(V2), q_var = var(V1), s_var = var(V2)), V3][
  , ":="(bias_q=abs(q-Y_q_true)/Y_q_true*100,
         bias_s=abs(s-Y_q_true)/Y_q_true*100,
         mse_q = q_var + (q-Y_q_true)^2,
         mse_s = s_var + (s-Y_q_true)^2)][
           , ":="(rmse_q=sqrt(mse_q), rmse_s=sqrt(mse_s),
                  rel = sqrt(mse_s)/sqrt(mse_q))][]

## causal inference




