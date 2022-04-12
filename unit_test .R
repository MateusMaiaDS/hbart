# This file is just to measure if the main code for BART is running well, isn't made for 
# evaluation purposes.

# Importing data and library
rm(list=ls())
library(purrr)
library(mlbench)
source("het_bart.R")
source("tree_manipulation_objects.R")
source("common_help_functions.R")

sim_data_het <- function(n, seed = NULL){
  
  # Setting the seed
  set.seed(seed)  
  
  # Creating the variables
  x <- sort(runif(n = n))
  y <- (4*x^2+0.2*exp(2*x)*rnorm(n = n))
  # y <- sin(x) + rnorm(n = n,sd = 0.25)
  return(list(x = x, y = y))
  
}

# Getting a simple train and test dataset
train_data <- sim_data_het(n = 200,seed = 42)
x_train <- train_data$x %>% as.matrix()
colnames(x_train) <- "x"
y_train <- train_data$y %>% as.matrix() #%>% normalize_bart() %>% as.matrix()
colnames(y_train) <- "y"

# test_data <- mlbench.friedman1(n = 100, sd = 1)
# x_test <- test_data$x
# y_test <- test_data$y

bart_mod <- h_bart(x = x_train,y = y_train,
                 control = list(node_min_size = 5,
                                scale_boolean = TRUE),
                 prior_parameter = list(a_tau = 10, # Parameters from the prior
                                       d_tau = 3,
                                       k_bart = 2,
                                       alpha = 0.95,
                                       beta = 2,
                                       prob_tau = 0.9),
                 number_trees_mu = 5,
                 number_trees_tau = 5,
                 init = list(tau = 1, # Initial values
                             mu = 0))

# bart_pred <- predict(bart_mod = bart_mod, newdata = x_test, type = "mean")

y_train_pred <- apply(bart_mod$y_hat_store,2,mean)
y_train_pred_sd <- (apply(bart_mod$precision_hat_store, 2,mean))^((-1/2))

plot(x_train,y_train %>% normalize_bart() , pch= 20)
lines(x_train,y_train_pred, pch= 20, col = "orange")
lines(x_train,y_train_pred+1.96*y_train_pred_sd, lty = "dashed")
lines(x_train,y_train_pred-1.96*y_train_pred_sd, lty = "dashed")
# lines(x_train,(sin(x_train)-min(train_data$y))/(max(train_data$y)-min(train_data$y))-0.5, col = "blue")
# lines(x_train,sin(x_train), col = "blue")

# Calling my bart function
# source("bart/bart.R")
# source("bart/common_help_functions.R")
# source("bart/tree_manipulation_objects.R")
# 
# my_bart_mod <- bart(x = x_train,y = y_train,number_trees = 5)
# 
# plot(x_train,y_train %>% normalize_bart(),pch=20)
# lines(x_train,apply(my_bart_mod$y_hat_store,2,mean))
# lines(x_train,apply(my_bart_mod$y_hat_store,2,function(x){quantile(x,probs = c(0.025))}),lty = "dashed")
# lines(x_train,apply(my_bart_mod$y_hat_store,2,function(x){quantile(x,probs = c(0.975))}),lty = "dashed")

# CALLING THE HBART LIBRARY
library(rbart)
hbart_mod <- rbart(x.train = x_train,y.train = y_train,x.test = x_train,
                   ntree = 5,ntreeh = 5,k = 2)

dbart_pred <- predict(hbart_mod, newdata = x_train) 

plot(x_train,y_train,pch=20)
# ci <-dbart_pred$mdraws %>% apply(2,function(x){quantile(x, probs = c(0.025,0.975))})
# lines(x_train,ci[1,],col="blue", lty = "dashed")
# lines(x_train,ci[2,],col="blue", lty = "dashed")

lines(x_train,dbart_pred$mdraws %>% colMeans(),col="blue")
lines(x_train,(dbart_pred$mdraws %>% colMeans()) + 1.96*(dbart_pred$sdraws %>% colMeans()),
      col="blue",lty = "dashed")
lines(x_train,(dbart_pred$mdraws %>% colMeans()) - 1.96*(dbart_pred$sdraws %>% colMeans()),
      col="blue",lty = "dashed")

dbart_pred$smean
# bart_pred
# y_test
# plot(bart_pred,dbart_pred)
# abline(a = 0 , b =1)
# (unlist(bart_mod$tau_store)/(max(y_train)-min(y_train))^2)^(-1/2) %>% plot(type = "l")
# dbart_mod$sigma %>% plot(type = "l")
