# A file for a simple version implementation of BART algorithm ---
# with all functions needed to it
# Calculate Tree Prior

# Load in functions which are common to GP-BART
source("common_help_functions.R")

# Function to calculate the tree complete conditional using BART model
tree_complete_conditional_het_bart_mu <- function(x,
                                           residuals_values,
                                           tree,
                                           tau_mu,
                                           precision_vector) {
  
  # Getting the number of observations of the data
  n <- nrow(x)
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]
  
  # Number of nodes
  n_node <- length(terminal_nodes)
  
  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))
  
  # Retrieving Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals_values[x$observations_index]
  })
  
  # Retrieving Residuals terminal nodes
  precision_terminal_nodes <- lapply(terminal_nodes, function(x) {
    precision_vector[x$observations_index]
  })
  
  # Retrieving Residuals terminal nodes
  precision_sum_terminal_nodes <- unlist(lapply(terminal_nodes, function(x) {
    sum(precision_vector[x$observations_index])
  }))
  
  
  # Calculating the  log(p(x_{i}))
  sum_log_precision_terminal_nodes <- unlist(lapply(terminal_nodes, function(x){
    sum(log(precision_vector[x$observations_index]))
  }))
  
  # \sum_p(x_i)*r^2
  sum_precision_residuals_squared <- unlist(mapply(precision_terminal_nodes,
                                                   residuals_terminal_nodes,
                                                   FUN = function(precision,residuals) {
    sum(precision*residuals^2)
  },SIMPLIFY = FALSE))
  
  # Calculating the "last term" \frac{(\sum_{p*ri})^2}{\sum_p + \tau_\mu}
  squared_sum_precision_residuals <- unlist( mapply(precision_sum_terminal_nodes,
                                             precision_terminal_nodes,
                                             residuals_terminal_nodes, FUN = function(precision_sum,
                                                                                      precision,
                                                                                      residuals){
                                        (sum(precision*residuals)^2)/(precision_sum+tau_mu)
                                      }, SIMPLIFY = FALSE))
  
  # Retrieve all nodes values and calculate all of them
  log_posterior <- 0.5*sum_log_precision_terminal_nodes+0.5*log(precision_sum_terminal_nodes+tau_mu)-
    0.5*sum_precision_residuals_squared-0.5*squared_sum_precision_residuals
  
  
  
  return(log_posterior)
}

# Function to calculate the tree complete conditional using het-BART model (FOR \tau)
tree_complete_conditional_het_bart_tau <- function(x,
                                                  precision_sq_residuals_values,
                                                  tree,
                                                  a_tau,
                                                  d_tau) {
  
  # Getting the number of observations of the data
  n <- nrow(x)
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]
  
  # Number of nodes
  n_node <- length(terminal_nodes)
  
  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))
  
  # Retrieving Residuals terminal nodes
  sum_precision_sq_residuals_terminal_nodes <- unlist(lapply(terminal_nodes, function(x) {
    sum(precision_sq_residuals_values[x$observations_index])
  }))
  

  # Retrieve all nodes values and calculate all of them
  log_posterior <- lgamma(0.5*nodes_size+a_tau) - (0.5*nodes_size+a_tau)*log(0.5*sum_precision_sq_residuals_terminal_nodes+d_tau)
  
  return(log_posterior)
}

# Function to calculate the tree complete conditional using het-BART model (FOR \tau)
update_tau_het_bart <- function(x,
                                precision_sq_residuals_values,
                                tree,
                                a_tau,
                                d_tau) {
  
    
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))
  
  # Getting the number of observations of the data
  n <- nrow(x)
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]
  
  # Number of nodes
  n_node <- length(terminal_nodes)
  
  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))
  
  # Retrieving Residuals terminal nodes
  sum_precision_sq_residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    sum(precision_sq_residuals_values[x$observations_index])
  })
  
  # Sampling tau values 
  tau_sample <- mapply(nodes_size, sum_precision_sq_residuals_terminal_nodes,
                       FUN = function(nodes_size, sum_res_sq){
                         rgamma(n = 1,shape = nodes_size*0.5+a_tau,
                                rate = 0.5*sum_res_sq+d_tau)
                       })
  
  # Adding the mu values calculated
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$tau <- tau_sample[[names_terminal_nodes[i]]]
  }
  
  return(tree)
}


# Update \mu using the BART simple version
update_mu_het_bart <- function(tree,
                           x,
                           tau,
                           tau_mu,
                           residuals,
                           precision_vector,
                           seed = NULL) {
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))
  
  # Getting the number of observations of the data
  n <- nrow(x)
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]
  
  # Number of nodes
  n_node <- length(terminal_nodes)
  
  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$observations_index)
  }, numeric(1))
  
  # Retrieving Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$observations_index]
  })
  
  # Retrieving Residuals terminal nodes
  precision_terminal_nodes <- lapply(terminal_nodes, function(x) {
    precision_vector[x$observations_index]
  })
  
  # Retrieving Residuals terminal nodes
  precision_sum_terminal_nodes <- unlist(lapply(terminal_nodes, function(x) {
    sum(precision_vector[x$observations_index])
  }))
  
  
  # Calculating the mean of \mu_post
  mu_post_mean <- mapply(residuals_terminal_nodes,
                         precision_terminal_nodes,
                         precision_sum_terminal_nodes, FUN = function(resid,prec,prec_sum){
                           sum(resid*prec)/(prec_sum+tau_mu)
                         },SIMPLIFY = FALSE)
  
  mu_post_var <- (precision_sum_terminal_nodes+tau_mu)^(-1)
  
  # Calculating mu values
  mu <- mapply(mu_post_mean, mu_post_var,
               FUN = function(x, y) {
                 rnorm(
                   n = 1,
                   mean = x,
                   sd = sqrt(mu_post_var)
                 )
               }, SIMPLIFY = FALSE)
  
  # Adding the mu values calculated
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$mu <- mu[[names_terminal_nodes[i]]]
  }
  return(tree)
}

# Getting the \mu vector from terminal nodes
get_mu_bart <- function(tree, x) {
  
  # New g (new vector prediction for g)
  predictions_new <- rep(NA, length(residuals))
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]
  
  # Getting the \mu_{j} vector
  mu_values <- vapply(terminal_nodes, "[[", numeric(1), "mu")
  
  # Adding the mu values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving g
    predictions_new[terminal_nodes[[i]]$observations_index] <- mu_values[[i]]
  }
  
  return(predictions_new)
}

# Getting the \tau vector from terminal nodes
get_tau_bart <- function(tree, x) {
  
  # New g (new vector prediction for g)
  predictions_new <- rep(NA, length(residuals))
  
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))
  
  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]
  
  # Getting the \tau_{j} vector
  tau_values <- vapply(terminal_nodes, "[[", numeric(1), "tau")
  
  # Adding the tau values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving tau
    predictions_new[terminal_nodes[[i]]$observations_index] <- tau_values[[i]]
  }
  
  return(predictions_new)
}


# Calculating the BART object model
h_bart <- function(x, # Covarariate matrix
                 y, # Target variable
                 number_trees_mu = 2, # Number of trees for \mu
                 number_trees_tau = 2, # Number of trees for \tau
                 control = list(node_min_size = 5,
                                scale_boolean = TRUE),
                 prior_parameter = list(a_tau = 10, # Parameters from the prior
                                        d_tau = 3,
                                        k_bart = 2,
                                        alpha = 0.95,
                                        beta = 2,
                                        prob_tau = 0.9),
                 mcmc_parameter = list(n_iter = 2000, # Parameters from MCMC
                                       burn = 1000,
                                       thin = 1),
                 init = list(tau = 1, # Initial values
                             mu = 0)) {
  
  # Setting the initial values for each one of them, 
  node_min_size <- control[["node_min_size"]]
  scale_boolean <- control[["scale_boolean"]]
  
  a_tau <- prior_parameter[["a_tau"]]
  d_tau <- prior_parameter[["d_tau"]]
  k_bart <- prior_parameter[["k_bart"]]
  alpha <- prior_parameter[["alpha"]]
  beta <- prior_parameter[["beta"]]
  prob_tau <- prior_parameter[["prob_tau"]]
  acc_ratio <- 0
  
  n_iter <- mcmc_parameter[["n_iter"]]
  burn <- mcmc_parameter[["burn"]]
  thin <- mcmc_parameter[["thin"]]
  
  tau <- init$tau
  mu <- init$mu
  
  # Giving names for col of "x" if they don't have it
  if(is.null(colnames(x))){
    colnames(x) <- paste0("x.", seq_len(ncol(x)))
  }
  
  # Create the storage for objects to be exported later
  store_size <- (n_iter - burn)/thin
  tree_mu_store <- list()
  tree_tau_store <- list()
  residuals_store <- list()
  y_hat_store <- matrix(NA, ncol = nrow(x),nrow = store_size)
  precision_hat_store <- matrix(NA, ncol = nrow(x),nrow = store_size)
  
  predictions_store <- list()
  precisions_store <- list()
  
  verb_store_mu_list <- list()
  verb_store_tau_list <- list()
  
  
  # Creating a object to store the urrent trees
  current_trees_mu <- list()
  current_trees_tau <- list()
  
  # Creating the list of initial stumps and predictions
  for(j in seq_len(number_trees_mu)) {
    current_trees_mu[[j]] <- stump_mu(x = x,
                                mu = mu)
  }
  
  for(j in seq_len(number_trees_tau)) {
    current_trees_tau[[j]] <- stump_tau(x = x,
                                      tau = tau)
    
  }
  
  # Naming the trees to get the objects name
  names(current_trees_mu) <- vapply(seq_len(number_trees_mu), function(x) paste0("tree_", x), character(1)) # Naming each tree
  names(current_trees_tau) <- vapply(seq_len(number_trees_tau), function(x) paste0("tree_", x), character(1)) # Naming each tree
  
  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running GP-Sum-Sampler..."
  )
  
  # Getting the maxium and the minimum from the data
  a_min <- min(y)
  b_max <- max(y)
  
  # Scaling the data or not (and changing the hyperparameters with respect to it)
  if(scale_boolean) {
    
    # Scaling the y between -0.5 and 0.5
    y_scale <- normalize_bart(y = y)
    
    # \tau_mu depends over the scale of y
    tau_mu  <- 4 * (k_bart^2) * number_trees_mu/((max(y_scale)-min(y_scale))^2)
    
    # Dividing by multiple trees
    a_tau_original <- a_tau
    a_tau <- a_tau_original^(1/number_trees_tau)
    
    # Getting the optimal tau value
    d_tau_original <- rate_tau(x = x,
                               y = y_scale,
                               prob = prob_tau,
                               shape = a_tau_original)
    d_tau <- d_tau_original^(1/number_trees_tau)
    
  } else {
    
    y_scale <- y
    
    # \tau_mu depends over the scale of y
    tau_mu <- 4 * (k_bart^2) * number_trees_mu/((max(y_scale)-min(y_scale))^2)
    
    # Dividing by multiple trees
    a_tau_original <- a_tau
    a_tau <- a_tau_original^(1/number_trees_tau)
    
    # Getting the optimal tau value
    d_tau_original <- rate_tau(x = x,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau_original)
    d_tau <- d_tau_original^(1/number_trees_tau)
  }
  
  # Creating the init values of predictions
  predictions <- matrix(0,
                        ncol = nrow(x),
                        nrow = number_trees_mu)
  
  precisions <- matrix(1, 
                       ncol = nrow(x),
                       nrow = number_trees_tau)
  
  # Having the partial residuals matrix and precisions
  current_partial_residuals_matrix <- predictions
  current_partial_precision_matrix <- precisions
  
  # Iterating over the MCMC sampling 
  for(i in seq_len(n_iter)) {
    
    utils::setTxtProgressBar(progress_bar, i)
    
    # Saving after the warmup step
    if ((i > burn) && ((i %% thin) == 0)) {
      
      curr <- (i - burn) / thin
      tree_mu_store[[curr]] <- current_trees_mu
      tree_tau_store[[curr]] <- current_trees_tau
      
      y_hat_store[curr,] <- colSums(predictions)
      precision_hat_store[curr,] <- apply(precisions,2,function(x){exp(sum(log(x)))})
      residuals_store[[curr]] <- current_partial_residuals_matrix
      predictions_store[[curr]] <- predictions
      precisions_store[[curr]] <- precisions
      
      verb_store_mu_list[[curr]] <- verb_store_mu
      verb_store_tau_list[[curr]] <- verb_store_tau
      
      
    }
    
    # Verb Store
    verb_store_mu <- data.frame(verb = rep(NA, number_trees_mu),
                             accepted = rep(NA, number_trees_mu),
                             identical = rep(NA, number_trees_mu))
    
    verb_store_tau <- data.frame(verb = rep(NA, number_trees_tau),
                                accepted = rep(NA, number_trees_tau),
                                identical = rep(NA, number_trees_tau))
    
    # Iterating over the \mu trees FIRST
    for(j in seq_len(number_trees_mu)) {
      
      # Making a exception for the case of only one tree
      if(number_trees_mu == 1) {
        
        # Getting the current partial values
        current_partial_residuals <- matrix(y_scale, ncol = length(y_scale))
        
      } else {
        
        # Getting the current partial values
        current_partial_residuals <- y_scale - colSums(predictions[-j,,drop = FALSE])
      }
      
      # Propose a new tree based on the verbs: grow/prune/change/swap
      verb <- sample(c("grow", "prune", "change", "swap"),
                     prob = c(0.25,0.25,0.4,0.1), size = 1)
      
      
      # Case of rotation
      if (i < max(floor(0.1 * burn), 10) || length(current_trees_mu[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
      
      
      # GETTING A NEW TREE
      new_trees_mu <- current_trees_mu # Creating new trees to updated as candidate
      
      new_trees_mu[[j]] <- update_tree_verb(
        tree = current_trees_mu[[j]],
        x = x,
        node_min_size = node_min_size,
        verb = verb
      )
      
      # Calculating the likelihood of the new tree
      likelihood_new <- tree_complete_conditional_het_bart_mu(
        tree = new_trees_mu[[j]], # Calculate the full conditional
        residuals_values = current_partial_residuals,
        x = x,tau_mu = tau_mu,
        precision_vector = apply(precisions,
                                 2,
                                 function(x){exp(sum(log(x)))})
      )
      
      # Calculating the likelihood of the old tree
      likelihood_old <- tree_complete_conditional_het_bart_mu(
        tree = current_trees_mu[[j]], # Calculate the full conditional
        residuals_values = current_partial_residuals,
        x = x,tau_mu = tau_mu,
        precision_vector = apply(precisions,
                                 2,
                                 function(x){exp(sum(log(x)))})
      )
      
      # Extracting only the likelihood
      l_new <- sum(likelihood_new) +
        tree_prior(
          tree = new_trees_mu[[j]], # Calculate the tree prior
          alpha = alpha,
          beta = beta
        )
      
      # Extracting only the likelihood
      l_old <- sum(likelihood_old) +
        tree_prior(
          tree = current_trees_mu[[j]], # Calculate the tree prior
          alpha = alpha,
          beta = beta
        )
      
      # (log) Probability of accept the new proposed tree
      acceptance <- l_new - l_old
      
      # In case of acceptance
      if(acceptance > 0 || acceptance > -rexp(1)) {
        
        # Counting acc ratio
        acc_ratio <- acc_ratio + 1
        
        # Make changes if accept
        current_trees_mu <- new_trees_mu
        
        # Create a data.frame with the verb and if it was accepted or not
        verb_store_mu[j,"verb"] <- verb
        verb_store_mu[j,"accepted"] <- TRUE
        
      } else {
        
        # Create a data.frame with the verb and if it was accepted or not
        verb_store_mu[j,"verb"] <- verb
        verb_store_mu[j,"accepted"] <- FALSE
        
      } # End of accept for MH sample
      
      # # # To update the mu values
      current_trees_mu[[j]] <- update_mu_het_bart(
        tree = current_trees_mu[[j]],
        x = x,
        residuals = current_partial_residuals, 
        tau_mu = tau_mu,
        precision_vector = apply(precisions,
                                 2,
                                 function(x){exp(sum(log(x)))})
        )
      
      # EQUATION FROM SECTION 4
      # ==== Using the prediction from R_star_bar
      predictions[j,] <- get_mu_bart(
        tree = current_trees_mu[[j]], x = x
      )
      current_partial_residuals_matrix[j,] <- current_partial_residuals
      
    } # End of iterations over mu trees
    
    # ========
    # ITERATING NOW OVER THE \TAU TREES 
    # ========
    # Iterating over the trees 
    for(j in seq_len(number_trees_tau)) {
      
      # Making a exception for the case of only one tree
      if(number_trees_tau == 1) {
        
        # Getting the current partial values
        current_partial_residuals_precision_squared <- matrix((y_scale^2)/c(precisions), ncol = length(y_scale))
        # MAYBE NEED TO REVIEW THIS LINE
        
      } else {
        
        # Getting the current partial values
        current_partial_residual_precision_squared <- ((y_scale - colSums(predictions[,,drop = FALSE]))^2)*apply(precisions[-j,,drop = FALSE],2, function(x){exp(sum(log(x)))})
      }
      
      # Propose a new tree based on the verbs: grow/prune/change/swap
      verb <- sample(c("grow", "prune", "change", "swap"),
                     prob = c(0.25,0.25,0.4,0.1), size = 1)
      
      
      # Case of rotation
      if (i < max(floor(0.1 * burn), 10) || length(current_trees_tau[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
      
      
      # GETTING A NEW TREE
      new_trees_tau <- current_trees_tau # Creating new trees to updated as candidate
      
      new_trees_tau[[j]] <- update_tree_verb(
        tree = current_trees_tau[[j]],
        x = x,
        node_min_size = node_min_size,
        verb = verb
      )
      
      # Calculating the likelihood of the new tree
      likelihood_new <- tree_complete_conditional_het_bart_tau(
        tree = new_trees_tau[[j]], # Calculate the full conditional
        precision_sq_residuals_values = current_partial_residual_precision_squared,
        x = x,a_tau = a_tau, d_tau = d_tau
      )
      
      # Calculating the likelihood of the old tree
      likelihood_old <- tree_complete_conditional_het_bart_tau(
        tree = current_trees_tau[[j]], # Calculate the full conditional
        precision_sq_residuals_values = current_partial_residual_precision_squared,
        x = x,a_tau = a_tau, d_tau = d_tau
      )
      
      # Extracting only the likelihood
      l_new <- sum(likelihood_new) +
        tree_prior(
          tree = new_trees_tau[[j]], # Calculate the tree prior
          alpha = alpha,
          beta = beta
        )
      
      # Extracting only the likelihood
      l_old <- sum(likelihood_old) +
        tree_prior(
          tree = current_trees_tau[[j]], # Calculate the tree prior
          alpha = alpha,
          beta = beta
        )
      
      # (log) Probability of accept the new proposed tree
      acceptance <- l_new - l_old
      
      # In case of acceptance
      if(acceptance > 0 || acceptance > -rexp(1)) {
        
        # Counting acc ratio
        acc_ratio <- acc_ratio + 1
        
        # Make changes if accept
        current_trees_tau <- new_trees_tau
        
        # Create a data.frame with the verb and if it was accepted or not
        verb_store_tau[j,"verb"] <- verb
        verb_store_tau[j,"accepted"] <- TRUE
        
      } else {
        
        # Create a data.frame with the verb and if it was accepted or not
        verb_store_tau[j,"verb"] <- verb
        verb_store_tau[j,"accepted"] <- FALSE
        
      } # End of accept for MH sample
      
      # # # To update the mu values
      current_trees_tau[[j]] <- update_tau_het_bart(
        tree = current_trees_tau[[j]],
        x = x,
        precision_sq_residuals_values = current_partial_residual_precision_squared,
        a_tau = a_tau,d_tau = d_tau)

      # EQUATION FROM SECTION 4
      # ==== Using the prediction from R_star_bar
      precisions[j,] <- get_tau_bart(
        tree = current_trees_tau[[j]], x = x
      )
      current_partial_precision_matrix[j,] <- current_partial_residual_precision_squared
      
    }
    
    
  } # End of the iteration over MCMC sampling
  
  results <- list(x = x, # Covariate matrix
                  y = y_scale, # Target variable
                  tree_mu_store = tree_mu_store, # Tree store
                  tree_tau_store  = tree_tau_store,
                  residuals_store = residuals_store,
                  predictions_store = predictions_store,
                  precision_hat_store = precision_hat_store,
                  precisions_store = precisions_store,
                  y_hat_store = y_hat_store,
                  number_trees_mu = number_trees_mu, # Number of trees mu
                  number_trees_tau = number_trees_tau, # Number of trees
                  control = list(node_min_size = node_min_size,
                                 scale_boolean = scale_boolean,
                                 a_min = a_min,
                                 b_max = b_max),
                  prior_parameter = list(a_tau = a_tau, # Parameters from the prior
                                         d_tau = d_tau,
                                         k_bart = k_bart,
                                         alpha = alpha,
                                         beta = beta,
                                         prob_tau = prob_tau),
                  mcmc_parameter = list(n_iter = n_iter, # Parameters from MCMC
                                        burn = burn,
                                        thin = thin),
                  verb_store_mu = verb_store_mu_list,
                  verb_store_tau = verb_store_tau_list)
  class(results) <- "gpbart_BART"
  return(results)
}

# Getting the predictions for each tree
get_predictions_tree <- function(tree,
                                 x) {
  
  new_tree <- tree
  pred_new <- numeric(nrow(x))
  # Setting the root node with the new observations
  new_tree[["node_0"]]$test_index <- seq_len(nrow(x))
  
  # Creating the list of nodes
  list_nodes <- names(new_tree)[-1]
  
  # IN CASE OF JUST ROOT NODE
  if (length(new_tree) == 1) {
    list_nodes <- "node_0"
  }
  
  # Iterating over the list of nodes
  for(i in seq_along(list_nodes)) {
    
    # Selecting the current node
    current_node_aux <- new_tree[[list_nodes[[i]]]]
    
    # To continuous covariates
    if(is.numeric(x[, current_node_aux$node_var])) {
      
      # Updating observations from the left node
      if (current_node_aux$left == 1) {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var]  < current_node_aux$node_var_split)] # Updating the left node
      } else {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] >= current_node_aux$node_var_split)]
      }
      
      # To categorical covariates
    } else {
      # Updating observations from the left node
      if(current_node_aux$left == 1) {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
      } else {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
      }
    } # Ending the else for type of variable
    
    # Checking if it is a terminal node or not
    if(new_tree[[list_nodes[[i]]]]$terminal == 1 && length(new_tree[[list_nodes[[i]]]]$test_index) > 0) {
      
      # Assign the mu value for that terminal node
      pred_new[new_tree[[list_nodes[[i]]]]$test_index] <- new_tree[[list_nodes[[i]]]]$mu
      
    } # End of check terminal
  } # End Terminal node iterations
  
  return(pred_new)
}

# Getting the predictions for each tree
get_tau_tree <- function(tree,
                         x) {
  
  new_tree <- tree
  pred_tau_new <- numeric(nrow(x))
  # Setting the root node with the new observations
  new_tree[["node_0"]]$test_index <- seq_len(nrow(x))
  
  # Creating the list of nodes
  list_nodes <- names(new_tree)[-1]
  
  # IN CASE OF JUST ROOT NODE
  if (length(new_tree) == 1) {
    list_nodes <- "node_0"
  }
  
  # Iterating over the list of nodes
  for(i in seq_along(list_nodes)) {
    
    # Selecting the current node
    current_node_aux <- new_tree[[list_nodes[[i]]]]
    
    # To continuous covariates
    if(is.numeric(x[, current_node_aux$node_var])) {
      
      # Updating observations from the left node
      if (current_node_aux$left == 1) {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var]  < current_node_aux$node_var_split)] # Updating the left node
      } else {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] >= current_node_aux$node_var_split)]
      }
      
      # To categorical covariates
    } else {
      # Updating observations from the left node
      if(current_node_aux$left == 1) {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] == current_node_aux$node_var_split)] # Updating the left node
      } else {
        new_tree[[list_nodes[i]]]$test_index <- new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index[which(x[new_tree[[paste0("node_", current_node_aux$parent_node)]]$test_index, current_node_aux$node_var] != current_node_aux$node_var_split)]
      }
    } # Ending the else for type of variable
    
    # Checking if it is a terminal node or not
    if(new_tree[[list_nodes[[i]]]]$terminal == 1 && length(new_tree[[list_nodes[[i]]]]$test_index) > 0) {
      
      # Assign the mu value for that terminal node
      pred_tau_new[new_tree[[list_nodes[[i]]]]$test_index] <- new_tree[[list_nodes[[i]]]]$tau
      
    } # End of check terminal
  } # End Terminal node iterations
  
  return(pred_tau_new)
}

# Doing a productory for multiple trees we would have
get_prediction_prod <- function(multiple_trees,x){
  
  # Creating a matrix to store all the trees
  n <- nrow(x)
  n_trees <- length(multiple_trees)
  tau_trees_store <- matrix(NA, nrow = n_trees, ncol = n)
  
  # Iterating over multiple trees
  for(i in 1:n_trees){
      # Storing each
      tau_trees_store[i,] <- get_tau_tree(tree = multiple_trees[[i]], x = x)
  }
  # Doing the productory of columns
  return(apply(tau_trees_store, 2, function(x) { exp(sum(log(x)))})) # Just to avoid overflows
}

# Predicting for new observation
predict.gpbart_BART <- function(bart_mod, newdata, type = c("all")) {
  
  # Getting the MCMC prediction values
  mcmc_post_pred <- matrix(NA,
                           nrow = length(bart_mod$tree_store),
                           ncol = nrow(newdata))
  
  # Setting up the colnames
  colnames(newdata) <- colnames(bart_mod$x)
  
  # Retrieving all information necessary to predict for new observations
  for(i in seq_along(bart_mod$tree_store)) {
    
    current_trees <- bart_mod$tree_store[[i]]
    
    # Creating a pred_aux matrix to store the predictions for new obs
    pred_aux <- matrix(0,
                       nrow = bart_mod$number_trees_mu,
                       ncol = nrow(newdata))
    
    for(j in seq_len(bart_mod$number_trees_mu)) {
      pred_aux[j,] <- get_predictions_tree(tree = current_trees[[j]],x = newdata)
    }
    
    # Adding up all trees
    mcmc_post_pred[i,] <- colSums(pred_aux)
  }
  
  # Getting the values
  out <- switch(type,
                all = unnormalize_bart(z = mcmc_post_pred, a = bart_mod$control$a_min, b = bart_mod$control$b_max),
                mean = unnormalize_bart(z = colMeans(mcmc_post_pred), a = bart_mod$control$a_min, b = bart_mod$control$b_max),
                median = unnormalize_bart(z = colMeans(mcmc_post_pred), a = bart_mod$control$a_min, b = bart_mod$control$b_max))
  return(out)
}