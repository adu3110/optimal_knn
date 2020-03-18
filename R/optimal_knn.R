#' Optimal KNN
#'
#' Optimal kNN function
#' @param formula : formula for the model
#' @param data : data frame containing input data
#' @param method : whether poisson or binomial distribution is to be used to determine optimal k
#' @param k_list: list of k values for grid search
#' @param num_repeats: Number of times k-NN algorthim should be bootstrapped for each k
#' @param bootstrap_fraction: fraction of data to be bottstrapped for each iteration
#' @return optimal_k: data frame containing optimal k value and other stats
#' @export

optimal_knn <- function(formula, data, method = c("poisson", "binomial"), k_list = c(3:9), 
                        num_repeats = 5, bootstrap_fraction=0.7){
  
  optimal_k <- list()
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  
  mf <- eval(expr = mf, envir = parent.frame())
  
  if(!(all(apply(mf[,2:ncol(mf)], 2, is.numeric)))){
    stop("Only numeric predictor variables are allowed")
  }
  
  X <- as.matrix(mf[,2:ncol(mf)])
  
  y <- as.matrix(mf[, 1])
  y_factors <- as.factor(y)
  
  y_levels <- levels(y_factors)
  
  if(length(y_levels) != 2){
    stop("Only two unique values allowed for y")
  }
  
  if(length(k_list) < 2){
    stop("provide at least 2 values for k")
  }
  
  if(num_repeats < 1){
    stop("num_repeats should be at least 1")
  }
  
  if(bootstrap_fraction <= 0 || bootstrap_fraction >= 1){
    stop("bootstrap_fraction should be at between 0 and 1")
  }
  
  num_obs <- nrow(x)
  M <- sum(y == y_levels[1])
  N <- num_obs - M
  k_error <- data.frame(k_values = k_list,
                        mean_error = rep(0, length(k_list)))
  for(k in k_list){
    error_fraction <- rep(0, num_repeats)
    for(num_repeat in num_repeats){
      M_star <- M + 1
      while(M_star > M){
        if(method == "poisson"){
          M_star <- rpois(1, M)
        }else if(method == "binomial"){
          M_star <- rbinom(1, M+N, M/(M+N))
        }
      }
      
      N_star <- N + 1
      while(N_star > N){
        if(method == "poisson"){
          N_star <- rpois(1, N)
        }else if(method == "binomial"){
          N_star <- rbinom(1, M+N, N/(M+N))
        }
      }
      
      M1_star <- as.integer(bootstrap_fraction * M_star)
      N1_star <- as.integer(bootstrap_fraction * N_star)
      
      if(M1_star <= k || N1_star <= k){
        stop("not enough samples")
      }
      
      x1_star = sample(x[y == y_levels[1], ], M_star)
      x2_star = sample(x[y == y_levels[2], ], N_star)
      
      x1_indices = sample(M_star, M1_star)
      x2_indices = sample(M_star, M1_star)
      
      x_train <- rbind(x1_star[x1_indices, ], x2_star[x2_indices, ])
      x_test <- rbind(x1_star[-x1_indices, ], x2_star[-x2_indices, ])
      
      y_train <- rbind(matrix(y_levels[1], nrow = M1_star, ncol = 1),
                       matrix(y_levels[2], nrow = N1_star, ncol = 1))
      y_test <- rbind(matrix(y_levels[1], nrow = M_star - M1_star, ncol = 1),
                      matrix(y_levels[2], nrow = N_star - N1_star, ncol = 1))
      distance_matrix <- apply(x_test, 1, 
                             function(test_row){
                               apply(sweep(x_train, 2, test_row, "-")^2, 
                                     1, function(train_row){sqrt(sum(train_row))})})
      distance_matrix <- apply(distance_matrix, 2, 
                               function(dist_col){sort(dist_col, index.return = TRUE)$ix})
      y_test_matrix <- apply(distance_matrix, 2, function(dist_col){y_test[dist_col<=k]})
      y_test_pred <- apply(y_test_matrix, 2, 
                           function(y_mat_col){ifelse(sum(y_mat_col == y_levels[1]) >= k/2, 
                                                      y_levels[1], y_levels[2])})
      error_fraction[num_repeat] <- sum(y_test != y_test_pred)/nrow(y_test)
    }
    k_error[k_error$k_values == k, ]$mean_error <- mean(error_fraction)
  }
  
  optimal_k$k_errors <- k_error
  optimal_k$k_opt_bootstrap <- k_error[which.min(k_error$mean_error), ]$k_values
  optimal_k$bootstrap_fraction <- bootstrap_fraction
  optimal_k$num_predictors <- ncol(x)
  optimal_k$k_optimal <- optimal_k$k_opt_bootstrap * (bootstrap_fraction)^((-4)/(4+ncol(x)))
  
  return(optimal_k)
  
}


