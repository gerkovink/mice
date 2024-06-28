#' Description of xgb function.
#' 
#' This function is called from \code{mice} function. Function uses XGBoost algorithm to impute missing data.There are three different match.types that can be used to obtain final predictions.
#' match.type = "predicted" uses actual predictions from XGBoost algorithm 
#' while match.type = "predicted.observed" or "original.observed" combines XGBoost predictions with donor based selection. 
#' For "predicted.observed" match.type, missing data prediction of each data point (returned from XGBoost algorithm) are matched with predictions obtained form training data (observed) data to select 5 donors. 
#' One donor is randomly picked as an imputed value.
#' For match.type = "original.obseved", donors are selected from observed data.
#' #' XGBoost uses a set of hyperparameters for tuning the algorithm. User has three options to select these parameters
#' 1) Using default parameters (adapted from Deng et al., 2023).
#' 2) Run xgb_param_calc function to let function randomly select one variable as a target and tune the parameters using only one target variable.
#' 3) Run xgb_param_calc function to iteratively select best set of hyperparameters for all variables.

#'
#'@aliases mice.impute.xgb
#' @param y: Vector of values of response variable
#' @param ry: Vector of true/false representing observed values in target variable. TRUE representing respective value was observed and FALSE represents missing value.
#' @param wy: Vector of true/false representing missing values in target variable. TRUE representing respective value was missing and FALSE represents observed value.
#' @param x: dataframe. Data of all variables except target variable.
#' @param xgb.params: list of XGBoost parameters obtained from xgb_param_calc function. User can also manually provide a list of parameters.
#' @param nthreads: Number of threads. By default, function uses all available - 1 cores.
#' @param match.type: Possible value can be "predicted.observed", "original.observed" or "predicted"
#' @return Imputed values for target values,same type as \code{y}, and of length  \code{sum(wy)}; and training error (actual - predicted)
#' 
#'mice.impute.xgb
#' @examples
#' Example 1
#' mice(nhanes, method = "xgb", param = NULL)
#' 
#' Example 2: Running mice.impute.xgb with xgb_param_calc
#' sigma <- matrix(data = c(1, 0.8, 0.8, 1), ncol = 2) # covariance matrix
#' expr = mvtnorm::rmvnorm(n = 1000, mean = c(4, 1), sigma = sigma) %>% 
#'                                    as_tibble() %>% # make into a tibble
#'                                     rename(x = V1, z = V2) %>% # rename columns
#'                                     mutate(y = as.numeric(x +  z  + rnorm(1000, sd = 2)))
#'  expr_mis <- ampute(expr, prop = 0.5, mech = "MAR", type = "RIGHT")
#' 
#' # Running algorithm for sequential parameter tuning for all variables x,y,z
#' 
#'  fixed_param <- list(x = NULL, y = NULL, z = NULL)
#'  fixed_param[["x"]]$parameter$subsample <-0.7
#'  fixed_param[["y"]]$parameter$subsample <-0.5
#'  fixed_param[["z"]]$parameter$subsample <-0.4
#' 
#' res1 <- xgb_param_calc(expr_mis,response = "all", select_features=NULL, iter = 50, fixed_params =fixed_param, Params_to_tune = c("eta", "alpha"), plots = FALSE)
#' mice_res <- mice(expr_mis,  method = "xgb", xgb.params = res1,match.type = "predicted.observed")
#' 
#' Example 3
#' 
#' fixed_param_2 <- list(subsample = 0.7)
#' res2 <- xgb_param_calc(expr_mis,response = "all", select_features=NULL, iter = 50, fixed_params =fixed_param_2, Params_to_tune = c("eta", "alpha"), plots = FALSE)
#' mice_res <- mice(expr_mis,  method = "xgb", xgb.params = res2, match.type = "predicted.observed")

#' @export


mice.impute.xgb<- function(y, ry, x, wy=NULL,xgb.params = NULL, nthreads = -1,match.type = "predicted", k = 5, early_stopping_rounds = NULL, ...) {
  
  ###########################################################################################################
  ############################          Parameters
  ###########################################################################################################
  
  # hardcoding general parameters
  nrounds = 100 # number of trees. This will be overridden if user has provided any number of trees in xgb.params
  nthread <- nthreads
  print_every_n = 50 # to make console less crowded but to also check if algorithm is also working.
  verbose = 0
  
 # Default set of parameters are adapted from Deng etal, 2023
  
  def_params <- list(tree_method = "hist", eta = 0.3, gamma = 0, max_depth = 3, 
                 min_child_weight = 1, subsample = 1, sampling_method = "uniform", 
                 colsample_bytree = 1, colsample_bylevel = 1, colsample_bynode = 1, 
                 lambda = 1, alpha = 0,max_bin = 256, num_parallel_tree = 1, nthread = -1) 
  
 
   #if no hyperparameters are provided, algorithm uses default set of parameters adapted from Deng etal, 2023
  # If only one set of hyperparameters are provided it is merged with def_params to fill in missing parameters
  # if separate set of hyperparameters are provided for each variable, we first extract parameters for the variable which is being imputed
  
  
  if(is.null(xgb.params)){
    
    params <- def_params
    
  } else if(all(colnames(x) %in% names(xgb.params))) {
    
    param_oneVariable <- xgb.params[[setdiff(names(xgb.params), colnames(x))]]$parameter
    params <- c(param_oneVariable, list(tree_method = "hist", num_parallel_tree = 1))
    if(!is.null(params$nrounds)){
    nrounds <- params$nrounds
    params <- subset(params, !names(params) %in% "nrounds")
    }
    
    
  } else {
    
    if (!is.list(xgb.params)) {
      stop("Parameters should be in a list object.")
    }
    def_params <- list(tree_method = "hist", eta = 0.3, gamma = 0, max_depth = 3, 
                      min_child_weight = 1, subsample = 1, sampling_method = "uniform", 
                      colsample_bytree = 1, colsample_bylevel = 1, colsample_bynode = 1, 
                      lambda = 1, alpha = 0,max_bin = 256, num_parallel_tree = 1, nthread = -1)
    
    params <- modifyList(def_params,xgb.params)
    
    if(!is.null(params$nrounds)){
      nrounds <- params$nrounds
      params <- subset(params, !names(params) %in% "nrounds")
    }
  }
  
  
  ############################################################################################################
  #########             Objective function and evaluation fuction according to data type       ###############
  ############################################################################################################
  
  
  var.type <- switch(class(y)[1],
                     numeric = "numeric",
                     integer = "numeric",
                     factor = ifelse(nlevels(y) == 2, "binary", "multiclass"),
                     ordered = ifelse(nlevels(y) == 2, "binary", "multiclass"),
                     stop(paste0("Unsupported data type: ", class(y)))
  )
  
  obj.type <- switch(var.type,
                     numeric = "reg:squarederror",
                     binary =  "binary:logistic",
                     multiclass = ifelse(match.type=="predicted", "multi:softmax", "multi:softprob"),
                     stop(paste0("Unsupported objective function: ", var.type))
  )
  
  eval.metric <- switch(var.type,
                        numeric = "MAPE",
                        binary = "logloss",
                        multiclass = "mlogloss",
                        stop(paste0("Unsupported data type: ", obj.type))
  )
  
  
  ############################################################################################################
  ####################              Defining matrices for observed and missing data.     #####################
  ############################################################################################################
  
  if (is.null(wy)) 
    wy <- !ry
  
  nmis <- sum(wy)                 # number of missing entries for target variable
  xobs <- x[ry, , drop = FALSE]   # All variable data (-target variable), where target variable was observed
  xmis <- x[wy, , drop = FALSE]   # All variable data (-target variable), where target variable was missing
  yobs <- y[ry]                   # Observed values of target variable
  
  
  if(var.type != "numeric"){
    yobs <- as.integer(yobs) - 1
    bin.t <- sort(table(yobs))
    
    
    if (is.na(bin.t[2])) {
      imp <- levels(yobs)[as.integer(names(bin.t[1])) + 1]
      msg <- paste("The variable", var, "in the data only have single class. Imputation models can't be built.")
      stop(msg)
    }
  }
  
  ###############################################################################################################
  ########################                    XGBoost Training                  #################################
  ###############################################################################################################
  
  dobs <- xgboost::xgb.DMatrix(data = xobs, label = yobs, nthread = nthreads)   #     Defining xgboost object
  dmis <- xgboost::xgb.DMatrix(data = xmis, nthread = nthreads)                 #     Defining xgboost object
  
  watchlist <- list(train = dobs)
  
  if(var.type!="multiclass"){
    xgb.fit <- xgboost::xgb.train(
      data = dobs, objective = obj.type, watchlist = watchlist,
      params = params, nrounds = nrounds, early_stopping_rounds = early_stopping_rounds, 
      print_every_n = print_every_n, verbose = verbose
    ) # training the model on observed data using parameters set earlier on this page
  } else if(var.type=="multiclass") {   #     xgboost for multi-class categorical data also needs number of classes to be given to the training process
    
    N.class <- length(levels(y))
    xgb.fit <- xgboost::xgb.train(
      data = dobs,num_class = N.class, objective = obj.type, watchlist = watchlist,
      params = params, nrounds = nrounds, early_stopping_rounds = early_stopping_rounds, 
      print_every_n = print_every_n, verbose = verbose
    ) # training the model on observed data using parameters set earlier on this page
    # For simple prediction, training will use multi:softmax objective function, but for other two options, multi:softprob will be used which returns raw probabilities,
  }
  
  
  ###############################################################################################################
  ########################                  XGBoost Prediction                  #################################
  ###############################################################################################################
  
  yhatmis <- predict(xgb.fit, dmis) #obtaining the prediction for missing values in target y variable
  yhatobs <- predict(xgb.fit, dobs) # obtaining the prediction for observed values in target y variable
  trainerror_obs <- yhatobs # To calculate training error at the end of the script
  
  if(var.type == "numeric"){
    
    if(match.type == "predicted" | is.null(match.type)) {
      imp <- yhatmis # will simply assign predictions to imputation output
    } else if(match.type == "predicted.observed"){
      idx <- matchindex(d = yhatobs, t = yhatmis, k = k) # will use predictions of missing data to match with predictions of observed data to find the index of closest match.     h
      imp <-yobs[idx] # retrieve observed values from indexes retrieved in idx adn assign these values as imputation output
    } else if(match.type == "original.observed"){
      idx <- matchindex(d = yobs, t = yhatmis, k = k) # will use predictions of missing data to match with original observed data to find the index of closest match.
      imp <-yobs[idx]   # retrieve observed values from these indices as imputation output
    }
    
  }else if(var.type == "binary"){
    
    if(match.type == "predicted"| is.null(match.type)){
      yhatmis <- ifelse(yhatmis >= 0.5, 1, 0) # Prediction output is continuous. To make it binary again,if prediction probability value is greater than 0.5, it will be converted to 1, else 0.
      yhatmis <- levels(y)[yhatmis + 1] #since values were changes to integer-1 to convert them to 0 and 1. Now, 1 is added to match prediction values with labels of target variable and convert them to corresponding labels.
      imp <- yhatmis
    } else if(match.type == "predicted.observed"){
      idx <- matchindex(d = yhatobs, t = yhatmis, k = k) # index of closest match of predicted missing value with predictions from observed data
      imp <- levels(y)[yobs[idx]+1] # retrieve corresponding observed values from these indices
    } else if(match.type == "original.observed"){
      yhatobs <- yobs #convert original observed value to numeric data. See quantify function from mice.impute.pmm
      idx <- matchindex(d = yhatobs, t = yhatmis, k = k) 
      imp <-levels(y)[yobs[idx]+1]
    }
    
  } else if(var.type == "multiclass"){
    if(match.type == "predicted"| is.null(match.type)){
      yhatmis <- levels(y)[yhatmis + 1] #convert numeric values to corresponding categories found in target variable originally
      imp <- yhatmis
    } else {
      if(match.type == "predicted.observed"){
        yhatobs <- predict(xgb.fit, dobs, reshape = TRUE) #obtain predictions for observed data. Reshape=TRUE is added so that each column represent probability of each category
      } else if(match.type == "original.observed"){
        yhatobs <- as.matrix(sparse.model.matrix(~y-1, data = as.data.frame(y)))[ry,]
      }
      yhatmis <- predict(xgb.fit, dmis, reshape = TRUE) #obtain predictions for missing entries. Reshape=TRUE is added so that each column represent probability of each category
      match.class <- Rfast::knn(xnew = yhatmis, y = y, x = yhatobs, k = k) # use knn to find closest match.
      imp <- levels(y)[match.class]
    }
  }
  ############################ Training error #############################################
  squared_errors <- yobs - trainerror_obs
  error <- mean(squared_errors)

  
  ############################ Training error #############################################
  
  return(list(imp, error))
  
}



