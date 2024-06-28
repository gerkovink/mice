#' Finds optimal parameters for XGBoost based imputation using Bayesian optimization
#'
#' This is a helper function and can be run before running \code{mice} function with XGBoost based imputation. 
#' The function accepts data in the form of dataframe and finds the best hyperparameters using bayesian optimization. 
#' xgb_pram_calc facilitate tuning of 9 most widely used hyperparameters. 
#' 
#' @aliases xgb_param_calc
#' @param data: Data in the form of matrix or dataframe with each column representing one feature
#' @param response: The default is \code{NULL}. Three options to input the target variable:
#' \itemize{
#'   \item \code{NULL}: The algorithm will randomly pick one target variable from the data.
#'   \item \code{all}: The algorithm will run sequentially, using each feature as the target variable and tuning separate values of parameters for each feature.
#'   \item User-specified index: The user can provide the index of the target variable (e.g., \code{1}), and the algorithm will return parameters tuned for that specific target variable.
#' }
#' @param select_features: The default is \code{NULL}. Three options to input training variable:
#' \itemize{
#'   \item \code{NULL}: same as "all". The algorithm will select all variable other than response variable to use as training feature.
#'   \item \code{all}: same as NULL. The algorithm will select all variable other than response variable to use as training feature.
#'   \item User-specified index: The user can provide the index of the target variable (e.g., \code{c(2,3,4)}), and the algorithm will only use these data of these variable for training of xgboost algorithm and cross-validation.
#' }
#' @param num_cores: The default is \code{1}. User can opt to run xgb_param_calc function with multiple cores or can also use \code{future} package to run bayesian optimization with multiple threads.
#' @param iter: The default is \code{50}. Number of iteration for Bayesian optimization.
#' @param fixed_params: The default is \code{NULL}. If a user wants to fix certain parameters to values other than the default and keep these values fixed while tuning other parameters, they can input those parameters along with their fixed values using the \code{fixed_param} parameter
#' \itemize{
#'   \item \code{response = NULL or index of a target variable}: In case of single target variable, provide a list of parameters with their fixed value (see example for details).
#'   \item \code{response = all}: In case of sequential tuning of response variable, provide a list of fixed parameter for each target variable (see example for details).
#' }
#' @param Params_to_tune: List of parameters to tune. Available parameters that can be tuned include eta, nrounds, max_depth, min_child_weight, subsample, colsample_bytree, alpha, lambda, gamma.
#' @param Plots: {TRUE/FALSE}; if true, returns hyperparameter tuning graphs, otherwise returns NA.
#' @return List with best performing parameters, plots for parameter optimization and response variable that was used to extract those parameters
#' @seealso \code{\link{mice}}
#' @examples
#' 
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
#'  res1 <- xgb_param_calc(expr_mis,response = "all", select_features=NULL, iter = 50, fixed_params =fixed_param, Params_to_tune = c("eta", "alpha"), plots = FALSE)
#' 
#' fixed_param_2 <- list(subsample = 0.7)
#' res2 <- xgb_param_calc(expr_mis,response = "all", select_features=NULL, iter = 50, fixed_params =fixed_param_2, Params_to_tune = c("eta", "alpha"), plots = FALSE)

#' 
#' @export
#' @import ParamHelpers
#' @import mlrMBO
#' @import xgboost
#' @import ggplot2
#' @import xgboost
#' @importFrom stats setNames
#' @importFrom Matrix sparse.model.matrix
#' @importFrom smoof getNumberOfParameters

xgb_param_calc<-function(data,response = NULL, select_features=NULL, num_cores = 1, iter=50, fixed_params=NULL, Params_to_tune, plots = TRUE){
  
  
  ##############################################################################
  ######                        Check 1:  Data checks
  ##############################################################################
  
  
  if(!is.data.frame(data) && !is.matrix(data)){
    stop("Data should be dataframe or matrix. Please check your data.")
  }
  
  data <- as.data.frame(data)
  num.cc <- sum(complete.cases(data))
  
  if (num.cc == 0) {
    stop("No complete cases in this dataset.")
  }
  
  if (num.cc > 0 & num.cc < 20) {
    stop("Less than 20 complete cases. Hyperparameter tuning is not possible")
  }
  
  cc.data <- data[complete.cases(data), ] # only extract complete data. Complete cases will be used for hyperparameter tuning.
  
  Names <- colnames(cc.data) # all feature names
  
  Types <- as.vector(apply(cc.data,2,class)) #type of each variable
  
  na.col <- which(colSums(is.na(data))!=0) #find features with missing values
  
  if(length(na.col)==0){
    stop("There is no varaible with missing data. Parameter optimization will not be performed")
  }
  
  ################################################################################################
  ############          Check 2: Target and training features          ###########################
  ################################################################################################
  
  
  # if user has not chosen from the given options, following message will be displayed
  
  if(!is.null(select_features) &!identical(tolower(select_features), "all") & !all(select_features %in% seq_along(Names))){
    stop("Please check the format of select_features. It can be either NULL, all or index of columns (e.g. c(1,3)) that you want to include as a feature")
  }
  
  if(!is.null(response) &!identical(tolower(response), "all") & !all(response %in% seq_along(Names))){
    stop("Please check the format of response variable. It can be either NULL, all or index of columns (e.g. 1) that you want to include as a response variable")
  }
  
  
  ################################################################################################
  ############          Training and target feature selection form input
  ################################################################################################
  
  
  ###########################################################################
  # response = NULL; different options for select_features
  ###########################################################################
  
  if (is.null(response)) { 
    
    if(is.null(select_features) | identical(tolower(select_features), "all")) {
      
      r.idx <- sample(na.col, size = 1) #randomly sample one feature from columns with missing data and assign it to response_var in next line
      response_var <- Names[r.idx]
      training_features <- setdiff(Names, response_var) # assign all other features (with/without missing data) to training features set
      
    }else if(is.numeric(select_features)){
      
      training_features <- Names[select_features] #features are based on user provided indexes
      response_var <- sample(setdiff(Names[na.col], training_features), size = 1) # response variable is randomly sampled from missing variables other than select_features
    
      } 
    
  } else if (is.numeric(response)) {
    
    ###########################################################################
    # response = index; different option of select_features
    ###########################################################################
    
    if(length(response)>1){ 
      message("You have selected more than one variables for parameter optimization. Please note that only the first index will be used for parameter estimation")
    } 
    
    response_var <- Names[response[1]] # in case of multiple response variables, select first feature from given response variables
    
    if(is.null(select_features) | identical(select_features, "all")) {
      training_features <- setdiff(Names, response_var) # all variables other than response variables
      
    } else if(is.numeric(select_features)){
      if(response %in% select_features){
        message("Response and feature variables should not be same. All variables other than response variable will be used as features")
        training_features <- setdiff(Names, response_var)
      } else {
        training_features <- Names[select_features]
      }
    } 
  } else if(response =="all"){
    

    ###########################################################################
    # response = all; different option of select_features
    ###########################################################################

    
    response_var <- Names[na.col]
    training_features <- NULL
    message("You have chosen all variables (with missing values) to be iteratively used as a target variable for parameter optimization. Please note that select_features will be automatically chosen for this process.")
  } else {
    stop("Please check the format of response variable It can be either NULL, \"all\" or index of column that you want to use for parameter optimization")
  }
  
  params = list()
  
  if(length(response_var)>1){
    
    response_var <- setNames(response_var, response_var)
    message("Performing bayesian optimization iteratively using each vairable as response variable")
    
    if(!is.null(fixed_params)) {
    if(names(fixed_param_2) %in% colnames(data)){
    params <- lapply(response_var,function(x) par_opt(cc.data, setdiff(Names, x), x, iter = iter, fixed_params[[x]]$parameter, param_list = Params_to_tune, plots)) # call to bayesian optimization using assigned features and response variable
    } else {
      stop("When response = all is selected, fixed parameters should be provided for each variable")
    }
    } else {
    params <- lapply(response_var,function(x) par_opt(cc.data, setdiff(Names, x), x, iter = iter, NULL, param_list = Params_to_tune, plots)) # call to bayesian optimization using assigned features and response variable
    }
  } else {
    message("Performing bayesian optimization using ", response_var, " as a response variable and ", training_features, " as features")
    params = par_opt(cc.data, training_features, response_var,iter = iter,fixed_params, param_list = Params_to_tune, plots)
  }
  
  return(params)
  
}

###########################################################################################################################
########################              Bayesian optimization function            ###########################################
###########################################################################################################################


par_opt <- function(cc.data, training_features, response_var, iter = iter,fixed_params, param_list, plots, nfold = 5,  early_stopping_rounds = 10){
  
  response_data <- cc.data[,response_var] # data from response variable 
  
  if(class(response_data)=="character"){
    response_data <- as.factor(response_data)
  } # for next steps in code
  
  response_class <- class(response_data)
  

  ######################################################################################################
  ########       Define variable, objective and evaluation metric type             ####################
  ######################################################################################################
  
  var.type <- switch(response_class[1],
                     numeric = "numeric",
                     integer = "numeric",
                     factor = ifelse(nlevels(response_data) == 2, "binary", "multiclass"),
                     ordered = ifelse(nlevels(response_data) == 2, "binary", "multiclass"),
                     stop(paste0("Unsupported data type: ", response_class[1]))
  )
  
  obj.type <- switch(var.type,
                     numeric = "reg:squarederror",
                     binary =  "binary:logistic",
                     multiclass = "multi:softmax",
                     stop(paste0("Unsupported objective function for ", var.type))
  )
  
  eval.metric <- switch(var.type,
                        numeric = "mape",
                        binary = "logloss",
                        multiclass = "mlogloss",
                        stop(paste0("Unsupported data type: ", obj.type))
  )
  
  ####################################################################################################
  
  #conditions for categorical/multi-level data ; convert categorical data to integers
  
  continue = TRUE
  
  if(var.type != "numeric"){
    response_data <- as.numeric(unlist(response_data)) - 1 # convert categorical data to integers because xgboost can only deal with numbers and integers as labels
    bin.t <- sort(table(response_data), decreasing = TRUE) # distribution of categories
    
    if (is.na(bin.t[2])) {
      msg <- paste("The variable", var, "in the data only have single cateogry. Optimization algorithm can not run for this variable. Default parameters are being assigned to ", response_var)
      message(msg)
      continue = FALSE
    }
  }
  
  ###############################################################################################
  ##########   !!!!!!!!!!           XGBoost object and cross validation             !!!!!!!!!!!!
  ###############################################################################################
  
  p <- length(training_features) + 1
  
  
  if (p == 2) { # convert training data to sparse matrix
    obs.data <- Matrix::sparse.model.matrix(reformulate(training_features, response_var), data = cc.data)
  } else {
    obs.data <- Matrix::sparse.model.matrix(reformulate(training_features, response_var), data = cc.data)[, -1]
  }
  
  dtrain <- xgboost::xgb.DMatrix(obs.data, label = response_data)
  
  if(var.type=="numeric"){
    
    cv_folds = rBayesianOptimization::KFold(response_data, # cross validation 
                                            nfolds= nfold,
                                            seed= TRUE)
    
  } else {
    
    cv_folds = rBayesianOptimization::KFold(response_data, 
                                            nfolds= nfold, stratified = TRUE, #so that distribution of different classes is similar in all folds
                                            seed= TRUE)
    
  }
  
  if(any(lapply(cv_folds, length)<5)){
    msg <- paste("Hyperparameter tuning is not possible for variable", response_var, 
                 " because it does not have sufficient data-points for cross-validation. Default parameters are being assigned to ", response_var)
    message(msg)
    continue = FALSE
  }
  
  ###############################################################################################
  ####################          Parameters definition       ################################
  ###############################################################################################
  
  available_params <- list(
    eta = makeNumericParam("eta", lower = 0, upper = 0.3),
    nrounds = makeIntegerParam("nrounds", lower = 100, upper = 1000),
    
    max_depth = makeIntegerParam("max_depth", lower = 2, upper = 12),
    min_child_weight = makeNumericParam("min_child_weight", lower = 0, upper = 25),
    subsample = makeNumericParam("subsample", lower = 0.5, upper = 1),
    colsample_bytree = makeNumericParam("colsample_bytree", lower = 0.5, upper = 1),
    
    gamma = makeNumericParam("gamma", lower = 1e-09, upper = 10),
    alpha = makeNumericParam("alpha", lower = 1e-09, upper = 100),  
    lambda = makeNumericParam("lambda", lower = 1e-09, upper = 100)  
    
  )
  
  ########################################################################################
  #       Parameters check and formating
  ########################################################################################
  
  browser()
  notfound <- setdiff(param_list, names(available_params))
  
  if(length(notfound) == length(param_list)) {
  
    stop(paste("xgb_param_calc only allows tuning the following parameters:", paste(names(available_params), collapse = ", "),
    "\nPlease provide parameters from this list, e.g. param_list = c(\"eta\", \"nrounds\")."))
  } else if(length(notfound)>0) {
    
    message(paste("xgb_param_calc only allows tuning the following parameters:", paste(names(available_params), collapse = ", "),
               "\nPlease provide parameters from this list, e.g. param_list = c(\"eta\", \"nrounds\"). Provided Parameter(s) ", 
               paste(notfound, collapse = ", "), "are not from the list and will not be used." ))
    
    param_list <- setdiff(param_list, notfound)
  }
  
  # Filter the available parameters based on param_list
  param_set <- do.call(makeParamSet, available_params[param_list])
  
  
  predefined_fixed_params <- list(
    booster = "gbtree",
    objective = obj.type,
    eval_metric = eval.metric
  )
  

  # Combine predefined fixed parameters with user-passed fixed parameters
  if(!is.null(fixed_params)){
    predefined_fixed_params <- c(predefined_fixed_params, fixed_params)
  }
  
  if(!("nrounds" %in% c(param_list, names(predefined_fixed_params)))) {predefined_fixed_params[["nrounds"]] = 100}
  

  ########################################################################################
  ####                Make single objective function
  ########################################################################################
  
  
  obj.fun <- smoof::makeSingleObjectiveFunction(
    name = "xgb_cv_bayes",
    fn = function(x, predefined_fixed_params) {
      # Dynamically build the params list from param_list and fixed_params
        params <- predefined_fixed_params
      for (param in names(x)) {
        params[[param]] <- x[[param]]
      }
        
      if (var.type != "multiclass") {
        
        cv <- xgb.cv(
          params = params,
          data = dtrain,
          nrounds = params[["nrounds"]],
          early_stopping_rounds = early_stopping_rounds,
          folds = cv_folds,
          prediction = FALSE,
          maximize = FALSE,
          showsd = FALSE,
          verbose = 0
        )
      } else {
        N <- length(levels(as.factor(cc.data[, response_var]))) # number of classes
        params[["num_class"]] <- N
        cv <- xgb.cv(
          params = params,
          data = dtrain,
          nrounds = params[["nrounds"]],
          early_stopping_rounds = early_stopping_rounds,
          folds = cv_folds,
          prediction = FALSE,
          maximize = FALSE,
          showsd = FALSE,
          verbose = 0
        )
      }
      eval_column <- names(cv$evaluation_log)[grep("test.*mean", names(cv$evaluation_log))] # extract column containing test metrics
      min(cv$evaluation_log[, eval_column]) # minimum test metrics
    },
    par.set = param_set,
    minimize = TRUE # minimize the evaluation metric
  ) # end of single objective function
  
  
  #################################################################################################################
  ###################          Call to objective function using bayesian optimization      ########################
  #################################################################################################################
  
  control = makeMBOControl(on.surrogate.error = "warn")
  
  num_design_points <- min(5 * getNumberOfParameters(obj.fun), 50)  # Adjust the design points value as needed
  num_focus_points <- 15* getNumberOfParameters(obj.fun)  # Adjust the focus points
  
  
  des <- generateDesign(n = num_design_points,
                        par.set = getParamSet(obj.fun),
                        fun = lhs::randomLHS)    # Generate a design matrix with 'num_design_points' rows for the parameter set of 'obj.fun'
                                                  # using a Latin Hypercube Sampling (LHS) method from the lhs package.
  
  
  control <- setMBOControlTermination(control, iters = iter)    # Set the termination criteria for the MBO (Model-Based Optimization) control object to stop after a specified number of iterations ('iter').
  control <- setMBOControlInfill(control, crit = makeMBOInfillCritEI(), opt.focussearch.points = num_focus_points, filter.proposed.points = TRUE)  # Configure the infill criteria for the MBO control object. Use the Expected Improvement (EI) criterion and set the number of focus search points to 'num_focus_points', and enable filtering of already proposed points.
  
  
  if(continue==TRUE){ # For multi-level data, if there are more than two categories, algorithm will run bayesian optimization. 
    run = mbo(fun = obj.fun,
              control = control,
              design = des,
              more.args = list(predefined_fixed_params=predefined_fixed_params)) # Run the Model-Based Optimization (MBO) process on 'obj.fun' using the specified control settings and initial design 'des'. Additional arguments, such as predefined fixed parameters, are passed via 'more.args'.
    
    
    # plot of iterations on x-axis and evaluation metric on x-axis
    if(plots==TRUE){
    plot_figs <- hyperparameter_plots(run, response_var, training_features,eval.metric, param_list) # returns three different plots
    } else {
    plot_figs <- NA # plotting graphs can take longer so user can opt to run optimization without plots
    }
    best_solution <- run$opt.path$env$path[which.min(run$opt.path$env$path$y),] # extract row which returned in minimum value for loss function
    best_eval <- best_solution$y
    best_parameters <- run$x # best performing parameters
    
    if(!is.null(fixed_params)){ # if user has provided some fixed parameters, than final output includes best performing parameters as well as fixed parameters.
    best_parameters <- c(best_parameters, fixed_params)
    }
  } else {
    
    best_parameters <- list("eta" = 0.3, "lambda" = 1, "alpha" = 0, "gamma" = 0,  # default parameters in xgboost
                            "max_depth" = 6, "min_child_weight" = 1, "nrounds" = 100, 
                            "subsample" = 1, "colsample_bytree" = 1)
  }
  
  list("parameter" = best_parameters, "fig" = plot_figs, "responseVariable" = response_var) # returns bet parameters, figures and response variable
} # end of bayesian optimization function






