#' Generate Hyperparameter Trace, Pairwise, and Convergence Plots
#'
#' This function generates three types of plots: hyperparameter trace plot, hyperparameter pairwise plot, and convergence plot.
#'
#' @param run An object containing the result from the `mbo` function.
#' @param response_var A string specifying the name of the response variable.
#' @param training_features A vector of strings specifying the names of the training features.
#' @param eval.metric A string specifying the evaluation metric used in the optimization.
#' 
#' @return A list containing three ggplot2 objects:
#' \describe{
#'   \item{Hyperparameter_trace_plot}{A plot showing the trace of hyperparameter values over iterations.}
#'   \item{Pairwise_parameter_plot}{A pairwise plot of the hyperparameters.}
#'   \item{Convergence_plot}{A plot showing the convergence of the objective function value over iterations.}
#' }
#' 
#' @import ggplot2
#' @import tidyr
#' @import GGally
#' 
#' @export

#' @examples
#' \dontrun{
#' result <- mbo(...)  # Run the Bayesian optimization to get the result
#' plots <- hyperparameter_plots(result, "response_var", c("feature1", "feature2"), "rmse")
#' print(plots$Hyperparameter_trace_plot)
#' print(plots$Pairwise_parameter_plot)
#' print(plots$Convergence_plot)
#' }


hyperparameter_plots <- function(run, response_var, training_features,eval.metric,param_list){

p = list()

###########################################################################################
#########                         Hyper-parameter trace plot
###########################################################################################


# Assuming `run` is the result from `mbo`
trace_data <- as.data.frame(run$opt.path$env$path)

# Assuming `param_list` contains the names of hyperparameters
trace_data_long <- trace_data %>%
  pivot_longer(cols = param_list, names_to = "Parameter", values_to = "Value")

# Add an Iteration column to use for the x-axis
trace_data_long$Iteration <- rep(seq_len(nrow(trace_data)), each = length(param_list))


p[[1]] <- ggplot(trace_data_long, aes(x = Iteration, y = Value, color = Parameter)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Parameter, scales = "free_y") +
  labs(title = sprintf("Response Variable: %s , Features: %s", 
                       response_var, paste(training_features,collapse = ","))) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")


###########################################################################################
#########                         Hyper-parameter pairwise plot
###########################################################################################


# Assuming `param_list` contains the names of hyperparameters
p[[2]] <- ggpairs(trace_data, columns = param_list) + 
  labs(title = sprintf("Response Variable: %s , Features: %s", 
                       response_var, paste(training_features,collapse = ","))) +
  theme(plot.title = element_text(hjust = 0.5))



###########################################################################################
#########                         Convergence plot
###########################################################################################

trace_data$Round <- as.numeric(rownames(trace_data))

p[[3]] <- trace_data %>% 
  mutate(Round = factor(Round, levels = 1:max(Round))) %>%
  ggplot(aes(x= Round, y= y)) + 
  geom_point() +
  labs(title = sprintf("Response Variable: %s , Features: %s", response_var, paste(training_features,collapse = ",")))+
  ylab(sprintf("Evaluation score - %s", eval.metric)) + theme(plot.title = element_text(hjust = 0.5),
                                                              axis.text.x = element_text(angle = 30, hjust = 1, size = 6)
  )


#combine all plots and assign names
parameter_plots <- list("Hyperparameter_trace_plot" = p[[1]], 
                        "Pairwise_parameter_plot" = p[[2]], 
                        "Convergence_plot" = p[[3]])


return(parameter_plots)

}
