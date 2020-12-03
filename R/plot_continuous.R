#' @title plot result of Gaussian mixture model
#'
#' @description This function plot result of Gaussian mixture model (type of data is continuous).
#'
#' @param x a data frame containing quantitative data. Rows
#' correspond to observations and columns correspond to variables.
#'
#' @param model model of Gaussian mixture model from function clustermixmod()
#'
#' @return a graph represents result of model
#'
#' @export plot_continuous


plot_continuous = function(x, model){
    pairs(x, lower.panel = NULL, col = model$z)
}
