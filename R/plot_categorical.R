#' @title plot result of Latent Class model
#'
#' @description This function plot result of Latent Class model (type of data is categorical).
#'
#' @param x a data frame containing quantitative,qualitative or heterogeneous data. Rows
#' correspond to observations and columns correspond to variables.
#'
#' @param model model of Latent Class Model from function clustermixmod()
#'
#' @return a graph represents result of model
#'
#' @export plot_categorical


plot_categorical = function(x, model){
    plot(x,col=model$z,cex=0.8,main='final partition')
}
