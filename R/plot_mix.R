#' @title plot result of miture model
#'
#' @description This function plot result of mixturemodel (type of data is both of continuous and categorical).
#'
#' @param x a data frame containing heterogeneous data. Rows
#' correspond to observations and columns correspond to variables.
#'
#' @param model mixture model from function clustermixmod()
#'
#' @return a graph represents result of model
#'
#' @export plot_mix


plot_mix = function(x, model){
    tsne_object <- Rtsne(x)
    tsne_df <- tsne_object$Y %>%
        data.frame() %>%
        setNames(c("X", "Y")) %>%
        mutate(cluster = factor(model$z))
    ggplot(aes(x = X, y = Y), data = tsne_df) +
        geom_point(aes(color = cluster))
}
