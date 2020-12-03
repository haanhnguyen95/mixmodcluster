## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.show='hold'---------------------------------------------------------
library(mixmodcluster)
library(bayess)
library(mvtnorm)

# Data quantitative
x1 = as.matrix(iris[, 1:4])
# Fit model 
model1 = clustermixmod(x1, 2, itermax = 30, init = "kmeans", datatype = "continuous")
# Plot the result
plot_continuous(x1, model1)

