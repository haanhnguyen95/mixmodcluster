# mixmodcluster
*mixmodcluster* is a R package for fitting a mixture model clustering of multivariate gaussian 
or multinomial components or the both to a given dataset.

**Author :** Ha Anh NGUYEN

The mixmodcluster package contains 2 types of function:

- The first one is function *clustermixmod()*. This function is used to fit mixture models of a given quantitative,qualitative or heterogeneous data. For quantitative data (which the parameter of function: `datatype = "continuous"`), we use Gaussian mixture model (GMM). For qualitative data (which the parameter of function: `datatype = "categorical"`), we use Latent class model (LCM). For heterogeneous data, we use mixture model of Gaussian and Latent class model.

- The second consists 3 plot functions which correspond to 3 types of model. For quantitative data (GMM), we use function *plot_continuous()*. For qualitative data (LCM), we use function *plot_categorical()*. For heterogeneous data, we use function *plot_mix()*.
