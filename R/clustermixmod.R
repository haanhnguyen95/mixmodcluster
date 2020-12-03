#' @title Mixture models for clustering
#'
#' @description This function fitting a mixture model clustering of multivariate gaussian or multinomial components or the both to a given dataset.
#'
#' @param x a data frame containing quantitative,qualitative or heterogeneous data. Rows
#' correspond to observations and columns correspond to variables.
#'
#' @param K numeric listing the number of clusters.
#'
#' @param itermax numeric. The max of iteration.
#'
#' @param init character. Type of initialization is "random" or "kmeans".
#'
#' @param datatype character. Type of data is "continuous", "categorical" or "mix".
#'
#' @return prop: posteriori probabilities for each individual to belong to each cluster.
#'
#' @return mu, sigma: mean and matrix correlation of each cluster for data "continuous".
#'
#' @return loglik: log-likelihood of model.
#'
#' @return iterfin: the last iteration.
#'
#' @return K:  numeric listing the number of clusters.
#'
#' @return z: cluster of each individual.
#'
#' @return BIC: criterion BIC of model.
#'
#' @return ICL: criterion ICL of model.
#'
#' @export clustermixmod


clustermixmod <- function(x,
                          K,
                          itermax = 100,
                          init = c("random", "kmeans"),
                          datatype = c("continuous", "categorical", "mix")) {
    #To ignore the warnings during usage
    options(warn=-1)
    options("getSymbols.warning4.0"=FALSE)
    #n: nombre de'observations, p: nombre de variables
    n <- nrow(x)
    p <- ncol(x)
    init <- match.arg(init)
    datatype <- match.arg(datatype)

    #Gaussian Mixture Model
    if (datatype == "continuous") {
        #initialisation des objets
        prop <- matrix(NA, itermax + 1, K)
        mu <- array(NA, dim = c(itermax + 1, K, p))
        sigma <- array(NA, dim = c(itermax + 1, K, p, p))
        loglik <- rep(0, itermax + 1)

        #initialisation de l'algo
        #init = random
        if (init == "random") {
            prop[1,] <- rdirichlet(1, par = rep(1, K))
            mu[1, ,] <- x[sample(1:n, K),]
            for (k in 1:K)
                sigma[1, k, ,] <- rWishart(1, 4, var(x))
        }
        #init = kmeans
        if (init == "kmeans") {
            z <- kmeans(x, K)$cluster
            for (k in 1:K) {
                prop[1, k] <- mean(z == k)
                mu[1, k,] <- colMeans(x[which(z == k),])
                sigma[1, k, ,] <- var(x[which(z == k),])
            }
        }
        #calcul de log-likelihood
        for (i in 1:n) {
            tmp <- 0
            for (k in 1:K) {
                tmp <- tmp + prop[1, k] * dmvnorm(x[i,], mean = mu[1, k,], sigma = sigma[1, k, ,])
            }
            loglik[1] <- loglik[1] + log(tmp)
        }
        #Algo EM
        for (iter in 1:itermax) {
            #E-step: calcul tik
            tik <- matrix(NA, n, K)
            for (k in 1:k) {
                tik[, k] <- prop[iter, k] * dmvnorm(x, mean = mu[iter, k,], sigma = sigma[iter, k, ,])
            }
            tik <- tik / rowSums(tik)
            #M-step: update mu, sigma, pk
            for (k in 1:K) {
                nk <- sum(tik[, k])
                prop[iter + 1, k] <- nk / n
                mu[iter + 1, k,] <- colSums(tik[, k] * x) / nk
                sigma[iter + 1, k, ,] <-
                    Reduce("+", lapply(1:n, function(m)
                        tik[m, k] * (x[m,] - mu[iter + 1, k,]) %*% t(x[m,] - mu[iter + 1, k,]))) /
                    nk
            }
            #calcul de log-likelihood
            for (i in 1:n) {
                tmp <- 0
                for (k in 1:K) {
                    tmp <- tmp + prop[iter + 1, k] * dmvnorm(x[i,], mean = mu[iter + 1, k,], sigma =
                                                                 sigma[iter + 1, k, ,])
                }
                loglik[iter + 1] <- loglik[iter + 1] + log(tmp)
            }

            fin <- itermax + 1
            if (abs(loglik[iter + 1] - loglik[iter]) < 1e-4) {
                fin <- iter + 1
                break
            }
        }
        #Cluster pour chaque d'individu
        z <- max.col(tik)
        #BIC
        BIC <- 2 * loglik[fin] -  (K * (p * (p + 1) / 2 + p + 1) - 1) * log(n)
        #ICL
        minus <- 0
        for (i in 1:n) {
            for (k in 1:K) {
                minus <- minus + tik[i, k] * log(tik[i, k])
            }
        }
        ICL <- BIC - minus

        propfn <- prop[fin,]
        mufn <- mu[fin, ,]
        sigmafn <- sigma[fin, , ,]
        loglikfn <- loglik[fin]
        res <- list(
            prop = propfn,
            mu = mufn,
            sigma = sigmafn,
            loglik = loglikfn,
            iterfin = fin,
            K = K,
            z = z,
            BIC = BIC,
            ICL = ICL
        )
        return(res)
    }

    #Latent Class Model
    if (datatype == "categorical") {
        #factorisation des variables
        for (name in names(x)) {
            x[[name]] <- as.factor(x[[name]])
        }
        m <-
            sapply(seq(p), function(p) {
                length(levels(x[, p]))
            })


        prop <- matrix(NA, itermax + 1, K)
        #alpha
        alpha <- array(NA, dim = c(itermax + 1, K, p))
        mode(alpha) <- "list"
        #log-vraissemblance
        loglik <- rep(0, itermax + 1)


        if (init == "random") {
            prop[1,] <- rdirichlet(1, par = rep(1, K))
            for (k in 1:K) {
                for (j in 1:p) {

                    alpha[1, k, j] <- list(rdirichlet(1, par = rep(1, m[j])))
                    names(alpha[1, k, j][[1]]) <- levels(x[, j])
                }
            }
        }

        #calcul de loglik
        for (i in 1:n) {
            tmp <- 0
            for (k in 1:K) {
                #calcul de fk(xi)
                fki <- 1
                for (j in 1:p) {
                    fki <- fki * alpha[1, k, j][[1]][x[i, j]]
                }

                tmp <- tmp + prop[1, k] * fki
            }
            loglik[1] <- loglik[1] + log(tmp)
        }

        # algo EM
        for (iter in 1:itermax) {
            #E step
            tik <- matrix(NA, n, K)
            for (k in 1:K) {
                #calcul des fk(xi)
                fki <- 1
                for (j in 1:p) {
                    fki <- fki * alpha[iter, k, j][[1]][x[, j]]
                }

                tik[, k] <- prop[iter, k] * fki
            }
            tik <- tik / rowSums(tik)


            for (k in 1:K) {
                nk <- sum(tik[, k])
                prop[iter + 1, k] <- nk / n
                for (j in 1:p) {
                    alpha[iter + 1, k, j] <- list(sapply(1:m[j], function(c) {
                        sum(tik[, k] * (x[, j] == levels(x[, j])[c])) / nk
                    }))
                    names(alpha[iter + 1, k, j][[1]]) <- levels(x[, j])
                }
            }

            #calcul de loglik
            for (i in 1:n) {
                tmp <- 0
                for (k in 1:K) {
                    #calcul de fk(xi)
                    fki <- 1
                    for (j in 1:p) {
                        fki <- fki * alpha[iter + 1, k, j][[1]][x[i, j]]
                    }

                    tmp <- tmp + prop[iter + 1, k] * fki
                }
                loglik[iter + 1] <- loglik[iter + 1] + log(tmp)
            }
            fin <- itermax + 1
            if (abs(loglik[iter + 1] - loglik[iter]) < 1e-4) {
                fin <- iter + 1
                break
            }
        }

        #Cluster pour chaque d'individu
        z <- max.col(tik)

        nbparam <- K * sum(sapply(1:p, function(i) m[i]-1))
        #BIC
        BIC <- loglik[fin] - nbparam * log(n)
        #ICL
        minus <- 0
        for (i in 1:n) {
            for (k in 1:K) {
                minus <- minus + tik[i, k] * log(tik[i, k])
            }
        }
        ICL <- BIC - minus

        propfn <- prop[fin,]
        loglikfn <- loglik[fin]
        res <- list(
            prop = propfn,
            K = K,
            loglik = loglikfn,
            z = z,
            BIC = BIC,
            ICL = ICL
        )
        return(res)
    }

    #Mixture model
    if (datatype == "mix") {
        discreteL <- function(x)
            length(unique(x)) < 0.05 * length(x)
        xcont <- as.matrix(x[,!sapply(x, discreteL)])
        xcat <- x[, sapply(x, discreteL)]
        pcat <- ncol(xcat) #nombre de parametres categorielles
        pcont <- ncol(xcont)
        #factorisation des variables categorielles
        for (name in names(xcat)) {
            xcat[[name]] <- as.factor(xcat[[name]])
        }
        #nombre de modalites par variable
        m <- sapply(seq(pcat), function(p) {
            length(levels(xcat[, p]))
        })
        #initialisation des objets
        prop <- matrix(NA, itermax + 1, K)
        mu <- array(NA, dim = c(itermax + 1, K, pcont))
        sigma <- array(NA, dim = c(itermax + 1, K, pcont, pcont))
        loglik <- rep(0, itermax + 1)
        alpha <- array(NA, dim = c(itermax + 1, K, pcat))
        mode(alpha) <- "list"
        #initialisation de l'algo
        #init <- random
        if (init == "random") {
            prop[1, ] <- rdirichlet(1, par = rep(1, K))
            mu[1, , ] <- xcont[sample(1:n, K), ]
            for (k in 1:K)
                sigma[1, k, , ] <- rWishart(1, 4, var(xcont))
            for (k in 1:K) {
                for (j in 1:pcat) {
                    #chaque alpha[,k,j] contient les parametres d'une loi multinomiale et les noms des modalites
                    alpha[1, k, j] <- list(rdirichlet(1, par = rep(1, m[j])))
                    names(alpha[1, k, j][[1]]) <- levels(xcat[, j])
                }
            }
        }
        #calcul de log-likelihood
        for (i in 1:n) {
            tmp <- 0
            for (k in 1:K) {
                #calcul de fk(xi)
                fki <- 1
                for (j in 1:pcat) {
                    fki <- fki * alpha[1, k, j][[1]][xcat[i, j]]
                }
                tmp <- tmp + prop[1, k] * fki * dmvnorm(xcont[i, ], mean = mu[1, k, ], sigma = sigma[1, k, , ])
            }
            loglik[1] <- loglik[1] + log(tmp)
        }

        #Algo EM
        for (iter in 1:itermax) {
            #E step
            tik <- matrix(NA, n, K)
            for (k in 1:K) {
                #calcul des fk(xi)
                fki <- 1
                for (j in 1:pcat) {
                    fki <- fki * alpha[iter, k, j][[1]][xcat[, j]]
                }
                tik[, k] <-
                    prop[iter, k] * fki * dmvnorm(xcont, mean = mu[iter, k, ], sigma = sigma[iter, k, , ])
            }
            tik <- tik / rowSums(tik)
            #M step
            #Mise a jour des parametres
            for (k in 1:K) {
                nk <- sum(tik[, k])
                prop[iter + 1, k] <- nk / n
                for (j in 1:pcat) {
                    alpha[iter + 1, k, j] <- list(sapply(1:m[j], function(c) {
                        sum(tik[, k] * (xcat[, j] == levels(xcat[, j])[c])) / nk
                    }))
                    names(alpha[iter + 1, k, j][[1]]) <-
                        levels(xcat[, j])
                }
                mu[iter + 1, k,] <- colSums(tik[, k] * xcont) / nk
                sigma[iter + 1, k, ,] <-
                    Reduce("+", lapply(1:n, function(m)
                        tik[m, k] * (xcont[m,] - mu[iter + 1, k,]) %*% t(xcont[m,] - mu[iter +
                                                                                            1, k,]))) / nk
            }
            #calcul de log-likelihood
            for (i in 1:n) {
                tmp <- 0
                for (k in 1:K) {
                    #calcul de fk(xi)
                    fki <- 1
                    for (j in 1:pcat) {
                        fki <- fki * alpha[iter + 1, k, j][[1]][xcat[i, j]]
                    }
                    tmp <- tmp + prop[iter + 1, k] * fki * dmvnorm(xcont[i,], mean =
                                                                       mu[iter + 1, k,], sigma = sigma[iter + 1, k, ,])
                }
                loglik[iter + 1] <- loglik[iter + 1] + log(tmp)
            }
            fin <- itermax + 1
            if (abs(loglik[iter + 1] - loglik[iter]) < 1e-4) {
                fin <- iter + 1
                break
            }
        }
        z <- max.col(tik)
        #nombre de parametres
        nbparam <- K * (p * (p + 1) / 2 + p + 1) - 1 + K * sum(sapply(1:pcat, function(i)
            m[i] - 1))
        #BIC
        BIC <- loglik[fin] - nbparam * log(n)
        #ICL
        minus <- rep(0, n)
        for (i in 1:n) {
            minus[i] <- sum(sapply(1:K, function(k) tik[i, k] * log(tik[i, k])))
        }
        ICL <- BIC - sum(minus)
        #resultat
        propfn <- prop[fin,]
        loglikfn <- loglik[fin]
        mufn <- mu[fin, ,]
        sigmafn <- sigma[fin, , ,]
        res <- list(
            prop = propfn,
            mu = mufn,
            sigma = sigmafn,
            loglik = loglikfn,
            iterfin = fin,
            K = K,
            z = z,
            BIC = BIC,
            ICL = ICL
        )
        return(res)
    }
}
