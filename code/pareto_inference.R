### score matching inference for the skewed Brown-Resnick and Truncated extremal-t Process ###
################################################################################################
fit.scoreMatching <- function (par2, obs, loc, vario.fun, weightFun = NULL, dWeightFun = NULL, nCores = 1L, ST = FALSE, ...){
    ellipsis <- list(...)
    if ("weigthFun" %in% names(ellipsis) && is.null(weightFun)) {
        weightFun <- ellipsis[["weigthFun"]]
        ellipsis[["weigthFun"]] <- NULL
    }
    if ("dWeigthFun" %in% names(ellipsis) && is.null(dWeightFun)) {
        dWeightFun <- ellipsis[["dWeigthFun"]]
        ellipsis[["dWeigthFun"]] <- NULL
    }
    if (is.null(weightFun)) {
        stop("`weightFun` argument missing, with no default value.")
    }
    if (is.null(dWeightFun)) {
        stop("`dWeightFun` argument missing, with no default value.")
    }
    if (!inherits(obs, "list") || length(obs) < 1 || !inherits(obs[[1]], 
        c("numeric", "integer"))) {
        stop("`obs` must be a list of vectors.")
    }
    if (!inherits(weightFun, "function")) {
        stop("`weightFun` must be a function.")
    }
    if (!inherits(dWeightFun, "function")) {
        stop("`dWeightFun` must be a function.")
    }
    if (!is.numeric(nCores) || nCores < 1) {
        stop("`nCores` must a positive number of cores.")
    }
    n <- length(obs)
    dim <- nrow(loc)
    if (!ST) {
        SigmaS = vario.fun(loc, par2)
        computeScores = function(i) {
            obs.i = .subset2(obs, i)
            ind = !is.na(obs.i)
            dim = sum(ind)
            obs.i = obs.i[ind]
            sigmaInv <- MASS::ginv(SigmaS[ind, ind])
            sigma <- diag(SigmaS[ind, ind])
            q <- rowSums(sigmaInv)
            A <- sigmaInv - q %*% t(q)/sum(q)
            zeroDiagA <- A
            diag(zeroDiagA) <- 0
            mtp <- 2 * q/(sum(q)) + 2 + sigmaInv %*% sigma - 
                (q %*% t(q) %*% sigma)/(sum(q))
            gradient <- -1/2 * ((A + t(A)) %*% log(obs.i)) * 
                (1/obs.i) - 1/2 * (1/obs.i) * mtp
            diagHessian <- -1/2 * diag(A + t(A)) * (1/obs.i^2) + 
                1/2 * ((A + t(A)) %*% log(obs.i)) * (1/obs.i)^2 + 
                1/2 * (1/obs.i)^2 * mtp
            weights <- do.call(what = "weightFun", args = c(ellipsis, 
                x = list(obs.i)))
            dWeights <- do.call(what = "dWeightFun", args = c(ellipsis, 
                x = list(obs.i)))
            sum(2 * (weights * dWeights) * gradient + weights^2 * 
                diagHessian + 1/2 * weights^2 * gradient^2)
        }
    }
    else {
        computeScores = function(i) {
            obs.i = .subset2(obs, i)
            ind = !is.na(obs.i)
            dim = sum(ind)
            obs.i = obs.i[ind]
            SigmaS = vario.fun(loc[ind, ], par2, i)
            sigmaInv <- tryCatch({
                MASS::ginv(SigmaS)
            }, warning = function(war) {
                print(war)
            }, error = function(err) {
                print(err)
                print(range(SigmaS, na.rm = TRUE))
                print(par2)
                browser()
            })
            sigma <- diag(SigmaS)
            q <- rowSums(sigmaInv)
            A <- sigmaInv - q %*% t(q)/sum(q)
            zeroDiagA <- A
            diag(zeroDiagA) <- 0
            mtp <- 2 * q/(sum(q)) + 2 + sigmaInv %*% sigma - 
                (q %*% t(q) %*% sigma)/(sum(q))
            gradient <- -1/2 * ((A + t(A)) %*% log(obs.i)) * 
                (1/obs.i) - 1/2 * (1/obs.i) * mtp
            diagHessian <- -1/2 * diag(A + t(A)) * (1/obs.i^2) + 
                1/2 * ((A + t(A)) %*% log(obs.i)) * (1/obs.i)^2 + 
                1/2 * (1/obs.i)^2 * mtp
            weights <- do.call(what = "weightFun", args = c(ellipsis, 
                x = list(obs.i)))
            dWeights <- do.call(what = "dWeightFun", args = c(ellipsis, 
                x = list(obs.i)))
            sum(2 * (weights * dWeights) * gradient + weights^2 * 
                diagHessian + 1/2 * weights^2 * gradient^2)
        }
    }
    if (nCores > 1) {
        scores <- parallel::mclapply(1:n, computeScores, mc.cores = nCores)
    }
    else {
        scores <- lapply(1:n, computeScores)
    }
    return(sum(unlist(scores))/n)
}

