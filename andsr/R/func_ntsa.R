#' @title func_ntsa
#' @description R function file for 'andsr' (Analysis of Nonlinear Dynamical Systems in R).
#'      This file contains functions related to nonlinear time-series analysis.
#'      Initially written on 20210706 by K.Kawatsu.
#'      Last update: 20220502.

#' Find knot positions in data
#'
#' \code{find_knot} returns a 2-column matrix, where start and end point of an i'th knot in ts.
#' @param ts data.frame or tibble contains the time-series data.
#' @param key_col Numeric or strings, which contains the labels of the data group.
#' @param time_col Numeric or strings, which specifies the time-point column.
#' @param diff Numeric, which sets the value of increment step in 'time_col'.
#' @export
find_knot <- function(ts, key_col = NULL, time_col = 1, diff = 1) {
    if(is.null(key_col)) {
        tmp <- ts |> pull(time_col) |> (\(.x) .x - diff)()
        knot <- cbind(s = which(tmp == min(tmp)), e = which(tmp == max(tmp)))
    } else {
        keys <- ts |> pull(key_col)
        knot <- foreach(k = unique(keys), .combine = rbind) %do% cbind(s = min(which(keys == k)), e = max(which(keys == k)))
    }

    return(knot)
}

#' Generate valid time index from knot data
#'
#' \code{gen_valid_idx} retunrs a vector, which contains positions of valid row positions in knot.
#' @param knot vector or 2-column matrix, which specifies the knot position in data matrix.
#' @param idx_all A vector contains valid position in the whole data matrix.
gen_valid_idx <- function(knot, idx_all) {
    tmp <- foreach(i = 1:nrow(knot), .combine = c) %do% knot[i, 1]:knot[i, 2]
    return(tmp[tmp %in% idx_all])
}

#' Embedding matrinx generation with multivariate and multi-timelags
#'
#' \code{gen_emat} return an embedding matrix.
#' @param ts vector, matrix, data.frame or tibble
#'      which contains time series to make an embedding matrix.
#' @param cols Numeric or charactor vector,
#'      which selects the ts's column for the embedding matrix reconstruction.
#' @param lags Integer vector, which sets
#'      the time delay value for each coordinate in embedding matrix.
#' @param knot vector or 2-column matrix, which specifies the knot position in ts.
#' @export
gen_emat <- function(ts, cols, lags, knot = matrix(c(1, nrow(ts)), nrow = 1)) {
    ### Check whether knot is provided appropriately
    if(!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if(ncol(knot) != 2) stop("Inappropriate style knot!")

    ### Check wheter ts is provided appropriately
    if(!is_tibble(ts)) {
        if(!(is.data.frame(ts) | is.matrix(ts))) ts <- matrix(ts, ncol = 1)
        ts <- ts %>% set_cnames(str_c("x", 1:ncol(.))) %>% as_tibble()
    }

    emat <- matrix(NA, nrow = max(knot) - min(knot) + 1, ncol = length(cols))

    for(i in 1:nrow(knot)) {
        idx <- knot[i, 1]:knot[i, 2]

        for(j in 1:length(cols)) {
            tmp <- ts %>% slice(idx) %>% pull(cols[j])
            emat[idx, j] <- shift(tmp, lags[j])
        }
    }

    return(emat)
}

#' Calculation of Multiview-Distance (MVD) matrix
#'
#' \code{get_mvd} returns MVD matrix calculated with an ensemble of euclidean
#'      distance matrix of the top-k embeeding in multivariate simplex projection.
#' Algorithm is adopted from Chang et al. (2021), Ecol Lett.
#' For more details, see also \code{\link[rEDM]{multiview}}.
#'
#' @inheritParams gen_emat
#' @param E Integer, which sets the embedding dimension.
#' @param Tp Integer, which sets the prediction-time horizon.
#' @param lmax Integer, which sets the maximum value of timelags.
#' @param threadNo Integer, which is the No. of cores used for parallel computing.
#' @export
get_mvd <- function(ts, cols, tar = 1, knot = matrix(c(1, nrow(ts)), nrow = 1), E, Tp = 1, lmax, threadNo = detectCores()) {
    ## Check whether knot is provided appropriately
    if(!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if(ncol(knot) != 2) stop("Inappropriate style knot provided!")

    ## Perform MVE (Multiview Embedding) with rEDM::Multiview
    columns <- str_c(cols, collapse = " ")
    target <- if_else(is.numeric(tar), columns[tar], tar)
    mve <- ts |> select(all_of(cols)) |> mutate(Index = 1:nrow(ts), .before = everything()) |>
        rEDM::Multiview(dataFrame = _, lib = knot, pred = knot, D = E, E = lmax, Tp = Tp, columns = columns, target = target, parameterList = TRUE, numThreads = threadNo)

    ## Summarize MVE result
    wc <- mve$View |> dplyr::transmute(rho = rho / sum(rho)) |> pull(rho)
    cmbs <- mve$View |> dplyr::transmute(across(starts_with("name"), ~str_sub(.x, 1, -6))) |> as.matrix()
    lags <- mve$View |> dplyr::transmute(across(starts_with("name"), ~as.numeric(str_sub(.x, -3, -2)))) |> as.matrix()

    ## Set environment for parallel computing
    cl <- makeCluster(spec = threadNo, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    mvd <- foreach(i = 1:nrow(cmbs), .combine = '+', .init = matrix(0, nrow = nrow(ts), ncol = nrow(ts))) %dopar% {
        mi <- ts |> andsr::gen_emat(cols = cmbs[i, ], lags = lags[i, ], knot = knot)
        di <- andsr::get_ned(mi)
        wc[i] * di
    }

    return(mvd)
}

#' R wrapper function of lmdSmplx in ands_ntsa.cpp
#'
#' \code{LMDsmplx} returns stats and predictions by LMD-based simplex projection.
#' @param ts A vector, matrix, data.frame or tibble, which contains time-series data.
#' @param col Integer, which specifides the position of prediction target in ts.
#' @param lib Vector or 2-column matrix, which sets the knot positions in library data.
#' @param pred Vector or 2-column matrix, which sets the knot position in test data.
#' @param E Integer, embedding dimension.
#' @param Tp Integer, which is the prediction-time horizon.
#' @param theta A matrix, where each element controls the locality of coordinates in ts.
#' @export
LMDsmplx <- function(ts, col, lib = matrix(c(1, nrow(ts)), nrow = 1), pred = NULL, E, Tp = 1, theta = NULL, both = TRUE) {
    ## Check whether lib is provided appropriately
    if(!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if(ncol(lib) != 2) stop("Inappropriate lib!")

    if(is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    if(is.null(theta)) theta <- diag(E)

    ## Implementation of Simplex Projection algorithm
    smplx <- function(tp) {
        X <- ts |> gen_emat(cols = rep(col, E), lags = 0:(sign(tp) * -(E - 1)), knot = knot)
        y <- ts |> gen_emat(cols = col, lags = tp, knot = knot)
        dmat <- get_lmd(X, theta)
        y_ <- rep(NA, nrow(X))

        idx_all <- cbind(y, X) |> complete.cases() |> which()
        idx_l <- gen_valid_idx(lib, idx_all)
        idx_p <- gen_valid_idx(pred, idx_all)

        tryCatch({
            y_[idx_p] <- lmdSmplx(X, y, dmat, idx_l - 1, idx_p - 1, E + 1)
            return(y_)
        }, error = function(e) return(y_))
    }

    fwd <- smplx(Tp) |> shift(-Tp)

    if (both) {
        bwd <- smplx(-Tp) |> shift(Tp)
        fwd <- cbind(fwd, bwd) |> rowMeans(na.rm = TRUE)
    }

    obs <- ts |> gen_emat(cols = col, lags = 0, knot = knot) |> as.numeric()
    output <- tibble(obs = obs, pred = fwd) |> list()
    tibble(rho = get_rho(obs, fwd), mae = get_mae(obs, fwd), rmse = get_rmse(obs, fwd), output = output)
}

#' Grid search of best embedding dimension with simplex projection
#'
#' \code{find_best_dim} returns the best embedding dimension.
#' @inheritParams LMDsmplx
#' @param range Integer vector, which sets the range of dimension search.
#' @param Tps Integer vector, which setse the prediction time horison for each variable.
#' @param criterion Strings, which switches the cost function (Rho, MAE or RMSE)
#'      for the dimension search (default = "rmse").
#' @param threadNo Integer, which sets the thread number of parallel computing.
#' @export
find_best_dim <- function(ts, cols, lib = matrix(c(1, nrow(ts)), nrow = 1), pred = NULL, range = 1:10,
                          Tps = rep(1, length(cols)), criterion = "rmse", both = TRUE, threadNo = detectCores()) {
    ## Check whether lib is provided appropriately
    if(!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if(ncol(lib) != 2) stop("Inappropriate lib!")

    if(is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    ## Implementation of Simplex Projection algorithm
    smplx <- function(col, E, Tp) {
        X <- ts |> gen_emat(cols = rep(col, E), lags = 0:(sign(Tp) * -(E - 1)), knot = knot)
        y <- ts |> gen_emat(cols = col, lags = Tp, knot = knot)
        dmat <- get_ned(X)
        y_ <- rep(NA, nrow(X))

        idx_all <- cbind(y, X) |> complete.cases() |> which()
        idx_l <- gen_valid_idx(lib, idx_all)
        idx_p <- gen_valid_idx(pred, idx_all)

        tryCatch({
            y_[idx_p] <- lmdSmplx(X, y, dmat, idx_l - 1, idx_p - 1, E + 1)
            return(tibble(E = E, rho = -get_rho(y, y_), mae = get_mae(y, y_), rmse = get_rmse(y, y_)))
        }, error = function(e) return(tibble(E = E, rho = NA, mae = NA, rmse = NA)))
    }

    ## Set environment for parallel computing
    cl <- makeCluster(spec = threadNo, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    ## Find best embedding dimension for each variable
    foreach(i = 1:length(cols), .combine = rbind) %:% foreach(E = range, .combine = rbind) %dopar% {
        fwd <- smplx(cols[i], E, Tps[i])

        if (both) {
            bwd <- smplx(cols[i], E, -Tps[i])
            fwd <- (fwd + bwd) / 2
        }

        fwd |> mutate(var = cols[i], Tp = Tps[i])
    } |> group_by(var) |> arrange(!!!rlang::syms(criterion)) %>% slice(1) %>%
        mutate(rho = -rho, var = factor(var, levels = cols)) |> ungroup() |> select(var, Tp, E, rho, mae, rmse) |> arrange(var)
}

#' S-map with elastic net (a.k.a. Regularized S-map)
#'
#' \code{smap_net} returns the (tibble) result of Regularized S-map.
#' Algorithm is adopted from Cenci et al. (2019), MEE.
#'
#' @inheritParams find_best_dim
#' @param X A matrix, which is constructed by the function \code{\link{gen_emat}}.
#' @param col Integer, which specifies the position of prediction target in X.
#' @param Tp Integer, which is the prediction-time horizon.
#' @param dmat Distance matrix to calculate local weights.
#'      If dmat is obtained by \code{\link{get_mvd}}, the procedure follows
#'      the MDR (MultiView-Distance Regularized) S-map algorithm (see Chang et al. 2021, Ecol Lett.)
#'      If NULL< then Normal-Eculidean Distance is instead used for the analysis.
#' @param range Numeric vector, which sets the range of theta (uniform localisation parameter) for S-map.
#' @param seed Integer, which sets the RNG seed.
#' @param s Strings, which swithces the criterion of the best theta in elastic net.
#'      (for more detail, see original function \code{\link[glmnet]{cv.glmnet}}).
#' @param lambda Numeric vector, which sets the sequence of penalty strength for glmnet.
#'      If NULL (default), then the adaptive search with cv.glmnet is adopted.
#' @param alpha Numeric, which controls the Elastic net parameter.
#'      \code{alpha = 0} corresponds to the ridge penalty and \code{alpha = 1} is the lasso one.
#'      (for more detail, see original function \code{\link[glmnet]{cv.glmnet}}).
#' @export
smap_net <- function(X, col = 1, lib = matrix(c(1, nrow(X)), nrow = 1), pred = NULL, Tp = 1, dmat = NULL, range = seq(0, 10, 1),
                     seed = NULL, threadNo = detectCores(), criterion = "rmse", s = "lambda.1se", lambda = NULL, alpha = 0.0) {
    ## Reinstate system seed after calculation
    if(!is.null(seed)) {
        sysSeed <- .GlobalEnv$.Random.seed
        set.seed(seed)

        on.exit({
            if(!is.null(sysSeed)) {
                .GlobalEnv$.Random.seed <- sysSeed
            } else {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        })
    }

    ## Check whether lib is provided appropriately
    if(!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if(ncol(lib) != 2) stop("Inappropraite lib style!")

    if(is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    if(is.null(dmat) | all(is.na(dmat))) dmat <- get_ned(X);

    ## Preparation for Elastic-net analysis with R package 'glmnet'
    dim <- dim(X)
    y <- X |> gen_emat(cols = col, lags = Tp, knot = knot)
    idx_all <- cbind(y, X) |> complete.cases() |> (\(.x) .x & apply(dmat, 1, \(r) sum(is.na(r)) != dim[1] - 1))() |> which()
    idx_l <- gen_valid_idx(lib, idx_all)
    idx_p <- gen_valid_idx(pred, idx_all)

    dbar <- rowMeans(dmat[, idx_l], na.rm = TRUE)

    ## Set environment for parallel computing
    cl <- makeCluster(spec = threadNo, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    seeds <- sample(32768, dim[1], replace = TRUE)

    ## grid-search for best theta
    op <- foreach(theta = range, .combine = rbind) %dopar% {
        tryCatch({
            ## Make weight matrix
            coef <- matrix(NA, nrow = dim[1], ncol = dim[2] + 1)
            wmat <- exp(-theta * dmat / dbar)
            diag(wmat) <- 0

            ## Sequential estimation of Jacobian
            coef[idx_p, ] <- foreach(t = idx_p, .combine = rbind) %do% {
                set.seed(seeds[t])
                fit <- cv.glmnet(x = X[idx_l, ], y = y[idx_l], weights = wmat[t, idx_l], lambda = lambda, alpha = alpha)
                coef(fit, s = s) |> as.numeric() |> (\(.x) .x[c(2:(dim[2] + 1), 1)])()
            }

            y_ <- getPred(coef, cbind(X, 1))
            output <- tibble(obs = as.numeric(y), pred = y_) |> list()
            coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()

            tibble(theta = theta, rho = -get_rho(y, y_), mae = get_mae(y, y_), rmse = get_rmse(y, y_),
                   output = output, coef = coef)
        }, error = function(e) {
            output <- tibble(obs = as.numeric(y), pred = NA) |> list()
            coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()
            tibble(theta = theta, rho = NA, mae = NA, rmse = NA, output = output, coef = coef)
        })
    } |> arrange(!!!rlang::syms(criterion)) |> slice(1) |> mutate(rho = -rho)

    return(op)
}

#' S-map with multivariate embedding (a.k.a. multivariate S-map)
#'
#' \code{smap_mulvar} returns the (tibble) result of multivariate S-map.
#' Algorithm is adopted from Deyle et al. (2016) PRSB.
#'
#' @inheritParams find_best_dim
#' @inheritParams smap_net
#' @export
smap_mulvar <- function(X, col = 1, lib = matrix(c(1, nrow(X)), nrow = 1), pred = NULL, Tp = 1, dmat = NULL,
                        range = seq(0, 10, 1), threadNo = detectCores(), criterion = "rmse") {
    ## Check whether lib is provided approximately
    if(!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if(ncol(lib) != 2) stop("Inappropriate lib style!")

    if(is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    ## Preparation for multivariate S-map analysis
    dim <- dim(X)
    y <- X |> gen_emat(cols = col, lags = Tp, knot = knot)
    idx_all <- cbind(y, X) |> complete.cases() %>% which()
    idx_l <- gen_valid_idx(lib, idx_all)
    idx_p <- gen_valid_idx(pred, idx_all)

    if(is.null(dmat) | all(is.na(dmat))) dmat <- get_ned(X);
    dbar <- rowMeans(dmat[, idx_l], na.rm = TRUE)

    ## Set environment for parallel computing
    cl <- makeCluster(spec = threadNo, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    ## Grid-search for best theta
    op <- foreach(theta = range, .combine = rbind) %dopar% {
        tryCatch({
            ## Make weight matrix
            coef <- matrix(NA, nrow = dim[1], ncol = dim[2] + 1)
            wmat <- exp(-theta * dmat / dbar)
            diag(wmat) <- 0

            ## Sequential estimation of Jacobian
            coef[idx_p, ] <- smapSVD(cbind(X, 1)[idx_l, ], y[idx_l], wmat[idx_p, idx_l])
            y_ <- getPred(coef, cbind(X, 1))
            output <- tibble(obs = as.numeric(y), pred = y_) |> list()
            coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()

            tibble(theta = theta, rho = -get_rho(y, y_), mae = get_mae(y, y_), rmse = get_rmse(y, y_),
                   output = output, coef = coef)
        }, error = function(e) {
            output <- tibble(obs = as.numeric(y), pred = NA) |> list()
            coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()
            tibble(theta = theta, rho = NA, mae = NA, rmse = NA, output = output, coef = coef)
        })
    } |> arrange(!!!rlang::syms(criterion)) |> slice(1) |> mutate(rho = -rho)

    return(op)
}

#' R wrapper function of lmdSMap in ands_ntsa.cpp
#'
#' \code{LMDmap} returns Jacobian elements estimated by LMD-based nonlinear regression.
#'
#' @inheritParams LMDsmplx
#' @param X An embedding matrix, which is constructed by the function \code{\link{gen_emat}}.
#' @param method Integer, which switches the LMD-based mapping solver.
#'      If col == 0, y = C X is solved with SVD, otherwise is solved with L2 norm penelization.
#' @param theta A matrix, where each element controls the locality of coordinates in X.
#' @param lambda Numeric, which controls the penalization strength in Ridge regression.
#' @export
LMDmap <- function(X, col = 1, lib = matrix(c(1, nrow(X)), nrow = 1), pred = NULL, Tp = 1,
                   method, theta, lambda) {
    ## Check whether lib is provided appropriately
    if(!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if(ncol(lib) != 2) stop("Inappropriate lib style!")

    if(is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    ## Preparation for LMD-based nonlinear regression
    dim <- dim(X)
    coef <- matrix(NA, nrow = dim[1], ncol = dim[2] + 1)
    y <- X |> gen_emat(cols = col, lags = Tp, knot = knot)
    idx_all <- cbind(y, X) |> complete.cases() |> which()
    idx_l <- gen_valid_idx(lib, idx_all)
    idx_p <- gen_valid_idx(pred, idx_all)

    tryCatch({
        coef[idx_p, ] <- lmdSMap(cbind(X, 1), y, idx_l - 1, idx_p - 1, method, theta, lambda);
        y_ <- getPred(coef, cbind(X, 1))
        output <- tibble(obs = as.numeric(y), pred = y_) |> list()
        coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()

        return(tibble(rho = get_rho(y, y_), mae = get_mae(y, y_), rmse = get_rmse(y, y_),
                      output = output, coef = coef))
    }, error = function(e) {
        output <- tibble(obs = as.numeric(y), pred = NA) |> list()
        coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()
        return(tibble(rho = NA, mae = NA, rmse = NA, output = output, coef = coef))
    })
}

#' R wrapper function of funcTPSA in ands_ntsa.cpp
#'
#' \code{find_best_theta} returns LMDr results at the best locality parameter.
#'      The best locality parameter is searched with TPSA (Temperature-Parallel Simulated Annealing),
#'      which is implemented wich C++ file (ands_ntsa.cpp).
#' 
#' @inheritParams LMDmap
#' @param gSeed Integer, which sets the RNG sed used in 'funcTPSA'.
#' @param threadNo Integer, which sets the No. of threads in 'funcTPSA'.
#' @param iterLim Integer, the maximum iteration of TPSA.
#' @param tsLength Integer, the time length in single SA run.
#' @param diag Integer, if 0, the algorithm seeks the optimal values for only diagonal elementse in Theta (default 0).
#' @param sigmas Numeric vector, which sets the SD of normal_distribution used in 'funcTPSA'.
#' @param temps Numeric vector, which sets the range of temperature used in 'funcTPSA'.
#' @export
find_best_theta <- function(X, col = 1, lib = matrix(c(1, nrow(X)), nrow = 1), pred = NULL, Tp = 1,
                            gSeed = NA, threadNo = detectCores(), iterLim = 100, tsLength = 100,
                            criterion = "rmse", method = 0, diag = 0, sigmas = c(1.0, 0.01), temps = c(1e-10, 1e-2)) {
    criterion <- case_when(criterion == "rmse" ~ 0, criterion == "mae" ~ 1, TRUE ~ 2)
    if(!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if(ncol(lib) != 2) stop("Inappropriate lib!");

    if(is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    ## Preparation for TPSA of LMD-based nonlinear regression
    dim <- dim(X)
    coef <- matrix(NA, nrow = dim[1], ncol = dim[2] + 1)
    y <- X |> gen_emat(cols = col, lags = Tp, knot = knot)
    idx_all <- cbind(y, X) |> complete.cases() |> which()
    idx_l <- gen_valid_idx(lib, idx_all)
    idx_p <- gen_valid_idx(pred, idx_all)

    tryCatch({
        op <- funcTPSA(gSeed = gSeed, threadNo = threadNo, iterLim = iterLim, tsLength = tsLength,
                       criterion = criterion, method = method, diag = diag, sigmas = sigmas, temps = temps,
                       cbind(X, 1), y, idx_l - 1, idx_p - 1)

        coef[idx_p, ] <- lmdSMap(cbind(X, 1), y, idx_l - 1, idx_p - 1, method, op$Theta, op$Lambda)
        y_ <- getPred(coef, cbind(X, 1))
        output <- tibble(obs = as.numeric(y), pred = y_) |> list()
        coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()

        return(tibble(theta = list(op$Theta), lambda = op$Lambda, rho = get_rho(y, y_),
                      mae = get_mae(y, y_), rmse = get_rmse(y, y_), output = output, coef = coef))
    }, error = function(e) {
        output <- tibble(obs = as.numeric(y), pred = NA) |> list()
        coef <- coef |> set_cnames(str_c("C", c(1:dim[2], 0))) |> as_tibble() |> list()
        return(tibble(theta = NA, lambda = NA, rho = NA, mae = NA, rmse = NA, output = output, coef = coef))
    })
}

