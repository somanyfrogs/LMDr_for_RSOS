#' @title func_surrogate
#' @description R function file for 'andsr' (Analysis of Nonlinear Dynamical Systems in R).
#'      This file contains functions related to surrogate data generation.
#'      Intially written on 20200227 by K.Kawatsu.
#'      Last update: 20220502.

#' FFT randomization
#'
#' \code{rand_fft} returns a time-series vector shuffled with
#'      FFT (Fast-Fourier Transformation) randomization.
#' @param amp Numeric, the amplitude of target time series.
#' @param n Integer, the length of target time series.
#' @param sigma Numeric, the standard deviation of target time series.
rand_fft <- function(amp, n, sigma) {
    n2 <- floor(n / 2)

    if (n %% 2 == 0) {
        thetas <- 2 * pi * runif(n2 - 1)
        angles <- c(0, thetas, 0, -rev(thetas))
        recf <- amp * exp(complex(imaginary = angles))
        recf[n2] <- complex(real = sqrt(2) * amp[n2] * cos(runif(1) * 2 * pi))
    } else {
        thetas <- 2 * pi * runif(n2)
        angles <- c(0, thetas, -rev(thetas))
        recf <- amp * exp(complex(imaginary = angles))
    }

    tmp <- Re(fft(recf, inverse = TRUE) / n)
    return(tmp / sd(tmp) * sigma)
}

#' Calculate recurrence matrix
#'
#' \code{gen_rec_mat} returns (binary) 0-1 matrix based on distance matrix.
#' @param dmat A distance matrix.
#' @param s Numeric, which sets neighbourhood radius based on the proportion 's' of dmat.
gen_rec_mat <- function(dmat, s = 0.100) {
    rad <- quantile(dmat, probs = s, na.rm = TRUE)
    ifelse(dmat <= rad, 1, 0)
}

#' Find the 'twin' status in embedding matrix
#'
#' \code{find_twins} returns a list preserving dynamical 'twins' for each time point.
#' Algorithm is adopted from Thiel et al. (2006) Europhys. Lett. and Ushio et al. (2018) Nature,
#'      with an extension that includes NN (nearest neighbour) twin method developed by K. Kawatsu.
#' @param dmat A distance matrix.
#' @param twin_op Strings, which can be "RM" or "NN".
#'      If "RM", twin search is based on recrrence matrix; otherwise is based on NN twin method.
#' @param time_idx Time index vector.
#' @param phase_lock Logical; if TRUE, twin search is limited to the points having timing modular
#'      (defined by 'period') same with the target (a.k.a. phase-locked twin surrogate, Ushio et al. 2018).
#' @param period Integer, the phase period (only works when phase_lock == TRUE).
#' @param lim Integer, the minimum number of twins (only works when twin_op == "RM").
#' @param range Vector, the range of threshold in recurrence matrix
#'      (only works when twin_op == "RM", see also \code{\link{gen_rec_mat}}).
#' @param nn Integer, the number of nearest neighbour points (only works when twin_op != "RM").
find_twins <- function(dmat, twin_op = "RM", time_idx, phase_lock = TRUE, period = 12, lim = 10, range, nn) {
    twins <- as.list(1:nrow(dmat))
    idx <- which(apply(dmat, 1, function(r) !all(is.na(r))))

    if (twin_op == "RM") {          # Twin search follows to recurrence matrix algorithm
        for (s in range) {
            rmat <- gen_rec_mat(dmat, s = s)

            for (i in idx) {
                tmp <- rmat %>% apply(1, function(r) all(r == rmat[i, ], na.rm = TRUE)) %>% which()
                if (phase_lock) tmp <- tmp[(time_idx[tmp] - time_idx[i]) %% period == 0]
                twins[[i]] <- tmp
            }

            if (length(unlist(twins)) - length(twins) >= lim) break;
        }
    } else if (twin_op == "NN") {   # Twin search follows to nearest neighbour algorithm
        tmp <- order(dmat[i, ], na.last = NA)
        twin_end <- nn
        if (phase_lock) tmp <- tmp[(time_idx[tmp] - tmp_idx[i]) %% period == 0]

        if (length(tmp) >= twin_end) {
            repeat {
                if (dmat[i, tmp[twin_end]] < dmat[i, tmp[twin_end + 1]]) break;
                twin_end <- twin_end + 1
            }

            tmp <- tmp[1:twin_en]
        }
    } else {
        warning("Inappropriate option for 'twin_op', NULL list returns.")
    }

    return(twins)
}

#' Calculate twin probabilities
#'
#' \code{twin_probs} returns a list of swapping probabilities of twins for each state.
#' @inheritParams find_twins
#' @param twins A list contains twins for each state (calculated by \code{\link{find_twins}}).
#' @param theta Numeric, the locality parameter (works if twin_op == "NN").
twin_probs <- function(twins, dmat, theta = 1, twin_op) {
    probs <- NULL
    
    if (twin_op == "RM") {
        probs <- lapply(twins, function(v) {rep(1 / length(v), length(v))})
    } else if (twin_op == "NN") {
        probs <- foreach(i = 1:length(twins)) %do% {
            di <- dmat[i, twins[i]]
            wi <- exp(-theta * di) / mean(di)
            wi / sum(wi)
        }
    } else {
        warning("Inappropriate option for 'twin_op', NULL list returns.")
    }

    return(probs)
}

#' Transition procedure to the next point from twins
#'
#' \code{point_nex} proceeds the time series with swapping between original and twin states.
#' If valid data does not exist, this function returns 0.
#' @inheritParams twin_probs
#' @param curr Integer, the position of current state.
#' @param knot Vector or 2-column matrix, which specifies the knot positions in the data.
#' @param probs List of twin probabilities (see \code{\link{twin_probs}}).
point_nex <- function(curr, twins, knot, probs) {
    if (curr == 0) return(0);
    if (!is.matrix(knot)) knot <- matrix(knot, nrow = 1);
    if (ncol(knot) != 2) stop("Inappropriate knot!");

    ## Replace with 0 if the next point is the end of the knot
    if (length(twins[[curr]]) > 1) {
        nex <- sample(x = twins[[curr]], size = 1, prob = probs[[curr]])
    } else {
        nex <- curr
    }

    nex <- if_else(nex %in% knot[, 2], 0, nex + 1)
    return(nex)
}

#' Surrogate data generation with random shuffling
#'
#' \code{gen_surr_rs} returns a column-wise piled surrogate data with RS algorithm.
#' @param X data matrix.
#' @param col Integer, the column position of surrogaate target variable in 'X'.
#' @param knot Vector or 2-column matrix, which specifies the knot positions in the data.
#' @param iter Integer, No. of surrogate generation.
#' @param seed Integer, RNG seed.
#' @export
gen_surr_rs <- function(X, col, knot = matrix(c(1, nrow(X)), nrow = 1), iter = 100, seed = NULL) {
    ## Reinstate system seed after calculation
    if (!is.null(seed)) {
        sysSeed <- .GlobalEnv$.Random.seed
        set.seed(seed)

        on.exit({
            if (!is.null(sysSeed)) {
                .GlobalEnv$.Random.seed
            } else {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        })
    }

    if (!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if (ncol(knot) != 2) stop("Inappropriate knot!")

    surr <- foreach(i = 1:nrow(knot), .combine = rbind) %do% {
        orig <- X[knot[i, 1]:knot[i, 2], col]
        len <- length(orig)
        idx <- which(!is.na(orig))

        foreach(j = 1:iter, .combine = cbind) %do% {
            tmp <- rep(NA, len)
            tmp[idx] <- orig[idx] %>% sample(size = len)
        }
    }

    return(surr)
}

#' Surrogate data generation with FFT
#'
#' \code{gen_surr_fs} returns a column-wise piled surrogate data with FS algorithm.
#' @inheritParams gen_surr_rs
#' @export
gen_surr_fs <- function(X, col, knot = matrix(c(1, nrow(X)), nrow = 1), iter = 100, seed = NULL) {
    ## Reinstate system seed after calculation
    if (!is.null(seed)) {
        sysSeed <- .GlobalEnv$.Random.seed
        set.seed(seed)

        on.exit({
            if (!is.null(sysSeed)) {
                .GlobalEnv$.Random.seed
            } else {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        })
    }

    if (!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if (ncol(knot) != 2) stop("Inappropriate knot!")

    surr <- foreach(i = 1:nrow(knot), .combine = rbind) %do% {
        orig <- X[knot[i, 1]:knot[i, 2], col]
        len <- length(orig)
        idx <- which(!is.na(orig))

        amp <- abs(fft(orig[idx]))
        n <- length(idx)
        sigma <- sd(orig[idx])

        foreach(j = 1:iter, .combine = cbind) %do% {
            tmp <- rep(NA, len)
            tmp[idx] <- rand_fft(amp, n, sigma)
            tmp - mean(tmp[idx]) + mean(orig[idx])
        }
    }

    return(surr)
}

#' Surrogate generation with twin-surrogate
#'
#' \code{gen_surr_ts} returns a column-wise piled surrogate data with TS algorithm.
#' @inheritParams gen_surr_rs
#' @inheritParams find_twins
#' @param theta Numeric, which determines the locality parameter (works if twin_op == "NN").
#' @export
gen_surr_ts <- function(X, col, knot = matrix(c(1, nrow(X)), nrow = 1), iter = 100, seed = NULL,
                        twin_op = "RM", time_idx, phase_lock = TRUE, period = 12, lim = 10,
                        range = c(0.125, seq(0.12, 0.05, -0.01), seq(0.15, 0.20, 0.01), 0.04), theta = 1) {
    ## Reinstate system seed after calculation
    if (!is.null(seed)) {
        sysSeed <- .GlobalEnv$.Random.seed
        set.seed(seed)

        on.exit({
            if (!is.null(sysSeed)) {
                .GlobalEnv$.Random.seed
            } else {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        })
    }

    if (!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if (ncol(knot) != 2) stop("Inappropriate knot!")

    dmat <- get_ned(X)
    twins <- find_twins(dmat, twin_op, time_idx, phase_lock, period, lim, range, ncol(X) + 1)

    if (max(sapply(twins, length)) == 1) {
        twins <- find_twins(dmat, "NN", time_idx, phase_lock, period, range, ncol(X) + 1)
        probs <- twin_probs(twins, dmat, theta, "NN")
    } else {
        probs <- twin_probs(twins, dmat, theta ,twin_op)
    }

    surr <- foreach(i = 1:nrow(knot), .combine = rbind) %do% {
        tmp <- foreach(j = 1:iter, .combine = cbind) %do% {
            if (phase_lock) {
                ids <- sample(which((time_idx - time_idx[knot[i, 1]]) %% period == 0), 1)
            } else {
                ids <- sample(1:nrow(X)[-knot[, 2]], 1)
            }

            for (k in 2:(knot[i, 2] - knot[i, 1] + 1)) ids <- c(ids, point_nex(ids[k - 1], twins, knot, probs))

            ## Return the surrogate if it reaches to the original data length
            if (all(ids != 0)) X[ids, col]
        }

        ## If cannot yeild surrogate, original data used instead
        ## also if the size of surrogate data does not reach the iteration number,
        ## then meet the condition by sampling the surrogate column with replace TRUE.
        if (is.null(tmp)) tmp <- X[knot[i, 1]:knot[i, 2], rep(col, iter)];
        if (iter == 1) tmp <- matrix(tmp, ncol = 1)
        if (ncol(tmp) < iter) tmp <- tmp[, sample(1:ncol(tmp), iter, replace = TRUE)]

        tmp
    }

    return(surr)
}
