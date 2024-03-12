#' @title func_common
#' @description R function file for 'andsr' (Analysis of Nonlinear Dynamical Systems in R).
#'      This file contains common functions throuhout \pkg{andsr}.
#'      Initially written on 20210319 by K.Kawatsu.
#'      Last update: 20230120.

#' Modified version of as_tibble
#'
#' @param x A data frame, list, matrix, or other object that could reasonably be coerced to a tibble.
#' @param coln A strings vector specifying column name.
#' @export
as_tibbler <- function(x, coln = NULL) {
    suppressMessages(tbl <- as_tibble(x, .name_repair = "unique"))
    if(!is.null(coln)) tbl <- tbl |> setNames(coln)
    return(tbl)
}

#' R wrapper function of getRMSE in ands_common.cpp
#'
#' \code{get_rmse} returns the RMSE (Root Mean Squared Error) between x and y.
#' @param x Numeric vector 1.
#' @param y Numeric vector 2.
#' @export
get_rmse <- function(x, y) {
    idx <- complete.cases(cbind(x, y))
    res <- if_else(any(idx), getRMSE(x[idx], y[idx]), NA_real_)
    return(res)
}

#' R wrapper function of getMAE in ands_common.cpp
#'
#' \code{get_mae} returns the MAE (Mean Absolute Error) between x and y.
#' @inheritParams get_rmse
#' @export
get_mae <- function(x, y) {
    idx <- complete.cases(cbind(x, y))
    res <- if_else(any(idx), getMAE(x[idx], y[idx]), NA_real_)
    return(res)
}

#' R wrapper function of getRho in ands_common.cpp
#'
#' \code{get_rho} returns the correlation coefficient between x and y.
#' @inheritParams get_rmse
#' @export
get_rho <- function(x, y) {
    idx <- complete.cases(cbind(x, y))
    res <- if_else(any(idx), getRho(x[idx], y[idx]), NA_real_)
    return(res)
}

#' Calculation of CV (Coefficient of Variance)
#'
#' \code{get_cv} returns the CV of input vector x.
#' @param x Numeric, input vector.
#' @export
get_cv <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE);

#' R wrapper function of getPred in ands_common.cpp
#'
#' \code{get_pred} returns a prediction vector calculated with Data matrix and model coefficient matrix.
#' @param C Model coefficient matrix.
#' @param X Data matrix.
#' @export
get_pred <- function(C, X) return(getPred(C, X))

#' R wrapper function of getNEDist in ands_common.cpp
#'
#' \code{get_ned} returns the NED (Normal Euclidean Distance) of the matrix X.
#' @param X Data matrix.
#' @export
get_ned <- function(X) {
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    return(getNEDist(X))
}

#' R wrapper function of getLEDist in ands_common.cpp
#'
#' \code{get_lmd} returns the LMD (Local Manifold Distance) of the matrix X.
#' @param X Data matrix.
#' @param theta A matrix, in which elements capture the nonlinearity of each column in X.
#' @export
get_lmd <- function(X, theta) {
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    if (!is.matrix(theta)) stop("Inappropriate type is used in theta!")
    if (ncol(X) != ncol(theta)) stop("Dimensions of X & theta are not matched!")

    return(getLMDist(X, theta))
}

#' r wrapper function of getoptht in ands_ntsa.cpp
#'
#' \code{get_trunc} returns the truncated rank of SVD with optimal hard threshold.
#' @param sv Numeric vectors for singular values.
#' @param m Integer, size of measurements.
#' @param n Integer, data size
#' @export
get_trunc <- function(sv, m, n) {
    getOptHT(sv, m, n)
}

#' Accessor function of list data
#'
#' @param lst A list to be accessed.
#' @param name A strings specifying the data to be accessed.
#' @export
pullist <- function(lst, name) lst[[name]]

#' Function to set column names within pipes
#'
#' @param x tibble, data.frame or matrix to be renamed
#' @param nm Strings vector for the column name
#' @export
set_cnames <- function(x, nm) {
    if (ncol(x) != length(nm)) stop("Dimension of x and length of nm are not matched!")
    colnames(x) <- nm
    return(x)
}

#' Shift vector elements forward or backward depending on the lag value
#'
#' \code{shift} shifts an input vector x bacward or forward.
#' @param x An input vector.
#' @param lag Integer to set the lag value for bwd/fwd shifts.
#' @export
shift <- function(x, lag) do.call(if_else(lag < 0, "lag", "lead"), args = list(x, abs(lag)));

#' Normalization function
#'
#' @param x Numeric vector to be normalized.
#' @param method Strings, which set the normalization method.
#'      if method == "0-1", x is normalized to range between [0, 1];
#'      otherwise, x is normalized with mean 0 and unit variance.
#' @export
normalize <- function(x, method = "0-1") {
    if (method == "0-1") {
        x_ <- (x - min(x, ma.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    } else {
        x_ <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    }

    return(x_)
}

#' ODE solver with 4-th order Runge-Kutta approximation
#'
#' @param x_ini A numeric vector, set the initial states.
#' @param fun A functional object specifying the type of differential equation.
#' @param params A list containing required parameters in 'fun'.
#' @param t_end A double specifying the end time at the termination of calculation.
#' @param h A double, specifying integration time in RK-approximation.
#' @export
odeRK4 <- function(x_ini, fun, params, input = NULL, t_end, h = 0.01) {
    x <- x_ini
    t_seq <- seq(0, t_end, h)
    if(is.null(input)) input <- matrix(0, nrow = length(t_seq), ncol = length(x))

    foreach(i = icount(length(t_seq) - 1), .combine = rbind, .init = x_ini) %do% {
        k1 <- h * fun(x = x, params = c(params, u = list(input[i, ])))
        k2 <- h * fun(x = x + k1 / 2, params = c(params, u = list(input[i, ])))
        k3 <- h * fun(x = x + k2 / 2, params = c(params, u = list(input[i, ])))
        k4 <- h * fun(x = x + k3, params = c(params, u = list(input[i, ])))

        x <- x + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    } |> as_tibbler(str_c("x", 1:length(x))) |> mutate(time = t_seq, .before = everything())
}

#' Calculate discrete-time time series data
#'
#' @inheritParams odeRK4
#' @param o_noise A numeric, observation error strength/
#' @export
gen_map <- function(x_ini, fun, params, input = NULL, t_end, o_noise = 0.00) {
    x <- x_ini
    if(is.null(input)) input <- matrix(0, nrow = t_end, ncol = length(x))

    foreach(i = icount(t_end), .combine = rbind, .init = x) %do% {
        if(any(!is.finite(x))) stop("Time-series data generation fault!")
        x <- fun(x = x, params = c(params, u = list(input[i, ])))
    } |> as_tibbler(str_c("x", 1:length(x))) |> mutate(across(starts_with("x"), \(.x) .x + rnorm(length(.x), 0, o_noise * sd(.x))), time = 1:length(x1), .before = everything())
}

