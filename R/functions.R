#' @title functions.R
#' @description An R file describing common functions for the paper entitled:
#'      "Local Manifold Distance regression: an estimation method for fluctuating biological interactions with empirical time series"
#'      Initialy written on 20230822.
#'      Last update: 20240530.

## Load dependent packages
library(andsr)
library(doParallel)
library(foreach)
library(ggraph)
library(igraph)
library(patchwork)
library(tidygraph)
library(tidyverse)
library(viridis)

getwd() |> setwd()

#' Calculate the dissimilarity metrics of two matrices column by column
#'
#' @param X A numeric matrix.
#' @param Y A numeric matrix.
compute_stats <- function(X, Y) {
    if(all(dim(X) != dim(Y))) stop("Different dimensions b/w X & Y")

    tibble(rho = colMeans(apply(X, 2, scale) * apply(Y, 2, scale), na.rm = TRUE),
           mae = colMeans(abs(X - Y), na.rm = TRUE),
           rmse = sqrt(colMeans((X - Y)^2, na.rm = TRUE)))
}

#' 5-species coupled food chain model from Polis et al. (2000).
#'
#' @param x A state vector in turn consisting of Resource, Consumer 1, 2 and Predator 1, 2.
#' @param params A list containing model parameters.
cfc <- function(x, params) with(params, {
    H1 <- A * rbind(
                c(1, 1 / (x[1] + Rstar), 1 / (x[1] + Rstar), 0, 0),
                c(1 / (x[1] + Rstar), 0, 0, 1 / (x[2] + Cstar[1]), 0),
                c(1 / (x[1] + Rstar), 0, 0, 0, 1 / (x[3] + Cstar[2])),
                c(0, 1 / (x[2] + Cstar[1]), 0, 0, 0),
                c(0, 0, 1 / (x[3] + Cstar[2]), 0, 0))

    as.numeric((r + H1 %*% x + u) * x)
})

#' Discretized Lotka-Volterra (DLV) model
#'
#' A: An interaction matrix; r: An intrinsic growth rate
#' @inheritParams cfc
dlv <- function(x, params) with(params, as.numeric(exp(r + A %*% x + u) * x))

#' Generalized Lotka-Volterra (GLV) model
#'
#' A: An interaction matrix; r: An intrinsic growth rate
#' @inheritParams cfc
glv <- function(x, params) with(params, as.numeric((r + A %*% x + u) * x))

#' Henon map
#'
#' a: Henon parameter 1.
#' b: Henon parameter 2.
#' @inheritParams cfc
henon <- function(x, params) with(params, c(1 - a * x[1]^2 + x[2] + u[1], b * x[1] + u[2]))

#' Lorenz63 system
#'
#' sigma: Lorenz63 parameter 1.
#' rho: Lorenz63 parameter 2.
#' beta: Lorenz63 parameter 3.
#' @inheritParams cfc
lorenz <- function(x, params) with(params, {
    M <- matrix(c(-sigma, rho, x[2], sigma, -1, 0, 0, -x[1], -beta), nrow = 3)
    as.numeric(M %*% x + u)
})

#' Calculate true Jacobian matrix of 5-species coupled food chain model.
#'
#' @inheritParams cfc
cfcJ <- function(x, params) with(params, {
    H1 <- A * rbind(
                c(1, 1 / (x[1] + Rstar), 1 / (x[1] + Rstar), 0, 0),
                c(1 / (x[1] + Rstar), 0, 0, 1 / (x[2] + Cstar[1]), 0),
                c(1 / (x[1] + Rstar), 0, 0, 0, 1 / (x[3] + Cstar[2])),
                c(0, 1 / (x[2] + Cstar[1]), 0, 0, 0),
                c(0, 0, 1 / (x[3] + Cstar[2]), 0, 0))

    H2 <- H1 + rbind(c(-sum(H1[1, 2:3]) / (x[1] + Rstar), 0, 0, 0, 0),
                     c(-H1[2, 1] * x[1] / (x[1] + Rstar), -H1[2, 4] * x[4] / (x[2] + Cstar[1]), 0, 0, 0),
                     c(-H1[3, 1] * x[1] / (x[1] + Rstar), 0, -H1[3, 5] * x[5] / (x[3] + Cstar[2]), 0, 0),
                     c(0, -H1[4, 2] * x[2] / (x[2] + Cstar[1]), 0, 0, 0),
                     c(0, 0, -H1[5, 3] * x[3] / (x[3] + Cstar[2]), 0, 0))

    diag(c(r + H1 %*% x)) + H2 * x
})

#' Calculate true Jacobian matrix of DLV dynamics
#'
#' @inheritParams cfc
dlvJ <- function(x, params) with(params, (A + diag(1 / x)) * dlv(x, params))

#' Calculate true Jacobian matrix of GLV dynamics
#'
#' @inheritParams cfc
glvJ <- function(x, params) with(params, diag(c(r + A %*% x)) + A * x)

#' Calculate true Jacobian matrix of Henon map
#'
#' @inheritParams cfc
henonJ <- function(x, params) with(params, matrix(c(-2 * a * x[1], b, 1, 0), nrow = 2))

#' Calculate true Jacobian matrix of Lorenz63 system
#'
#' @inheritParams cfc
lorenzJ <- function(x, params) with(params, matrix(c(-sigma, rho - x[3], x[2], sigma, -1, x[1], 0, -x[1], -beta), nrow = 3))

#' Generate a set of true Jacobians for a discrete map or discretized flow
#'
#' @param X An embedding matrix.
#' @param fun function object to be used for data generation.
#' @param params A list of model parameters.
#' @param h A double, specifying integration time.
gen_lstJ <- function(X, fun, params, h = NULL) {
    H <- list()
    
    for(i in 1:nrow(X)) {
        J <- fun(x = X[i, ], params = params)
        H <- c(H, list(J))
    }

    if(!is.null(h)) H <- lapply(H, \(.x) matExp(h * .x))
    return(H)
}

#' Calculate LMD (Local Manifold Distance) from the target point
#'
#' @param M An embedding matrix.
#' @param Theta A locality matrix.
#' @param tar A numeric vector of the target point.
LMDist <- function(M, Theta, tar) {
    tmp <- foreach(i = 1:nrow(M), .combine = rbind) %do% c(M[i, ] - tar)
    rowSums((tmp %*% Theta)^2) |> sqrt() |> as.numeric()
}

#' Calculate matrix exponential.
#'
#' @param M A real/complex matrix.
matExp <- function(M) {
    eig <- eigen(M)
    P <- eig$vectors
    tmp <- P %*% diag(exp(eig$values)) %*% solve(P)
    if(!is.complex(M)) tmp <- Re(tmp)
    return(tmp)
}

#' One-prey-two-predator system with predator-dependent functional response (Beddington 1975; DeAngelis et al. 1975)
#'
#' @inheritParams cfc
predatorFR <- function(x, params) with(params, {
    denom12 <- 1 + k * a12 * x[1] + g12 * x[2]
    denom13 <- 1 + k * a13 * x[1] + g13 * x[3]
    x1 <- x[1] * (r1 * (1 - x[1]) - a12 * x[2] / denom12 - a13 * x[3] / denom13)
    x2 <- x[2] * (r2 * (1 - x[2]) + e * a12 * x[1] / denom12)
    x3 <- x[3] * (r3 * (1 - x[3]) + e * a13 * x[1] / denom13)
    
    x <- c(x1, x2, x3)
})

#' Calculate true Jacobian matrix of predator FR dynamics
#'
#' @inheritParams cfc
predatorFRJ <- function(x, params) with(params, {
    denom12 <- 1 + k * a12 * x[1] + g12 * x[2]
    denom13 <- 1 + k * a13 * x[1] + g13 * x[3]

    rbind(c(r1 * (1 - x[1]) - a12 * x[2] / denom12 - a13 * x[3] / denom13 + x[1] * (-r1 + k * a12^2 * x[2] / denom12^2 + k * a13^2 * x[3] / denom13^2), a12 * x[1] / denom12 * (g12 * x[2] / denom12 - 1), a13 * x[1] / denom13 * (g13 * x[3] / denom13 - 1)),
          c(e * a12 * x[2] / denom12 * (1 - k * a12 * x[1] / denom12), r2 * (1 - x[2]) + e * a12 * x[1] / denom12 - x[2] * (r2 + e * a12 * g12 * x[1] / denom12^2), 0),
          c(e * a13 * x[3] / denom13 * (1 - k * a13 * x[1] / denom13), 0, r3 * (1 - x[3]) + e * a13 * x[1] / denom13 - x[3] * (r3 + e * a13 * g13 * x[1] / denom13^2)))
})

#' Modified version of read_csv
#'
#' @param path A strings, setting the file path to be loaded from.
read_csvr <- function(path, ...) read_csv(path, progress = FALSE, show_col_types = FALSE, ...)

#' Modified version of write_csv
#'
#' @param x A tibble or data.frame to be written to.
#' @param path A strings, setting the file path to be written to.
write_csvr <- function(x, path) write_csv(x, path, progress = FALSE)

#' Modified version of ggsave
#'
#' @param gp ggplot object.
#' @param name A strings (file name).
#' @param width A numeric, figure width with cm scale.
#' @param height A numeric, figure height with cm scale.
#' @param ext A strings, specifying the file type (default = "eps").
#' @param rs An integer, specifying the value of fallback_resolution.
ggsaver <- function(gp, name, width, height, ext = "eps", rs = 900) {
    file <- str_c("fig/", name, ".", ext)

    if(ext == "eps") {
        ggsave(file, gp, device = cairo_ps, fallback_resolution = rs, family = "Times", width = width, height = height, units = "cm")
    } else if (ext == "pdf") {
        ggsave(file, gp, device = cairo_pdf,fallback_resolution = rs, family = "Times", width = width, height = height, units = "cm")
    } else {
        ggsave(file, gp, family = "Times", width = width, height = height, units = "cm")
    }
}

#' Original ggplot theme1
#'
#' @param just A numeric vector of legend justification.
#' @param pos A numeric vector of legend position.
#' @param lunit A numeric specifying the line unit.
theme_st <- function(just = NULL, pos = NULL, lunit = 5.5) {
    mgn <- margin(0.5 * lunit, 0.5 * lunit, 0.5 * lunit, 0.5 * lunit)

    th <- theme_light() +
        theme(axis.text = element_text(size = 6),
              axis.title = element_text(size = 8),
              legend.background = element_blank(),
              legend.key.size = unit(0.25, "cm"),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 6),
              panel.grid.minor = element_blank(),
              plot.tag = element_text(size = 8),
              plot.title = element_text(size = 8),
              strip.text = element_text(size = 6, color = "black", margin = mgn))

    if(!is.null(just)) th <- th + theme(legend.justification = just)
    if(!is.null(pos)) th <- th + theme(legend.position = pos)
    return(th)
}

#' Original ggplot theme2
#'
#' @inheritParams theme_st
theme_wt0 <- function(just = NULL, pos = NULL) {
    th <- theme_st(just, pos) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.title = element_blank(),
              panel.background = element_rect(color = "black", fill = "white"),
              panel.grid = element_blank())
    return(th)
}

#' Original ggplot theme3
#'
#' @inheritParams theme_st
theme_wt1 <- function(just = NULL, pos = NULL) theme_wt0(just, pos) + theme(axis.text = element_blank(), axis.title = element_blank())

#' Original ggplot theme4
#'
theme_wt2 <- theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
                   panel.grid = element_blank(), panel.border = element_blank(), plot.margin = margin(0, 0, 0, 0))

