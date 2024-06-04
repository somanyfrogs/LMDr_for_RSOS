#' @title Figure_3.R
#' @description Making figure 3 & 4 for the Paper entitled:
#'      "Local-manifold-distance-based regression: An estimation method for quantifying dynamic biological interactions"
#'      Initially written on 20230822.
#'      Last update: 20240531.

## Load R functions
source("R/functions.R")

## Figure 2: Comparison in forecasting performance between LMDr and S-map variants
## S-map variants:
## Normal S-map (Deyle et al. 2016, PRSE)
## Regularized S-map (Cenci et al. 2019, MEE)
## MDR S-map (Chang et al. 2021, Ecol Lett)
## Geodesic-distance S-map (Qu et al. 2023, Ecol Indic)
seed <- 188
ths <- 4
cores <- floor(detectCores() / ths)
noises <- expand_grid(iter = 1:2, p_noise = seq(0.0, 0.1, length.out = 3) + 1e-3, o_noise = seq(0.0, 0.1, length.out = 3) + 1e-3)

## One-prey-two-predator system with predator-dependent functional response (discrete-time system with three variables)
fun <- predatorFR; funJ <- predatorFRJ; funName <- "predatorFR"
size <- 3
len <- 200

## Set the model parameters
r1 <- 3.75; r2 <- 3.61; r3 <- 3.65
k <- 0.1; e <- 0.25
a12 <- 0.4; a13 <- 0.1
g12 <- 0.2; g13 <- 0.2

params <- list(r1 = r1, r2 = r2, r3 = r3, k = k, e = e, a12 = a12, a13 = a13, g12 = g12, g13 = g13)

## Generate original time series
ts_orig <- gen_map(x_ini = rep(0.1, size), fun = fun, params = params, t_end = 300)
sd_orig <- ts_orig |> summarize(across(starts_with("x"), sd))
x_ini <- last(ts_orig) |> select(!time) |> as.numeric()

## Perform simulation
cl <- makeCluster(spec = ths, type = "PSOCK", outfile = "")
clusterSetRNGStream(cl, seed)
registerDoParallel(cl)

system.time(foreach(tbl = iter(noises, by = "row"), .combine = bind_rows,
                    .packages = c("andsr", "doParallel", "foreach", "iterators", "tidyverse")) %dopar% {
    ## Set the noise parametrs
    iter <- tbl$iter
    p_noise <- tbl$p_noise
    o_noise <- tbl$o_noise
    op <- NULL

    while(is.null(op)) {
        cat(str_c("iterNo = ", iter, "; p_noise = ", p_noise, "; o_noise = ", o_noise, "\n"))

        try({
            ## Make time-series with process noise
            input <- sd_orig |> reframe(across(everything(), \(.x) rnorm(len + 1, 0, p_noise * .x))) |> as.matrix()
            ts <- gen_map(x_ini = x_ini, fun = fun, params = params, input = input, t_end = len)

            ## Calculate original Jacobian matrices
            H_true <- ts |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size)) |> gen_lstJ(fun = funJ, params = params)
            X <- ts |> mutate(across(starts_with("x"), \(.x) .x + rnorm(length(.x), 0, o_noise * sd(.x)))) |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size))
            dmat <- X |> get_ned()

            op <- foreach(i = 1:size, .combine = rbind) %do% {
                J_true <- foreach(j = 1:length(H_true), .combine = rbind) %do% {H_true[[j]][i, ]} |>
                    as_tibbler(str_c("C", 1:size)) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "True")

                ned <- smap_mulvar(X = X, col = i, dmat = dmat, range = seq(0, 25, 0.1), threadNo = cores) |> select(!theta) |> mutate(dist = "Euclidean", .before = everything())
                gdd <- smap_geodes(X = X, col = i, dmat = dmat, seed = seed, threadNo = cores, iter = 250, k_range = 3:20, sigmas = c(0.05, 0.05)) |> select(!c(k, delta, theta)) |> mutate(dist = "Geodesic", .before = everything())
                mvd <- get_mvd(X, cols = 1:size, tar = i, E = 3, lmax = 3, threadNo = cores) |> smap_mulvar(X = X, col = i, dmat = _, range = seq(0, 25, 0.1), threadNo = cores) |> select(!theta) |> mutate(dist = "Multiview", .before = everything())
                lmd <- find_best_theta(X = X, col = i, gSeed = 123, threadNo = cores, method = 1) |> select(!theta) |> mutate(dist = "LMD", .before = everything())
                comb <- bind_rows(ned, gdd, mvd, lmd)

                coef <- foreach(d = c("Euclidean", "Geodesic", "Multiview", "LMD"), .combine = rbind) %do% {
                    J_est <- comb |> filter(dist == d) |> pull(coef) |> pullist(1) |> select(!C0) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "Estimate")
                    left_join(x = J_true, y = J_est, by = join_by(seq, Jacobian)) |> mutate(dist = d, .before = everything())
                } |> group_by(dist, Jacobian) |> summarize(rho_coef = get_rho(True, Estimate), mae_coef = get_mae(True, Estimate), rmse_coef = get_rmse(True, Estimate), .groups = "drop")

                comb |> select(!c(output, coef, lambda)) |> left_join(coef, by = join_by(dist)) |> mutate(tar = str_c("x", i), .before = Jacobian)
            } |> (\(.x) bind_rows(op, .x))()
        })
    }

    op |> mutate(system = funName, p_noise = p_noise, o_noise = o_noise, .before = everything())
} |> write_csvr("output/sim_1_1.csv"))

stopCluster(cl)

## Lorenz63 system (continuous-time system with three variables)
fun <- lorenz; funJ <- lorenzJ; funName <- "Lorenz63"
size <- 3
h <- 0.01
freq <- 5
len <- 300

## Set the model parameters
sigma <- 10
beta <- 8 / 3
rho <- 28

params <- list(sigma = sigma, beta = beta, rho = rho)

## Generate original time series
ts_orig <- odeRK4(x_ini = rep(0.1, size), fun = fun, params = params, t_end = 300 * h * freq, h = h) |> slice(seq(1, length(x1), freq))
sd_orig <- ts_orig |> summarize(across(starts_with("x"), sd))
x_ini <- last(ts_orig) |> select(!time) |> as.numeric()

## Perform simulation
cl <- makeCluster(spec = ths, type = "PSOCK", outfile = "")
clusterSetRNGStream(cl, seed)
registerDoParallel(cl)

system.time(foreach(tbl = iter(noises, by = "row"), .combine = bind_rows,
                    .packages = c("andsr", "doParallel", "foreach", "iterators", "tidyverse")) %dopar% {
    ## Set the noise parameters
    iter <- tbl$iter
    p_noise <- tbl$p_noise
    o_noise <- tbl$o_noise
    op <- NULL

    system.time(while(is.null(op)) {
        cat(str_c("iterNo = ", iter, "; p_noise = ", p_noise, "; o_noise = ", o_noise, "\n"))

        try({
            ## Make time-series with process noise
            input <- sd_orig |> reframe(across(everything(), \(.x) rnorm(len * freq + 1, 0, p_noise * .x))) |>
                mutate(time = 1:length(x1), across(starts_with("x"), \(.x) ifelse(time %% freq == 1, .x, 0))) |> select(!time) |> as.matrix()
            ts <- odeRK4(x_ini = x_ini, fun = fun, params = params, input = input, t_end = len * h * freq, h = h) |> slice(seq(1, length(x1), freq))

            ## Calculate original Jacobian matrices
            H_true <- ts |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size)) |> gen_lstJ(fun = funJ, params = params, h = h * freq)
            X <- ts |> mutate(across(starts_with("x"), \(.x) .x + rnorm(length(.x), 0, o_noise * sd(.x)))) |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size))
            dmat <- X |> get_ned()

            op <- foreach(i = 1:size, .combine = rbind) %do% {
                J_true <- foreach(j = 1:length(H_true), .combine = rbind) %do% {H_true[[j]][i, ]} |>
                    as_tibbler(str_c("C", 1:size)) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "True")

                ned <- smap_mulvar(X = X, col = i, dmat = dmat, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Euclidean", .before = everything())
                gdd <- smap_geodes(X = X, col = i, dmat = dmat, seed = seed, iter = 250, k_range = 3:20, sigmas = c(0.05, 0.05)) |> select(!c(k, delta, theta)) |> mutate(dist = "Geodesic", .before = everything())
                mvd <- get_mvd(X, cols = 1:size, tar = i, E = 3, lmax = 3) |> smap_mulvar(X = X, col = i, dmat = _, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Multiview", .before = everything())
                lmd <- find_best_theta(X = X, col = i, gSeed = seed, method = 1) |> select(!theta) |> mutate(dist = "LMD", .before = everything())
                comb <- bind_rows(ned, gdd, mvd, lmd)

                coef <- foreach(d = c("Euclidean", "Geodesic", "Multiview", "LMD"), .combine = rbind) %do% {
                    J_est <- comb |> filter(dist == d) |> pull(coef) |> pullist(1) |> select(!C0) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "Estimate")
                    left_join(x = J_true, y = J_est, by = join_by(seq, Jacobian)) |> mutate(dist = d, .before = everything())
                } |> group_by(dist, Jacobian) |> summarize(rho_coef = get_rho(True, Estimate), mae_coef = get_mae(True, Estimate), rmse_coef = get_rmse(True, Estimate), .groups = "drop")

                comb |> select(!c(output, coef, lambda)) |> left_join(coef, by = join_by(dist)) |> mutate(tar = str_c("x", i), .before = Jacobian)
            } |> (\(.x) bind_rows(op, .x))()
        })
    })

    op |> mutate(system = funName, p_noise = p_noise, o_noise = o_noise, .before = everything())
} |> write_csvr("output/sim_1_2.csv"))

stopCluster(cl)

## A five species coupled Logistic map (discrete-time system with five variables, from Hsieh et al. 2008, Am Nat)
fun <- glv; funJ <- glvJ; funName <- "5CLM"
size <- 5
len <- 200

## Set the model aprameters
set.seed(seed)
r <- rnorm(size, 3.65, 0.05)
A <- matrix(rnorm(size^2, 0.00, 0.05), nrow = size) |> (\(.x) .x - diag(diag(.x) + r))()

params <- list(r = r, A = A)

## Generate original time series
ts_orig <- gen_map(x_ini = rep(0.1, size), fun = fun, params = params, t_end = 300)
sd_orig <- ts_orig |> summarize(across(starts_with("x"), sd))
x_ini <- last(ts_orig) |> select(!time) |> as.numeric()

## Perform simulation
cl <- makeCluster(spec = ths, type = "PSOCK", outfile = "")
clusterSetRNGStream(cl, seed)
registerDoParallel(cl)

system.time(foreach(tbl = iter(noises, by = "row"), .combine = bind_rows,
                    .packages = c("andsr", "doParallel", "foreach", "iterators", "tidyverse")) %dopar% {
    ## Set the noise parameters
    iter <- tbl$iter
    p_noise <- tbl$p_noise
    o_noise <- tbl$o_noise
    op <- NULL

    while(is.null(op)) {
        cat(str_c("iterNo = ", iter, "; p_noise = ", p_noise, "; o_noise = ", o_noise, "\n"))

        try({
            ## Make time series with process noise
            input <- sd_orig |> reframe(across(everything(), \(.x) rnorm(len + 1, 0, p_noise * .x))) |> as.matrix()
            ts <- gen_map(x_ini = x_ini, fun = fun, params = params, input = input, t_end = len)

            ## Calculate original Jacobian matrices
            H_true <- ts |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size)) |> gen_lstJ(fun = funJ, params = params)
            X <- ts |> mutate(across(starts_with("x"), \(.x) .x + rnorm(length(.x), 0, o_noise * sd(.x)))) |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size))
            dmat <- X |> get_ned()

            op <- foreach(i = 1:size, .combine = rbind) %do% {
                J_true <- foreach(j = 1:length(H_true), .combine = rbind) %do% {H_true[[j]][i, ]} |>
                    as_tibbler(str_c("C", 1:size)) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "True")

                ned <- smap_mulvar(X = X, col = i, dmat = dmat, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Euclidean", .before = everything())
                gdd <- smap_geodes(X = X, col = i, dmat = dmat, seed = seed, iter = 250, k_range = 3:20, sigmas = c(0.05, 0.05)) |> select(!c(k, delta, theta)) |> mutate(dist = "Geodesic", .before = everything())
                mvd <- get_mvd(X, cols = 1:size, tar = i, E = 3, lmax = 3) |> smap_mulvar(X = X, col = i, dmat = _, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Multiview", .before = everything())
                lmd <- find_best_theta(X = X, col = i, gSeed = seed, method = 1) |> select(!theta) |> mutate(dist = "LMD", .before = everything())
                comb <- bind_rows(ned, gdd, mvd, lmd)

                coef <- foreach(d = c("Euclidean", "Geodesic", "Multiview", "LMD"), .combine = rbind) %do% {
                    J_est <- comb |> filter(dist == d) |> pull(coef) |> pullist(1) |> select(!C0) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "Estimate")
                    left_join(x = J_true, y = J_est, by = join_by(seq, Jacobian)) |> mutate(dist = d, .before = everything())
                } |> group_by(dist, Jacobian) |> summarize(rho_coef = get_rho(True, Estimate), mae_coef = get_mae(True, Estimate), rmse_coef = get_rmse(True, Estimate), .groups = "drop")

                comb |> select(!c(output, coef, lambda)) |> left_join(coef, by = join_by(dist)) |> mutate(tar = str_c("x", i), .before = Jacobian)
            }
        }) |> (\(.x) bind_rows(op, .x))()
    }

    op |> mutate(system = funName, p_noise = p_noise, o_noise = o_noise, .before = everything())
} |> write_csvr("output/sim_1_3_2.csv"))

## A five species coupled food chain model (Post et al. 2000, Ecology)
fun <- cfc; funJ <- cfcJ; funName <- "5CFC"
size <- 5
h <- 0.01
freq <- 200
len <- 300

## Set the model parameters
v1 <- 0.10; v2 <- 0.07
lambda1 <- 3.2; lambda2 <- 2.9
mu1 <- 0.15; mu2 <- 0.15
kappa1 <- 2.5; kappa2 <- 2.0
k <- 1.2

params <- list(A = rbind(c(-1 / k, -mu1 * kappa1, -mu2 * kappa2, 0, 0),
                         c(mu1 * kappa1, 0, 0, -v1 * lambda1, 0),
                         c(mu2 * kappa2, 0, 0, 0, -v2 * lambda2),
                         c(0, v1 * lambda1, 0, 0, 0),
                         c(0, 0, v2 * lambda2, 0, 0)),
               r = c(1, -mu1, -mu2, -v1, -v2),
               Rstar = 0.3, Cstar = c(0.5, 0.5))

## Generate original time series
ts_orig <- odeRK4(x_ini = rep(0.1, size), fun = fun, params = params, t_end = 300 * h * freq, h = h) |> slice(seq(1, length(x1), freq))
sd_orig <- ts_orig |> summarize(across(starts_with("x"), sd))
x_ini <- last(ts_orig) |> select(!time) |> as.numeric()

## Perform simulation
cl <- makeCluster(spec = ths, type = "PSOCK", outfile = "")
clusterSetRNGStream(cl, seed)
registerDoParallel(cl)

system.time(foreach(tbl = iter(noises, by = "row"), .combine = bind_rows,
                    .packages = c("andsr", "doParallel", "foreach", "iterators", "tidyverse")) %dopar% {
    ## Set the noise parameters
    iter <- tbl$iter
    p_noise <- tbl$p_noise
    o_noise <- tbl$o_noise
    op <- NULL

    system.time(while(is.null(op)) {
        cat(str_c("iterNo = ", iter, "; p_noise = ", p_noise, "; o_noise = ", o_noise, "\n"))

        try({
            ## Make time-series with process noise
            input <- sd_orig |> reframe(across(everything(), \(.x) rnorm(len * freq + 1, 0, p_noise * .x))) |>
                mutate(time = 1:length(x1), across(starts_with("x"), \(.x) ifelse(time %% freq == 1, .x, 0))) |> select(!time) |> as.matrix()
            ts <- odeRK4(x_ini = x_ini, fun = fun, params = params, input = input, t_end = len * h * freq, h = h) |> slice(seq(1, length(x1), freq))

            ## Calculate original Jacobian matrices
            H_true <- ts |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size)) |> gen_lstJ(fun = funJ, params = params, h = h * freq)
            X <- ts |> mutate(across(starts_with("x"), \(.x) .x + rnorm(length(.x), 0, o_noise * sd(.x)))) |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size))
            dmat <- X |> get_ned()

            op <- foreach(i = 1:size, .combine = rbind) %do% {
                J_true <- foreach(j = 1:length(H_true), .combine = rbind) %do% {H_true[[j]][i, ]} |>
                    as_tibbler(str_c("C", 1:size)) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "True")

                ned <- smap_mulvar(X = X, col = i, dmat = dmat, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Euclidean", .before = everything())
                gdd <- smap_geodes(X = X, col = i, dmat = dmat, seed = seed, iter = 250, k_range = 3:20, sigmas = c(0.05, 0.05)) |> select(!c(k, delta, theta)) |> mutate(dist = "Geodesic", .before = everything())
                mvd <- get_mvd(X, cols = 1:size, tar = i, E = 3, lmax = 3) |> smap_mulvar(X = X, col = i, dmat = _, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Multiview", .before = everything())
                lmd <- find_best_theta(X = X, col = i, gSeed = seed, method = 1) |> select(!theta) |> mutate(dist = "LMD", .before = everything())
                comb <- bind_rows(ned, gdd, mvd, lmd)

                coef <- foreach(d = c("Euclidean", "Geodesic", "Multiview", "LMD"), .combine = rbind) %do% {
                    J_est <- comb |> filter(dist == d) |> pull(coef) |> pullist(1) |> select(!C0) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "Estimate")
                    left_join(x = J_true, y = J_est, by = join_by(seq, Jacobian)) |> mutate(dist = d, .before = everything())
                } |> group_by(dist, Jacobian) |> summarize(rho_coef = get_rho(True, Estimate), mae_coef = get_mae(True, Estimate), rmse_coef = get_rmse(True, Estimate), .groups = "drop")

                comb |> select(!c(output, coef, lambda)) |> left_join(coef, by = join_by(dist)) |> mutate(tar = str_c("x", i), .before = Jacobian)
            } |> (\(.x) bind_rows(op, .x))()
        })
    })

    op |> mutate(system = funName, p_noise = p_noise, o_noise = o_noise, .before = everything())
} |> write_csvr("output/sim_1_4.csv"))

stopCluster(cl)

## 
fun <- dlv; funJ <- dlvJ; funName <- "MulInt"
size <- 20
len <- 200

## Set the model parameters
set.seed(seed)
r <- rnorm(size, 3.0, 0.05)
A <- matrix(rnorm(size^2, 0, 0.05), nrow = size, ncol = size) |> (\(.x) .x - r * diag(size))()

params <- list(r = r, A = A)

## Generate original time series
ts_orig <- gen_map(x_ini = rep(0.1, size), fun = fun, params = params, t_end = 300)
sd_orig <- ts_orig |> summarize(across(starts_with("x"), sd))
x_ini <- last(ts_orig) |> select(!time) |> as.numeric()

## Perform simulation
cl <- makeCluster(spec = ths, type = "PSOCK", outfile = "")
clusterSetRNGStream(cl, seed)
registerDoParallel(cl)

system.time(foreach(tbl = iter(noises, by = "row"), .combine = bind_rows,
                    .packages = c("andsr", "doParallel", "foreach", "iterators", "tidyverse")) %dopar% {
    ## Set the noise parameters
    iter <- tbl$iter
    p_noise <- tbl$p_noise
    o_noise <- tbl$o_noise
    op <- NULL

    system.time(while(is.null(op)) {
        cat(str_c("iterNo = ", iter, "; p_noise = ", p_noise, "; o_noise = ", o_noise, "\n"))

        try({
            ## Make time-series with process noise
            input <- sd_orig |> reframe(across(everything(), \(.x) rnorm(len + 1, 0, p_noise * .x))) |> as.matrix()
            ts <- gen_map(x_ini = x_ini, fun = fun, params = params, input = input, t_end = len)

            ## Calculate original Jacobian matrices
            H_true <- ts |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size)) |> gen_lstJ(fun = funJ, params = c(params, u = 0))
            X <- ts |> mutate(across(starts_with("x"), \(.x) .x + rnorm(length(.x), 0, o_noise * sd(.x)))) |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size))
            dmat <- X |> get_ned()

            op <- foreach(i = 1:5, .combine = rbind) %do% {
                J_true <- foreach(j = 1:length(H_true), .combine = rbind) %do% {H_true[[j]][i, ]} |>
                    as_tibbler(str_c("C", 1:size)) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "True")

                ned <- smap_mulvar(X = X, col = i, dmat = dmat, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Euclidean", .before = everything())
                gdd <- smap_geodes(X = X, col = i, dmat = dmat, seed = seed, iter = 250, k_range = 3:20, sigmas = c(0.05, 0.05)) |> select(!c(k, delta, theta)) |> mutate(dist = "Geodesic", .before = everything())
                mvd <- get_mvd(X, cols = 1:size, tar = i, E = 3, lmax = 2) |> smap_mulvar(X = X, col = i, dmat = _, range = seq(0, 25, 0.1)) |> select(!theta) |> mutate(dist = "Multiview", .before = everything())
                lmd <- find_best_theta(X = X, col = i, gSeed = seed, method = 1) |> select(!theta) |> mutate(dist = "LMD", .before = everything())
                comb <- bind_rows(ned, gdd, mvd, lmd)

                coef <- foreach(d = c("Euclidean", "Geodesic", "Multiview", "LMD"), .combine = rbind) %do% {
                    J_est <- comb |> filter(dist == d) |> pull(coef) |> pullist(1) |> select(!C0) |> mutate(seq = 1:length(C1)) |> pivot_longer(!seq, names_to = "Jacobian", values_to = "Estimate")
                    left_join(x = J_true, y = J_est, by = join_by(seq, Jacobian)) |> mutate(dist = d, .before = everything())
                } |> group_by(dist, Jacobian) |> summarize(rho_coef = get_rho(True, Estimate), mae_coef = get_mae(True, Estimate), rmse_coef = get_rmse(True, Estimate), .groups = "drop")

                comb |> select(!c(output, coef, lambda)) |> left_join(coef, by = join_by(dist)) |> mutate(tar = str_c("x", i), .before = Jacobian)
            } |> (\(.x) bind_rows(op, .x))()
        })
    })

    op |> mutate(system = funName, p_noise = p_noise, o_noise = o_noise, .before = everything())
} |> write_csvr("output/sim_1_5.csv"))

stopCluster(cl)

## Make figure 3-4
sim_1_1 <- read_csvr("output/sim_1_1.csv")
sim_1_2 <- read_csvr("output/sim_1_2.csv") |> filter(rho > 0.975)
sim_1_3 <- read_csvr("output/sim_1_3.csv")
sim_1_4 <- read_csvr("output/sim_1_4.csv")
sim_1_5 <- read_csvr("output/sim_1_5.csv")

sim <- bind_rows(sim_1_1, sim_1_2, sim_1_3, sim_1_4, sim_1_5) |> mutate(system = factor(system, levels = unique(system)), dist = factor(dist, levels = c("LMD", "Multiview", "Geodesic", "Euclidean")))

gp1 <- sim |> ggplot(aes(x = factor(o_noise), y = rho_coef, fill = dist)) + facet_grid(system ~ p_noise, scales = "free_y") + geom_boxplot(outlier.size = 0.1, linewidth = 0.25) +
    scale_fill_manual(values = ggsci::pal_jama()(4)[4:1]) + xlab(expression(paste("Observation error ", (sigma[obs])))) + ylab(expression(paste("Correlation coefficient ", rho))) + theme_st(lunit = 2)

gp2 <-sim |>  ggplot(aes(x = factor(o_noise), y = rho, fill = dist)) + facet_grid(system ~ p_noise, scales = "free_y") + geom_boxplot(outlier.size = 0.1, linewidth = 0.25) +
   scale_fill_manual(values = ggsci::pal_jama()(4)[4:1]) + xlab(expression(paste("Observation error ", (sigma[obs])))) + ylab(expression(paste("Correlation coefficient ", rho))) + theme_st(lunit = 2)

gp1 |> ggsaver("fig03", width = 14.0, height = 11.0, ext = "pdf")
gp2 |> ggsaver("fig04", width = 14.0, height = 11.0, ext = "pdf")
