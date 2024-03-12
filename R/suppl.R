#' @title Figure_S1.R
#' @description Making figure S1 for the Paper entitled:
#'      "Local-Manifold-Distance-based regression: An estimation method for quantifying dynamic biological interactions"
#'      Initially written on 20240308.
#'      Last update: 20240310.

## Load R functions
source("R/functions.R")

seed <- 123
fun <- glv; funJ <- glvJ
size <- 5

lens <- expand_grid(iter = 1:100, length = seq(50, 250, 50))

set.seed(seed)
r <- rnorm(size, 3.65, 0.05)
A <- matrix(rnorm(size^2, 0.00, 0.05), nrow = size) |> (\(.x) .x - diag(diag(.x) + r))()

params <- list(r = r, A = A)

## GEnerate original time series
ts_orig <- gen_map(x_ini = rep(0.1, size), fun = fun, params = params, t_end = 300)
sd_orig <- ts_orig |> summarize(across(starts_with("x"), sd))
x_ini <- last(ts_orig) |> select(!time) |> as.numeric()

o_noise <- 0.1001

system.time(foreach(tbl = iter(lens, by = "row"), .combine = bind_rows) %do% {
    repNo <- tbl$iter
    len <- tbl$length
    
    cat(str_c("iterNo = ", repNo, "; length = ", len, "\n"))

    try({
        ## Generate time series data
        ts <- gen_map(x_ini = x_ini, fun = fun, params = params, t_end = len)
        H_true <- ts |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size)) |> gen_lstJ(fun = funJ, params = params)
        X <- ts |> mutate(across(starts_with("x"), \(.x) .x + rnorm(length(.x), 0, o_noise * sd(.x)))) |> gen_emat(cols = str_c("x", 1:size), lags = rep(0, size))

        ## Perofmr LMDr and MDR S-map
        t_lmd <- proc.time(); lmd <- find_best_theta(X = X, gSeed = 123, method = 1); t_lmd <- proc.time() - t_lmd
        t_mvd <- proc.time(); mvd <- get_mvd(X, cols = 1:size, tar = 1, E = 5, lmax = 3) |> smap_mulvar(X = X, col = 1, dmat = _, range = seq(0, 25, 0.1)); t_mvd <- proc.time() - t_mvd

        ## Compare true and estimated Jacobian coefficients
        J_true <- foreach(i = 1:length(H_true), .combine = rbind) %do% {H_true[[i]][1, ]}
        J_lmd <- lmd$coef[[1]] |> select(!C0) |> as.matrix()
        J_mvd <- mvd$coef[[1]] |> select(!C0) |> as.matrix()

        lmd_stat <- foreach(i = 1:size, .combine = rbind) %do% tibble(method = "LMD", rho_coef = get_rho(J_true[, i], J_lmd[, i]), mae_coef = get_mae(J_true[, i], J_lmd[, i]), rmse_coef = get_rmse(J_true[, i], J_lmd[, i]))
        mvd_stat <- foreach(i = 1:size, .combine = rbind) %do% tibble(method = "MVD", rho_coef = get_rho(J_true[, i], J_mvd[, i]), mae_coef = get_mae(J_true[, i], J_mvd[, i]), rmse_coef = get_rmse(J_true[, i], J_mvd[, i]))
        lmd <- lmd |> select(!c(theta, lambda, output, coef)) |> mutate(method = "LMD", len = len, repNo = repNo, time = t_lmd[3], .before = everything()) |> left_join(lmd_stat, by = join_by(method))
        mvd <- mvd |> select(!c(theta, output, coef)) |> mutate(method = "MVD", len = len, repNo = repNo, time = t_mvd[3], .before = everything()) |> left_join(mvd_stat, by = join_by(method))

        bind_rows(lmd, mvd)
    })
} |> write_csvr("output/suppl.csv"))

sim <- read_csvr("output/suppl.csv")

gp1 <- sim |> pivot_longer(c(rmse, rmse_coef), names_to = "var", values_to = "RMSE") |> mutate(var = ifelse(var == "rmse", "Prediction", "Jacobian") |> factor(levels = c("Prediction", "Jacobian"))) |>
    ggplot(aes(x = factor(len), y = RMSE, fill = method)) + facet_wrap(. ~ var, scales = "free_y") + geom_boxplot(outlier.size = 0.5) + ggsci::scale_fill_npg() +
    labs(tag = expression((italic(a)))) + xlab(expression(paste("Data length ", (italic(m))))) + ylab("RMSE") + theme_st(lunit = 2)
gp2 <- sim |> distinct(method, len, repNo, time) |> ggplot(aes(x = len, y = time, color = method)) + geom_point(position = position_jitter(width = 1), shape = 16, alpha = 0.5) +
    geom_smooth() + ggsci::scale_color_npg() + labs(tag = expression((italic(b)))) + xlab(expression(paste("Data length ", (italic(m))))) + ylab("Computational time (sec)") + theme_st(lunit = 2)

(gp1 / gp2) |> ggsaver("figS1", width = 11.0, height = 8.8, ext = "pdf")
