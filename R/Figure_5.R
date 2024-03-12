#' @title Figure_5.R
#' @description Making figure 5, 6 for the Paper entitled:
#'      "Local-manifold-distance-regression: an estimation method for fluctuating biological interactions with empirical time series"
#'      Initially written on 20230825.
#'      Last update: 20240308.

## Load R functions
source("R/functions.R")

protozoan <- read_csvr("data/protozoan.csv")
seed <- 123
cols <- c("Pa", "Dn")
iterNo <- 100
threadNo <- 4
cores <- floor(detectCores() / threadNo)

cl <- makeCluster(spec = threadNo, type = "PSOCK", outfile = "")
clusterSetRNGStream(cl, seed)
registerDoParallel(cl)

system.time(foreach(run = icountn(c(iterNo, 3)), .combine = rbind, .packages = c("andsr", "foreach", "tidyverse")) %dopar% {
    ## Generate time series data of protozoan predator-prey system 
    cat(str_c("exp = ", run[2], "; run = ", run[1], "\n"))
    ts <- protozoan |> filter(Exp == run[2]) |>
        mutate(Pa = Pa + rnorm(length(Pa), 0, 0.417 * sqrt(Pa)), Dn = Dn + rnorm(length(Dn), 0, 0.165 * sqrt(Dn)))

    ## Estimate embedding dimension for Paramecium and Didinium, then reconstructed delay embeddings and multiview distance
    Es <- ts |> find_best_dim(cols = cols, range = 2:10) |> pull(E)
    Xp <- ts |> mutate(across(all_of(cols), \(.x) andsr::normalize(.x, method = "scale"))) |> gen_emat(cols = c(cols, rep("Pa", Es[1] - 2)), lags = c(0, -(0:(Es[1] - 2))))
    Xd <- ts |> mutate(across(all_of(cols), \(.x) andsr::normalize(.x, method = "scale"))) |> gen_emat(cols = c(cols, rep("Dn", Es[2] - 2)), lags = c(0, -(0:(Es[2] - 2))))

    ## Perform LMDr and MDR S-map analysis
    lmdr_p <- Xp |> find_best_theta(col = 1, Tp = 1, gSeed = 123, method = 1)
    lmdr_d <- Xd |> find_best_theta(col = 2, Tp = 1, gSeed = 123, method = 1)
    smap_p <- ts |> get_mvd(cols = cols, tar = 1, E = Es[1], Tp = 1, lmax = 3) |> smap_net(X = Xp, col = 1, dmat = _, range = seq(0, 20, 1))
    smap_d <- ts |> get_mvd(cols = cols, tar = 2, E = Es[2], Tp = 1, lmax = 3) |> smap_net(X = Xd, col = 2, dmat = _, range = seq(0, 20, 1))

    op_lmdr_p <- lmdr_p$coef[[1]] |> select(C1, C2) |> mutate(method = "LMDr", tar = "Pa", exp = run[2], run = run[1], rho = lmdr_p$rho, mae = lmdr_p$mae, rmse = lmdr_p$rmse, Time = ts$Time, Pa = ts$Pa, Dn = ts$Dn, .before = everything())
    op_lmdr_d <- lmdr_d$coef[[1]] |> select(C1, C2) |> mutate(method = "LMDr", tar = "Dn", exp = run[2], run = run[1], rho = lmdr_d$rho, mae = lmdr_d$mae, rmse = lmdr_d$rmse, Time = ts$Time, Pa = ts$Pa, Dn = ts$Dn, .before = everything())
    op_smap_p <- smap_p$coef[[1]] |> select(C1, C2) |> mutate(method = "MDR", tar = "Pa", exp = run[2], run = run[1], rho = smap_p$rho, mae = smap_p$mae, rmse = smap_p$rmse, Time = ts$Time, Pa = ts$Pa, Dn = ts$Dn, .before = everything())
    op_smap_d <- smap_d$coef[[1]] |> select(C1, C2) |> mutate(method = "MDR", tar = "Dn", exp = run[2], run = run[1], rho = smap_d$rho, mae = smap_d$mae, rmse = smap_d$rmse, Time = ts$Time, Pa = ts$Pa, Dn = ts$Dn, .before = everything())
    bind_rows(op_lmdr_p, op_lmdr_d, op_smap_p, op_smap_d)
} |> write_csvr("output/sim_2"))

## Make figure 
sim <- read_csvr("output/sim_2.csv") |> mutate(tar = factor(tar, levels = c("Pa", "Dn")), run = 100 * (exp - 1) + run, exp = factor(exp, labels = c("1: CC0.500", "2: CC0.375", "3: CC0.500")), SI = ifelse(tar == "Pa", C2, C1))

gp1 <- sim |> group_by(method, tar, exp, run) |> summarize(across(everything(), \(.x) mean(.x, na.rm = TRUE)), .groups = "drop") |>
    ggplot(aes(x = exp, y = rmse, color = method)) + facet_wrap(. ~ tar, scales = "free_y") + ggbeeswarm::geom_quasirandom(method = "pseudorandom", shape = 16, size = 1.0, alpha = 0.5) +
    stat_summary(aes(fill = method), geom = "point", fun = mean, shape = 23, color = "black") + ggsci::scale_color_npg() + ylab("Prediction skill (RMSE)") +
    labs(tag = expression((italic(a))), color = "method", fill = "method") + theme_st(lunit = 2) + theme(axis.title.x = element_blank())

gp2 <- sim |> group_by(method, tar, exp, run) |> summarize(across(everything(), \(.x) mean(.x, na.rm = TRUE)), .groups = "drop") |>
    ggplot(aes(x = exp, y = SI, color = method)) + facet_wrap(. ~ tar, scales = "free_y") + ggbeeswarm::geom_quasirandom(method = "pseudorandom", shape = 16, size = 1.0, alpha = 0.5) +
    stat_summary(aes(fill = method), geom = "point", fun = mean, shape = 23, color = "black") + ggsci::scale_color_npg() + ylab("Estimated Jacobian") +
    labs(tag = expression((italic(b))), color = "method", fill = "method") + theme_st(lunit = 2) + theme(axis.title.x = element_blank())

((gp1 / gp2) + plot_layout(guides = "collect")) |> ggsaver("fig05", width = 11, height = 10, ext = "pdf")

gp3 <- sim |> filter(tar == "Pa") |> ggplot(aes(x = Pa, y = SI)) + facet_wrap(. ~ method) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_smooth(aes(fill = exp, group = factor(run)), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), color = NA, alpha = 0.01) +
    geom_smooth(aes(color = exp), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), linewidth = 0.5, se = FALSE) + ggsci::scale_color_npg() + ggsci::scale_fill_npg() +
    scale_y_continuous(limits = c(-1.0, 0.5)) + xlab(expression(paste("Scaled density ", (italic(Paramecium))))) + ylab("Effect of Dn on Pa") + labs(tag = expression((italic(a)))) + theme_st(lunit = 2)

gp4 <- sim |> filter(tar == "Dn") |> ggplot(aes(x = Pa, y = SI)) + facet_wrap(. ~ method) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_smooth(aes(fill = exp, group = factor(run)), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), color = NA, alpha = 0.01) +
    geom_smooth(aes(color = exp), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), linewidth = 0.5, se = FALSE) + ggsci::scale_color_npg() + ggsci::scale_fill_npg() +
    scale_y_continuous(limits = c(0.0, 1.5)) + xlab(expression(paste("Scaled density ", (italic(Paramecium))))) + ylab("Effect of Pa on Dn") + labs(tag = expression((italic(b)))) + theme_st(lunit = 2)

gp5 <- sim |> filter(tar == "Pa") |> ggplot(aes(x = Dn, y = SI)) + facet_wrap(. ~ method) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_smooth(aes(fill = exp, group = factor(run)), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), color = NA, alpha = 0.01) +
    geom_smooth(aes(color = exp), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), linewidth = 0.5, se = FALSE) + ggsci::scale_color_npg() + ggsci::scale_fill_npg() +
    scale_y_continuous(limits = c(-1.0, 0.5)) + xlab(expression(paste("Scaled density ", (italic(Didinium))))) + ylab("Effect of Dn on Pa") + theme_st(lunit = 2)

gp6 <- sim |> filter(tar == "Dn") |> ggplot(aes(x = Dn, y = SI)) + facet_wrap(. ~ method) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_smooth(aes(fill = exp, group = factor(run)), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), color = NA, alpha = 0.01) +
    geom_smooth(aes(color = exp), method = "lm", formula = y ~ poly(x, degree = 2, raw = TRUE), linewidth = 0.5, se = FALSE) + ggsci::scale_color_npg() + ggsci::scale_fill_npg() +
    scale_y_continuous(limits = c(0.0, 1.5)) + xlab(expression(paste("Scaled density ", (italic(Didinium))))) + ylab("Effect of Pa on Dn") + theme_st(lunit = 2)

(((gp3 / gp5) | (gp4 / gp6)) + plot_layout(guides = "collect")) |> ggsaver("fig06", width = 17.3, height = 11, ext = "pdf")

