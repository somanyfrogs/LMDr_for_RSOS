#' @title Figure_2.R
#' @description Making figure 2 for the Paper entitled:
#'      "Local-manifold-distance-based regression: An estimation method for quantifying dynamic biological interactions with empirical time series"
#'      Initially written on 20230822.
#'      Last update: 20240531.

## Load R functions
source("R/functions.R")

## Figure 2: A schematic diagram of LMDr analysis
## Ex. one-prey-two-predator system with predator-dependent functional response (Beddington 1975; DeAngelis et al 1975)
## Case 1: multifaceted state dependency caused by different interaction strengths

## Set the model parameters
r1 <- 3.75; r2 <- 3.61; r3 <- 3.65
k <- 0.1; e <- 0.25

params <- list(r1 = r1, r2 = r2, r3 = r3, k = k, e = e,
               a12 = 0.4, a13 = 0.1, g12 = 0.2, g13 = 0.2)

## Generate time-series data and make embedding matrix
ts <- gen_map(x_ini = rep(0.2, 3), fun = predatorFR, params = params, t_end = 250) |> filter(time > 50)
X <- ts |> gen_emat(cols = str_c("x", 1:3), lags = rep(0, 3))

## Perform LMDr and S-map variants for the dataset
system.time(ned_res <- smap_mulvar(X = X, range = seq(0, 50, 0.1), drop_dist = FALSE))
system.time(gdd_res <- smap_geodes(X = X, seed = 123, iter = 500, k_range = 3:25, sigmas = c(0.05, 0.05), drop_dist = FALSE))
system.time(mvd_res <- get_mvd(X, cols = 1:3, tar = 1, E = 3, lmax = 4) |> smap_mulvar(X = X, dmat = _, range = seq(0, 50, 0.1), drop_dist = FALSE))
system.time(lmd_res <- find_best_theta(X = X, gSeed = 123, method = 1, diag = 1))

## Make figure 2a
tmp1 <- X[26:50, ] |> as_tibbler(str_c("x", 1:3)) |> mutate(time = 1:length(x1)) |> pivot_longer(!time, names_to = "sp", values_to = "val") |>
    ggplot(aes(x = sp, y = time, fill = val)) + geom_tile(color = "grey50") + scale_fill_viridis(option = "cividis", guide = "none") +
    scale_x_discrete(labels = c(expression(italic(x[1])), expression(italic(x[2])), expression(italic(x[3]))), expand = rep(0, 2)) + scale_y_continuous(expand = rep(0, 2)) +
    theme_st() + theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(0, 0, 0, 0))

tmp2 <- lmd_res$theta[[1]] |> as_tibbler(1:3) |> mutate(row = 1:3) |> pivot_longer(!row, names_to = "col", values_to = "theta") |>
    mutate(across(1:2, as.numeric)) |> ggplot(aes(x = col, y = rev(row), fill = theta)) + geom_tile(color = "grey50") +
    scale_fill_viridis(option = "cividis", guide = "none") + scale_x_continuous(expand = rep(0, 2)) + scale_y_continuous(expand = rep(0, 2)) +
    theme_st() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(0, 0, 0, 0))

tmp3 <- (X[26:50, ] %*% lmd_res$theta[[1]]) |> as_tibbler(str_c("x", 1:3)) |> mutate(time = 1:length(x1)) |> pivot_longer(!time, names_to = "sp", values_to = "val") |>
    ggplot(aes(x = sp, y = time, fill = val)) + geom_tile(color = "grey50") + scale_fill_viridis(option = "cividis", guide = "none") +
    scale_x_discrete(labels = c("", "", ""), expand = rep(0, 2)) + scale_y_continuous(expand = rep(0, 2)) +
    theme_st() + theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(0, 0, 0, 0))

gp1 <- tibble(x = 0, y = 0) |> ggplot(aes(x = x, y = y)) + labs(tag = expression((italic(a)))) + theme_st() + theme_wt2 + coord_cartesian(clip = "off") +
    scale_x_continuous(expand = rep(0, 2), limits = c(-1, 1)) + scale_y_continuous(expand = rep(0, 2), limits = c(-1, 1)) +
    annotate(geom = "text", x = c(-0.72, -0.12), y = c(1.0, 1.0), label = c(expression(bold(X)), expression(bold(Theta))), size = 3, family = "Times New Roman") +
    annotate(geom = "text", x = 0.28, y = 0.6, label = "=", family = "Helvetica") +
    annotation_custom(ggplotGrob(tmp1), xmin = -1.0, xmax = -0.5, ymin = -1.1, ymax = 0.9) +
    annotation_custom(ggplotGrob(tmp2), xmin = -0.4, xmax =  0.1, ymin =  0.27, ymax = 0.9) +
    annotation_custom(ggplotGrob(tmp3), xmin =  0.4, xmax =  0.9, ymin = -1.1, ymax = 0.9)

## Make figure 2b
## Calculate Geodesic distance and LMD based on the above results
## Obtain Distance matrix for each method
ned <- ned_res$dist[[1]]
gdd <- gdd_res$dist[[1]]
mvd <- mvd_res$dist[[1]]
lmd <- get_lmd(X, lmd_res$theta[[1]])

## Make figure 2b
nedw <- exp(-ned_res$theta * ned / rowMeans(ned, na.rm = TRUE)) |> as_tibbler(1:ncol(ned)) |> mutate(row = 1:nrow(ned), method = "Euclidean") |>
    pivot_longer(!c(row, method), values_to = "weight", names_to = "col") |> mutate(col = as.numeric(col))
gddw <- exp(-gdd_res$theta * gdd / rowMeans(gdd, na.rm = TRUE)) |> as_tibbler(1:ncol(gdd)) |> mutate(row = 1:nrow(gdd), method = "Geodesic") |>
    pivot_longer(!c(row, method), values_to = "weight", names_to = "col") |> mutate(col = as.numeric(col))
mvdw <- exp(-mvd_res$theta * mvd / rowMeans(mvd, na.rm = TRUE)) |> as_tibbler(1:ncol(mvd)) |> mutate(row = 1:nrow(mvd), method = "Multiview") |>
    pivot_longer(!c(row, method), values_to = "weight", names_to = "col") |> mutate(col = as.numeric(col))
lmdw <- exp(-lmd / rowMeans(ned, na.rm = TRUE)) |> as_tibbler(1:ncol(lmd)) |> mutate(row = 1:nrow(lmd), method = "LMD") |>
    pivot_longer(!c(row, method), values_to = "weight", names_to = "col") |> mutate(col = as.numeric(col))

gp2 <- bind_rows(nedw, gddw, mvdw, lmdw) |> mutate(method = factor(method, levels = c("LMD", "Multiview", "Geodesic", "Euclidean")), weight = if_else(weight <= 1e-10, log(1e-10), log(weight))) |>
    filter(between(col, 26, 50) & between(row, 26, 50)) |> ggplot(aes(x = col, y = rev(row), fill = weight)) + facet_wrap(. ~ method, nrow = 1) + geom_tile() +
    scale_fill_viridis(option = "inferno", breaks = c(0, log(1e-2), log(1e-4), log(1e-8)), labels = c(1, expression(1^-2), expression(1^-4), expression(1^-8))) +
    scale_x_continuous(expand = rep(0, 2)) + scale_y_continuous(expand = rep(0, 2)) + labs(tag = expression((italic(b))), fill = expression(log(italic(w[ij])))) + theme_st(lunit = 2) + theme_wt2

## Make figure 2c
ts |> filter(between(x2, 0.58, 0.62) & between(x3, 0.70, 0.75))
x_tar <- ts |> slice(51) |> select(!time) |> as.numeric()
J_true <- foreach(i = 1:nrow(ts), .combine = rbind) %do% {ts |> slice(i) |> select(!time) |> as.numeric() |> predatorFRJ(params = params) |> (\(.x) .x[1, ])()}

gp3 <- ts |> mutate(J11 = J_true[, 1], J12 = J_true[, 2], J13 = J_true[, 3]) |> filter(between(x2, x_tar[2] - 0.2, x_tar[3] + 0.2) & between(x3, x_tar[3] - 0.2, x_tar[3] + 0.2)) |>
    ggplot(aes(x = x2, y = x3)) + geom_segment(aes(xend = x2 + 0.3 * J12, yend = x3 + 0.3 * J13, color = J11), arrow = arrow(length = unit(0.5, "mm"), type = "closed")) +
    geom_point(shape = 21, size = 0.75, fill = "white") + scale_color_viridis(option = "inferno") + annotate(geom = "point", x = x_tar[2], y = x_tar[3], shape = 21, fill = ggsci::pal_npg()(3)[1]) +
    xlab(expression(italic(x[2]))) + ylab(expression(italic(x[3]))) + labs(tag = expression((italic(c))), color = expression(italic(J[11]))) +
    theme_st() + theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.text = element_text(hjust = 1))

J_tar <- predatorFRJ(x = x_tar, params = params)[1, ]

vf_sim <- foreach(i = 1:nrow(ts), .combine = c) %do% {
    J_i <- slice(ts, i) |> select(!time) |> as.numeric() |> predatorFRJ(params = params)
    sqrt(sum((J_i[1, ] - J_tar)^2))
}

tmp <- ts |> mutate(vf_sim = vf_sim,
                    Euclidean = exp(-ned_res$theta * ned[51, ] / mean(ned[51, ], na.rm = TRUE)),
                    Geodesic = exp(-gdd_res$theta * gdd[51, ] / mean(gdd[51, ], na.rm = TRUE)),
                    Multiview = exp(-mvd_res$theta * mvd[51, ] / mean(mvd[51, ], na.rm = TRUE)),
                    LMD = exp(-0.25 * lmd[51, ] / mean(ned[51, ], na.rm = TRUE)))

gp4 <- tmp |> pivot_longer(c(Euclidean, Geodesic, Multiview, LMD), names_to = "method", values_to = "weights") |>
    mutate(method = factor(method, levels = c("LMD", "Multiview", "Geodesic", "Euclidean"))) |> ggplot(aes(x = vf_sim, y = log(weights))) + facet_wrap(. ~ method, nrow = 1) +
    geom_smooth(method = "lm", color = NA, fill = ggsci::pal_npg()(3)[1], alpha = 0.75) + geom_point(shape = 16, size = 0.25) +
    xlab(expression(paste("Dissimilarity ||", italic(J[1][i](t)) - italic(J[1][i])(tar), "||"[2]))) + ylab(expression(log(italic(w[ij])))) + theme_st(lunit = 2)

gp <- (gp1 + gp2 + plot_layout(width = c(1, 4))) / (gp3 + gp4 + plot_layout(width = c(1, 3))) + plot_layout(height = c(4, 5))
gp |> ggsaver("fig02", width = 17.3, height = 8.8, ext = "pdf")

