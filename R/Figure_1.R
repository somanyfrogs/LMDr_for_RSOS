#' @title Figure_1.R
#' @description Making figure 1 for the Paper entitled:
#'      "Local-manifold-distance-based regression: An estimation method for quantifying dynamic biological interactions with empirical time series"
#'      Initially written on 20230822.
#'      Last update: 20240227.

## Load R functions
source("R/functions.R")

## Figure 1: A multifaceted state dependency of ecological networks
## Ex. one-prey-two-predator system with predator-dependent functional response (Beddington 1975; DeAngelis et al 1975)
## Case 1: multifaceted state dependency caused by different interaction strengths
## Case 2: multifaceted state dependency caused by different functional responses

## Set the model parameters
r1 <- 3.75; r2 <- 3.61; r3 <- 3.65
k <- 0.1; e <- 0.25

param1 <- list(r1 = r1, r2 = r2, r3 = r3, k = k, e = e,
               a12 = 0.4, a13 = 0.1, g12 = 0.2, g13 = 0.2)
param2 <- list(r1 = r1, r2 = r2, r3 = r3, k = k, e = e,
               a12 = 0.4, a13 = 0.4, g12 = 0.2, g13 = 1.0)

## Make figure 1a (Upper: Network graph)
nodes <- tibble(id = str_c("x", 1:3), type = c(TRUE, FALSE, FALSE))
edges <- tibble(From = c(2, 3), tp = c(1, 1))

gp1 <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE) |> ggraph(layout = "igraph", algorithm = "bipartite") +
    geom_edge_link(arrow = arrow(length = unit(2, "mm"), type = "closed"), end_cap = circle(4, "mm")) +
    geom_node_point(size = 4, color = "antiquewhite") + geom_node_text(aes(label = id), size = 2) +
    scale_x_continuous(expand = rep(0.1, 2)) + scale_y_continuous(expand = rep(0.1, 2)) + labs(tag = expression((italic(a)))) +
    theme_wt1() + theme(panel.background = element_rect(color = "white"), panel.border = element_rect(color = "white"))

gp2 <- bind_rows(with(param1, tibble(case = "case1", xj = seq(0.0, 2.0, length.out = 101), J12 = a12 * xj / (1 + k * a12 * 0.5 + g12 * xj), J13 = a13 * xj / (1 + k * a13 * 0.5 + g13 * xj))),
                 with(param2, tibble(case = "case2", xj = seq(0.0, 2.0, length.out = 101), J12 = a12 * xj / (1 + k * a12 * 0.5 + g12 * xj), J13 = a13 * xj / (1 + k * a13 * 0.5 + g13 * xj)))) |>
    pivot_longer(!c(case, xj)) |> ggplot(aes(x = xj, y = value, color = name)) + facet_wrap(. ~ case, nrow = 1) + geom_line() +
    ggsci::scale_color_npg(label = c(expression(italic(x[2])), expression(italic(x[3]))), guide = guide_legend(title = NULL)) + xlab("Predator density") + ylab("Consumption rate") +
    theme_st(lunit = 2) + theme(axis.text = element_blank(), axis.ticks = element_blank())

## Make figure 1b (Vector field in the x2-x3 plane)
Jtar1 <- predatorFRJ(x = rep(0.5, 3), params = param1)
Jtar2 <- predatorFRJ(x = rep(0.5, 3), params = param2)
tbl_tar <- tibble(case = c("case1", "case2"), x1 = 0.5, x2 = 0.5, x3 = 0.5,
                  J11 = c(Jtar1[1, 1], Jtar2[1, 1]), J12 = c(Jtar1[1, 2], Jtar2[1, 2]), J13 = c(Jtar1[1, 3], Jtar2[1, 3]))

ful <- expand_grid(x1 = 0.5, x2 = seq(0.0, 1.0, length.out = 101), x3 = seq(0.0, 1.0, length.out = 101))
sub <- expand_grid(x1 = 0.5, x2 = seq(0.1, 1.0, length.out = 6), x3 = seq(0.06, 0.94, length.out =  6))

grid_ful <- foreach(tbl = iter(ful, by = "row"), .combine = bind_rows) %do% {
    J1 <- predatorFRJ(x = as.numeric(tbl), params = param1)
    J2 <- predatorFRJ(x = as.numeric(tbl), params = param2)
    op1 <- tbl |> mutate(case = "case1", J11 = J1[1, 1], J12 = J1[1, 2], J13 = J1[1, 3], dissimilarity = sqrt((Jtar1[1, 1] - J11)^2 + (Jtar1[1, 2] - J12)^2 + (Jtar1[1, 3] - J13)^2))
    op2 <- tbl |> mutate(case = "case2", J11 = J2[1, 1], J12 = J2[1, 2], J13 = J2[1, 3], dissimilarity = sqrt((Jtar2[1, 1] - J11)^2 + (Jtar2[1, 2] - J12)^2 + (Jtar2[1, 3] - J13)^2))
    bind_rows(op1, op2)
}

grid_sub <- foreach(tbl = iter(sub, by = "row"), .combine = bind_rows) %do% {
    J1 <- predatorFRJ(x = as.numeric(tbl), params = param1)
    J2 <- predatorFRJ(x = as.numeric(tbl), params = param2)
    op1 <- tbl |> mutate(case = "case1", J11 = J1[1, 1], J12 = J1[1, 2], J13 = J1[1, 3], dissimilarity = sqrt((Jtar1[1, 1] - J11)^2 + (Jtar1[1, 2] - J12)^2 + (Jtar1[1, 3] - J13)^2))
    op2 <- tbl |> mutate(case = "case2", J11 = J2[1, 1], J12 = J2[1, 2], J13 = J2[1, 3], dissimilarity = sqrt((Jtar2[1, 1] - J11)^2 + (Jtar2[1, 2] - J12)^2 + (Jtar2[1, 3] - J13)^2))
    bind_rows(op1, op2)
}

gp3 <- grid_ful |> ggplot(aes(x = x2, y = x3)) + facet_wrap(. ~ case, nrow = 1) +
    geom_raster(aes(fill = dissimilarity), interpolate = TRUE) + geom_contour(aes(z = dissimilarity), linewidth = 0.1) +
    geom_segment(data = grid_sub, aes(xend = x2 + 0.3 * J12, yend = x3 + 0.3 * J13), linewidth = 0.1, arrow = arrow(length = unit(1.0, "mm"), type = "closed")) +
    geom_segment(data = tbl_tar, aes(xend = x2 + 0.3 * J12, yend = x3 + 0.3 * J13), linewidth = 0.25, arrow = arrow(length = unit(1.0, "mm"), type = "closed")) +
    geom_point(data = tbl_tar, shape = 21, fill = ggsci::pal_npg()(3)[1]) + scale_fill_viridis(direction = -1, breaks = c(0.0, 0.1, 0.2, 0.3), limits = c(0, 0.35)) +
    scale_x_continuous(expand = rep(2e-3, 2)) + scale_y_continuous(expand = rep(2e-3, 2)) + labs(tag = expression((italic(b)))) +
    xlab(expression(italic(x[2]))) + ylab(expression(italic(x[3]))) + theme_st(lunit = 2) + theme(axis.text = element_blank(), axis.ticks = element_blank())

gp <- (gp1 / gp2 | gp3) + plot_layout(width = c(1, 3))
gp |> ggsaver("fig01", width = 17.3, height = 6, ext = "pdf")
