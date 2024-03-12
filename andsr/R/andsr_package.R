#' @title Applications of analysis of nonlinear dynamical systems in R.
#' 
#' @description The \pkg{andsr} package is a new implementation of andsr methods,
#'      including the functions of nonlinear time-series analysis.
#'      those of lyapunov exponent/vector analysis (ands_lyap) and
#'      The \pkg{andsr} will be applicable to diverse purposes, including
#'      nonlinear time-series forecasting, system identification and stability analysis.
#'
#' @name andsr
#' @docType package
#' @useDynLib andsr, .registration=TRUE
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr arrange
#' @importFrom dplyr case_when
#' @importFrom dplyr group_by
#' @importFrom dplyr if_else
#' @importFrom dplyr lag
#' @importFrom dplyr lead
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr slice
#' @importFrom dplyr starts_with
#' @importFrom dplyr ungroup
#' @importFrom foreach foreach
#' @importFrom foreach %:%
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom glmnet cv.glmnet
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rcpp sourceCpp
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_sub
#' @importFrom tibble as_tibble
#' @importFrom tibble is_tibble
#' @importFrom tibble tibble
NULL
