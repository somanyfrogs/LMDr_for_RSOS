#' @title zzz
#' @description Set the environment on load 'andsr'.
#'      Initially written on 20220501 by K.Kawatsu.
#'      Last update: 20220501.

.onLoad <- function(...) {
    if (!interactive()) return()

    intro_message <- paste("Development of 'andsr' is on going, and further improvements will be made in the future.\nPlease check my GitHub account (https://github.com/somanyfrogs/andsr)")
    packageStartupMessage(intro_message)
}
