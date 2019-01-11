#' sigmoid
#'
#'
#' Takes in a vector and applies the sigmoid function elementwise to that vector
#'
#' @name sigmoid
#' @param x : numeric vector
#' @return sigmoid(x)
#' @keywords internal

sigmoid <- function(x) {
    1/(1 + exp(-x))
}
