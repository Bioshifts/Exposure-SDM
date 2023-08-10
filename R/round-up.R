#' Round number a la Unix
#'
#' @param x a `numeric` of length 1.
#' 
#' @param digits an `integer` of length 1 indicating the number of decimal 
#'   places
#'
#' @return A `numeric` of length 1.
#' @export
#' 
#' @examples
#' round_up(1.4, digits = 0)
#' round_up(1.5, digits = 0)
#' round_up(1.45, digits = 1)
#' round_up(-1.45, digits = 1)
#' round_up(-1, digits = 2)
#' round_up(1, digits = 2)

round_up <- function(x, digits) {
    
    posneg = sign(x)
    z = abs(x)*10^digits
    z = z + 0.5 + sqrt(.Machine$double.eps)
    z = trunc(z)
    z = z/10^digits
    z*posneg
    
}
