#' Check OpenMP
#' 
#' Print a message from each thread to check if OpenMP works.     
#' Code written by Drew Schmidt.
#' @export

c_hello <- function() {
  invisible(.Call("c_hello", PACKAGE = "pmdr"))
}


#' Parallel Manhattan distance
#' 
#' Parallel Manhattan distance using C and OpenMP.
#' This is a basic distance computation for integer values only.
#' It has been written in order to compare different looping strategies
#' and their impact on computation time.
#' This function is experimental and has a limited scope. It is not meant 
#' to be used for a serious purpose.
#' 
#' @param x An integer matrix. Missing values are not allowed.
#' @param loop Possible values are "standard", "colwise", and "diagwise".
#' With "standard", the code is run using two loops (for the two rows being 
#' compared). With "colwise" and "diagwise", the code is run using a single 
#' loop on an index then used to compute the row indices, in a 
#' "col-wise" and "diag-wise" way, respectively.
#' 
#' @return
#' An integer vector of distances for all the pairwise combinations.
#' Note that \code{loop = "diagwise"} will return values in a different order.
#' 
#' @export

distance <- function(x, loop = c("standard", "colwise", "diagwise")) {

  loop <- match.arg(loop)

  if (!is.integer(x)) {
    stop("The input matrix must contain integer values.")
  }
  if (anyNA(x)) {
    stop("Missing values are not allowed in the input matrix.")
  }

  switch(loop,
         standard = .Call("c_dist_two_loops", x, PACKAGE = "pmdr"),
         colwise  = .Call("c_dist_one_loop_colwise", x, PACKAGE = "pmdr"),
         diagwise = .Call("c_dist_one_loop_diagwise", x, PACKAGE = "pmdr"))

}

