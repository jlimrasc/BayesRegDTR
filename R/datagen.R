#' Generate Dataset
#'
#' Generate a dataset to test functions on
#'
#' @param n             Number of samples/individuals to generate
#' @param num_stages    Total number of stages per individual
#' @param p_list        Vector of dimension for each stage
#' @param num_treats    Vector of number of treatment options at each stage
#'
#' @returns Observed data organised as a list of \eqn{\{y, X_1, X_2..., X_{num\_states}, A\}} where y is a
#' vector of the final outcomes, \eqn{X_1, X_2..., X_{num\_states}} is a list of matrices
#' of the intermediate covariates and A is an \eqn{n \times num\_stages}{n x num_stages} matrix of the
#' assigned treatments
#' @export
#'
#' @examples
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' n           <- 5000
#' num_stages  <- 3
#' p_list      <- rep(2, num_stages)
#' num_treats  <- rep(3, num_stages)
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' Data        <- generate_dataset(n, num_stages, p_list, num_treats)
generate_dataset <- function(n, num_stages, p_list, num_treats) {
    if (any(p_list > 1))
        return (generate_dataset_uvt(n = n, num_stages = num_stages, num_treats = num_treats))
    else
        return (generate_dataset_mvt(n = n, num_stages = num_stages, p_list = p_list, num_treats = num_treats))
}
