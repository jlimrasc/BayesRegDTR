#' Generate a toy dataset in the right format for testing BayesLinRegDTR.model.fit
#'
#' Generates a toy dataset simulating observed data of treatments over time with
#' final outcomes and intermediate covariates. Follows the method outlined in
#' \href{https://github.com/jlimrasc/BayesRegDTR/blob/main/man/Toy_Datagen.pdf}{Toy-Datagen on Github}
#'
#' @param n             Number of samples/individuals to generate
#' @param num_stages    Total number of stages per individual
#' @param p_list        Vector of dimension for each stage
#' @param num_treats    Vector of number of treatment options at each stage
#'
#' @returns Observed data organised as a list of \eqn{\{y, X_1, X_2..., X_{num\_stages}, A\}} where y is a
#' vector of the final outcomes, \eqn{X_1, X_2..., X_{num\_stages}} is a list of matrices
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
#' p_list_uvt  <- rep(1, num_stages)
#' p_list_mvt  <- c(1, 3, 3)
#' num_treats  <- rep(3, num_stages)
#'
#' # -----------------------------
#' # Main
#' # -----------------------------
#' Data_uvt    <- generate_dataset(n, num_stages, p_list_uvt, num_treats)
#' Data_mvt    <- generate_dataset(n, num_stages, p_list_mvt, num_treats)
generate_dataset <- function(n, num_stages, p_list, num_treats) {
    if (any(p_list > 1))
        return (generate_dataset_mvt(n = n, num_stages = num_stages, p_list = p_list, num_treats = num_treats))
    else
        return (generate_dataset_uvt(n = n, num_stages = num_stages, num_treats = num_treats))
}
