% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesRegDTR-package.R
\docType{package}
\name{BayesRegDTR-package}
\alias{BayesRegDTR}
\alias{BayesRegDTR-package}
\title{BayesRegDTR: Optimal Dynamic Treatment Regime Inference using Bayesian Likelihood-Based Regression Method}
\description{
Methods to estimate optimal dynamic treatment regimes using Bayesian likelihood-based regression approach as described in Yu, W., & Bondell, H. D. (2023). Uses backward induction and dynamic programming theory for computing expected values. Offers options for future parallel computing.
}
\references{
References: Yu, W., & Bondell, H. D. (2023), “Bayesian Likelihood-Based Regression
for Estimation of Optimal Dynamic Treatment Regimes”, Journal of the Royal
Statistical Society Series B: Statistical Methodology, 85(3), 551-574.
\doi{doi:10.1093/jrsssb/qkad016}
}
\seealso{
\code{\link[=generate_dataset]{generate_dataset()}} for generating a toy dataset to test the model fitting on

\code{\link[=BayesLinRegDTR.model.fit]{BayesLinRegDTR.model.fit()}} for obtaining an estimated posterior
distribution of the optimal treatment option at a user-specified prediction stage

Useful links:
\itemize{
\item \url{https://github.com/jlimrasc/BayesRegDTR}
\item Report bugs at \url{https://github.com/jlimrasc/BayesRegDTR/issues}
}
}
\author{
\strong{Maintainer}: Jeremy Lim \email{jeremylim23@gmail.com}

Authors:
\itemize{
  \item Weichang Yu \email{weichang.yu@unimelb.edu.au} (\href{https://orcid.org/0000-0002-0399-3779}{ORCID})
}

}
