\name{BayesRegDTR-package}
\alias{BayesRegDTR-package}
\alias{BayesRegDTR}
\docType{package}
\title{
\packageTitle{BayesRegDTR}
}
\description{
\packageDescription{BayesRegDTR}
}
\author{
\packageAuthor{BayesRegDTR}

Maintainer: \packageMaintainer{BayesRegDTR}
}
\references{
References: Yu, W., & Bondell, H. D. (2023), “Bayesian Likelihood-Based Regression
for Estimation of Optimal Dynamic Treatment Regimes”, Journal of the Royal
Statistical Society Series B: Statistical Methodology, 85(3), 551-574. Page 3
\doi{doi:10.1093/jrsssb/qkad016}
}
\keyword{ package }
\seealso{
\code{\link[=generate_dataset]{generate_dataset()}} for generating a test dataset

\code{\link[=BayesLinRegDTR.model.fit]{BayesLinRegDTR.model.fit()}} for obtaining an estimated posterior
distribution of the optimal treatment option at a user-specified prediction stage

}
\examples{

#Data   <- generate_dataset(n, num_stages, p_list, num_treats)
#result <- BayesLinRegDTR.model.fit(Dat.train, Dat.pred, n.train, n.pred, num_stages, num_treats,
#p_list, t, R, numCores, tau, B, nu0, V0, alph, gam)
}
