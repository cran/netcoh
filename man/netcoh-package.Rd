\name{netcoh-package}
\alias{netcoh-package}
\alias{netcoh}
\docType{package}
\title{
Statistical models with network cohesion.
}
\description{
The package provides the methods to build a class of statistical models
on networks with network cohesion effects being taken into consideration.
}
\details{
\tabular{ll}{
Package: \tab netcoh\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2016-04-29\cr
License: \tab GPL (>= 2)\cr
}
The major function \code{rncreg} accepts \code{A} network adjacency matrix, \code{X} design matrix and
response \code{Y} for regression models and then produces the fitted modelfor the tuning
parameter \code{lambda}.

Another function \code{predict.rncReg} accepts the model fitted by
\code{rncreg} and the full design matrix as well as the full adjacency
matrix (with both training and test data).
}
\author{
Tianxi Li, Elizaveta Levina, Ji Zhu\cr

Maintainer: Tianxi Li <tianxili@umich.edu>
}
\references{
Tianxi Li, Elizaveta Levina and Ji Zhu. (2016)
  \emph{Regression with network cohesion},
  \url{http://arxiv.org/pdf/1602.01192v1.pdf}
}
\keyword{ package }
\keyword{models}
\keyword{regression}
\keyword{network}
