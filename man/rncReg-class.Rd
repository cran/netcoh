\name{rncReg-class}
\Rdversion{0.1}
\docType{class}
\alias{rncReg-class}

\title{Class \code{"rncReg"}}
\description{
A generic class for regression objects with network cohesion.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rncReg", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{The individual effects of the regression.}
    \item{\code{beta}:}{The fixed effects or covariate coefficients of
      the regression.}
    \item{\code{A}:}{The network adjacency matrix for which cohession is
      assumed.}
    \item{\code{lambda}:}{Parameter for cohesion penalty. }
    \item{\code{X}:}{Covariate matrix.}
    \item{\code{Y}:}{Response matrix. }
    \item{\code{dt}:}{The response data frame with the first column
      being the observed time and the second column being the event indicator. }
    \item{\code{gamma}:}{Regularization parameter for graph Laplacian.}
    \item{\code{cv}:}{Number of folds in cross-validation. }
    \item{\code{cv.loss}:}{Cross-validated prediciton loss. It is MSE
      for linear regression, binomial deviance for logistic regression
      and test partial loglikelihood for Cox's model (see reference
      paper).}
    \item{\code{cv.sd}:}{Standard deviation of cross-validation
      loss. It can be used for cross-validation by 1 sigma rule. It is
      more robust to noises.}
    \item{\code{model}:}{The specific regression model used.}
  }
}
\section{Methods}{
No methods defined with class "rncReg" in the signature.
}
\references{
Tianxi Li, Elizaveta Levina and Ji Zhu. (2016)
  \emph{Regression with network cohesion},
  \url{http://arxiv.org/pdf/1602.01192v1.pdf}\cr
  Verweij, Pierre JM, and Hans C. Van Houwelingen. (1993)
  \emph{Cross-validation in survival analysis}, Statistics in medicine 12, no. 24: 2305-2314.
}
\author{
Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@umich.edu}
}
\seealso{\code{rncreg}}
\keyword{classes}
