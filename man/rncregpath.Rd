\name{rncregpath}
\alias{rncregpath}
\title{Fits regression models with network cohesion effects for a
sequence of cohesion levels.}
\description{Fits models by \code{rncreg} for a sequence of lambda values.
}
\usage{
rncregpath(A,lambdaseq,Y=NULL,X=NULL,dt=NULL,gamma=0.05,
     model=c("linear","logistic","cox"), max.iter=50,tol=1e-4,
     init=NULL, cv=NULL,cv.seed=999,low.dim=NULL,verbose=FALSE)
}
\arguments{
  \item{A}{An nxn symmetric adjacency matrix for the network.}

  \item{lambdaseq}{
A vector of lambda values.
}
  \item{Y}{
An nx1 matrix of response, if the model to fit is linear or logistic. It
will not be used if one fits Cox's model.
}
  \item{X}{An nxp covariate matrix with each row being the covariates
    for one individual. If one want to fits a model without using any
    covariate, it can be empty.
  }
  \item{dt}{
Only used to fit Cox's model. An nx2 data.frame such that the first
column is the observed time while the second column is the event
indicator which is 1 for truely observed events and 0 for censored events.
}
  \item{gamma}{
    The amount of diagonal regularization added to graph Laplacian.
  }
  \item{model}{
Can only be one of "linear", "logistic" or "cox".
}
  \item{max.iter}{
The maximum number of newton steps to iterate. Only used for logistic model or Cox's model.
}
  \item{tol}{
The tolerance level for convergence. Only used for logistic model or Cox's model.
}
\item{init}{
  The initial point for newton algoritm. It should be an (n+p)x1 matrix
  that stacks alpha and beta. Only used for logistic model or Cox's
  model. If not specified, zeros will be used.
}
  \item{cv}{
Number of folds for cross-validation. If unspecified, then no
cross-validation will be done.
}
\item{cv.seed}{Random number generator seed for cross-validation.
}
\item{low.dim}{Only used for linear model. If the probelm is a low
dimensional problem such that n>>p, then using low.dim=TRUE is
potentially faster.}

\item{verbose}{If TRUE, the log likelihood in each newton step will be
  printed. Only used for logistic model and Cox's model.}

}

\details{
  The function repeatly calls \code{rncreg} for the values of lambda.
}
\value{ A list will be returned with the following slots:
	\item{models}{A list of objects from \code{rncReg} class, each of which corresponds to a specific lambda value in lambdaseq. }
	\item{cv.seq}{The sequence of cross-validation loss
	  corresponding to lambdaseq.}
	\item{cv.sd}{The sequence of corssivalidation standard
	  deviation}
	\item{cv.min.index}{The index of model that corresponds to the
	  one with miminum cv loss.}
	\item{cv.1sd.index}{The index of model that corresponds to the
  one selected by 1 sigma rule.}
}
\references{Tianxi Li, Elizaveta Levina and Ji Zhu. (2016)
  \emph{Regression with network cohesion},
  \url{http://arxiv.org/pdf/1602.01192v1.pdf}
}
\author{Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@umich.edu}
}
\seealso{\code{rncReg},\code{predict.rncReg}
}
\examples{
set.seed(100)

A <- matrix(0,200,200)
A[1:100,1:100] <- 1
A[101:200,101:200] <- 1
diag(A) <- 0

alpha <- c(rep(1,100),rep(-1,100)) + rnorm(200)*0.5
A <- A[c(1:50,101:150,51:100,151:200),c(1:50,101:150,51:100,151:200)]
alpha <- alpha[c(1:50,101:150,51:100,151:200)]

beta <- rnorm(2)

X <- matrix(rnorm(400),ncol=2)

Y <- X%*%beta + alpha + rnorm(200)

A1 <- A[1:100,1:100]
X1 <- X[1:100,]
Y1 <- matrix(Y[1:100],ncol=1)


## If one wants to regularize the Laplacian
## by using gamma > 0 in rncreg,
## we suggest use centered data.
#mean.x <- colMeans(X1)
#mean.y <- mean(Y1)
#Y1 <- Y1-mean.y
#X1 <- t(t(X1)-mean.x)
#Y <- Y-mean.y
#X <- t(t(X)-mean.x)



m.list <- rncregpath(A=A1,X=X1,Y=Y1,model="linear",
                     lambdaseq=exp(seq(0,log(200),length.out=20)),gamma=0,cv=5)
m.list$cv.seq
  }
\keyword{models}
\keyword{networks}
\keyword{regression}


