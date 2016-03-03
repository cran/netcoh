\name{predict.rncReg}
\alias{predict.rncReg}
\title{make predictions from a "rncReg" object.
}
\description{The generic prediction function that makes predictions
  according an "rncReg" object, which is a regression model fitted with
  network cohesion.
}
\usage{
\method{predict}{rncReg}(object,full.X=NULL,full.A,...)
}
\arguments{
  \item{object}{An object returned by the function \code{rncreg}.}

  \item{full.X}{
    A matrix with first n rows being the original training covariates
    and the last m rows being the new test covariates. If not provided,
    covariates will not be used in prediction.
}
  \item{full.A}{
An adjacency matrix of the complete network with both training and test
    samples. The first n vertices are representing the training samples
    (in the original order) and the last m vertices are reprenting the
    test samples.
  }
  \item{...}{further arguments passed to or from other methods.}
}

\details{
  The function first predicts individual effects on test sample by
  minimizing network cohesion penalty and then incorporates the covariate effects.
  For full details, please check the reference paper. The predicted test
  sample individual effects will be returned, as well as the
  corresponding linear terms. For linear
  regression model, the predicted response \code{y} is also given, which
  is exactly the same as the linear term. For logistic regression, the
  predicted probability is also given.
}
\value{ A list with following slots:
  \item{\code{terms}}{the linear term in the model.}
  \item{\code{alpha}}{the predicted individual effects.}
  \item{\code{y}}{the predicted responses in linear model.}
  \item{\code{p}}{the predicted probabilities in logistic regression.}
  \item{\code{model}}{the model used in prediction. This is the same as
  in \code{object}.}
}
\references{Tianxi Li, Elizaveta Levina and Ji Zhu. (2016)
  \emph{Regression with network cohesion},
  \url{http://arxiv.org/pdf/1602.01192v1.pdf}
}
\author{Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@umich.edu}
}
\seealso{\code{rncreg}
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


## If one wants to regularize the Laplacian by 
## using gamma > 0 in rncreg, we suggest use 
## centered data.
#mean.x <- colMeans(X1)
#mean.y <- mean(Y1)
#Y1 <- Y1-mean.y
#X1 <- t(t(X1)-mean.x)
#Y <- Y-mean.y
#X <- t(t(X)-mean.x)



m <- rncreg(A=A1,X=X1,Y=Y1,model="linear",lambda=10,gamma=0,cv=5)

p <- predict(m,full.A=A,full.X=X)
  }
\keyword{models}
\keyword{networks}
\keyword{regression}


