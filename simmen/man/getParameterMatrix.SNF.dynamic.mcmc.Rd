\name{getParameterMatrix.SNF.dynamic.mcmc}
\alias{getParameterMatrix.SNF.dynamic.mcmc}
\title{getParameterMatrix.SNF.dynamic.mcmc}
\usage{
\method{getParameterMatrix}{SNF.dynamic.mcmc}(x, tail, ...)
}
\arguments{
  \item{x}{SNF.dynamic.mcmc object}

  \item{tail}{iteration to be used. Negative value:
  Removing the first \code{tail} iterations. Positive
  value: keep the last \code{tail} iterations. If -1<
  code{tail}< 1, it represent the percentage of iterations.
  '}

  \item{...}{not used}
}
\value{
A matrix
}
\description{
Get a matrix of parameter
}
\author{
TszKin Julian Chan \email{ctszkin@gmail.com}
}

