\name{getParameterMatrix.SNF.static.mcmc}
\alias{getParameterMatrix.SNF.static.mcmc}
\title{getParameterMatrix.SNF.static.mcmc}
\usage{
\method{getParameterMatrix}{SNF.static.mcmc}(x, tail, ...)
}
\arguments{
  \item{x}{SNF.static.mcmc}

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
Description Get a matrix of parameter
}
\author{
TszKin Julian Chan \email{ctszkin@gmail.com}
}

