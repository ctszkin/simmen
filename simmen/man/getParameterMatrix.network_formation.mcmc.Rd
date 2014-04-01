\name{getParameterMatrix.network_formation.mcmc}
\alias{getParameterMatrix.network_formation.mcmc}
\title{getParameterMatrix.network_formation.mcmc}
\usage{
\method{getParameterMatrix}{network_formation.mcmc}(x, tail, ...)
}
\arguments{
  \item{x}{network_formation}

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

