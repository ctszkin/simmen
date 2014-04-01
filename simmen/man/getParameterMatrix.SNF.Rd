\name{getParameterMatrix.SNF}
\alias{getParameterMatrix.SNF}
\title{getParameterMatrix.SNF}
\usage{
\method{getParameterMatrix}{SNF}(x, tail, ...)
}
\arguments{
  \item{x}{SNF object}

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

