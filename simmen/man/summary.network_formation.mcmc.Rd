\name{summary.network_formation.mcmc}
\alias{summary.network_formation.mcmc}
\title{summary.network_formation.mcmc}
\usage{
\method{summary}{network_formation.mcmc}(object, ...)
}
\arguments{
  \item{object}{network_formation object}

  \item{...}{tail: iteration to be used. Negative value:
  Removing the first \code{tail} iterations. Positive
  value: keep the last \code{tail} iterations. If -1<
  code{tail}< 1, it represent the percentage of
  iterations.}
}
\value{
A summary table
}
\description{
Description Create a summary table
}
\author{
TszKin Julian Chan \email{ctszkin@gmail.com}
}

