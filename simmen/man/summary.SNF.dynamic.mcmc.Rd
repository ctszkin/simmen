\name{summary.SNF.dynamic.mcmc}
\alias{summary.SNF.dynamic.mcmc}
\title{summary.SNF.dynamic.mcmc}
\usage{
\method{summary}{SNF.dynamic.mcmc}(object, ...)
}
\arguments{
  \item{object}{SNF.dynamic.mcmc object}

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
Create a summary table
}
\author{
TszKin Julian Chan \email{ctszkin@gmail.com}
}

