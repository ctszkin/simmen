\name{summary.SNF}
\alias{summary.SNF}
\title{summary.SNF}
\usage{
\method{summary}{SNF}(object, ...)
}
\arguments{
  \item{object}{SNF object}

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

