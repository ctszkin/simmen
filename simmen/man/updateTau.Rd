\name{updateTau}
\alias{updateTau}
\title{updateTau}
\usage{
updateTau(tau, update_rate, lower_bound = 0.3, upper_bound = 0.4,
  optim_rate = 0.35, min_rate = 0.001)
}
\arguments{
  \item{tau}{previous tau}

  \item{update_rate}{update rate of last estimation}

  \item{lower_bound}{update if the update_rate <
  lower_bound}

  \item{upper_bound}{update if the update_rate >
  upper_bound}

  \item{optim_rate}{targeted updated rate}

  \item{min_rate}{minimum value of tau}
}
\value{
value
}
\description{
Update the value of tau (scale of the sampling
distribution)
}
\author{
TszKin Julian Chan \email{ctszkin@gmail.com}
}

