\name{simmen}
\alias{simmen}
\title{Social Interaction Model with Multiple Endogenous Networks}
\usage{
simmen(m, data, last_out, allow_correlation = TRUE, update_sign = TRUE,
  initial_seed = 100, seed_for_chain = round(runif(1) * 10000),
  start_value_delta_diff_e = 0, tau_e = 1)
}
\arguments{
  \item{m}{Number of iteration}

  \item{data}{data}

  \item{last_out}{previous estimation}

  \item{allow_correlation}{Logical. Allow correlation
  between network}

  \item{update_sign}{Logical. Whether to update sign of e}

  \item{initial_seed}{seed for the first round}

  \item{seed_for_chain}{generate seed for each chain}

  \item{start_value_delta_diff_e}{start value of the
  coefficient of abs(ei-ej)}

  \item{tau_e}{sd of sampling distribution}
}
\value{
value
}
\description{
Social Interaction Model with Multiple Endogenous Networks
}
\author{
TszKin Julian Chan \email{ctszkin@gmail.com}
}

