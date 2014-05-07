library(Jmisc)
library(SESHK2011)
library(Formula)
library(mvtnorm)
library(compiler)
library(parallel)
library(foreach)
library(MASS)
library(rbenchmark)
library(microbenchmark)
library(MCMCpack)
library(maxLik)
library(RcppEigen)
library(plyr)
library(truncnorm)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(inline)
library(Matrix)

sourceCpp("simmen/src/function.cpp")

library(simmen)

data = genModelData(genSpecStandardMath())
data_friends = genModelData(genSpecStandardMath("friends"))
data_studymates = genModelData(genSpecStandardMath("studymates"))

a = network_formation(data, method="maxLik")
b = network_formation(data_friends, method="maxLik")
c = network_formation(data_studymates)

d = network_formation(data,method="mcmc",50)
e = network_formation(data_friends,method="mcmc",50)
f = network_formation(data_studymates,method="mcmc",50)

e = network_formation(data_friends,method="mcmc",50,last=e)
d = network_formation(data_friends,method="mcmc",50,last=d)


summary(a)
summary(b)
summary(c)
summary(d)
summary(e)
summary(f)



a = SNF(50, data)
b = SNF(50, data, last = a)

summary(a)
summary(b)


multi_endo = multi_endogenous_mcmc(100, data)

multi_endo = multi_endogenous_mcmc(100, data_friends)

