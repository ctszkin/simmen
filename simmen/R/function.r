#' @import Jmisc SESHK2011 MASS foreach mvtnorm parallel maxLik plyr truncnorm Matrix MCMCpack Formula
NULL

#' @useDynLib simmen
NULL

#' Flexible version of tail
#' @name extractTail
#' @detail if  |x| < 1, x will be considered as a percentage, otherwise, it will be return tail(x,tail)
#' @aliases extractTail
#' @title extractTail
#' @param x x
#' @param tail default = 0 
#' @return value
#' @seealso link \code{\link{tail}}
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
extractTail = function(x, tail=0){
  if (tail==0){
    return(x)
  }
  if( abs(tail) < 1){
    tail= round( tail*NROW(x) )
  }
  tail(x,tail)
}

#' computeSummaryTable
#' @name computeSummaryTable
#' @aliases computeSummaryTable
#' @title computeSummaryTable
#' @detail computeSummaryTable
#' @param x x
#' @param tail tail
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
computeSummaryTable = function(x, tail){
  if (!is.matrix(x)){
    x= getParameterMatrix(x,tail)
  }
  # x = extractTail(x, tail)

  mean = round(colMeans(x),6)
  sd = round(apply(x,2,sd),6)
  out = cbind(mean,sd)
  name = rownames(out)
  generateSignificance(out,name)
}

#' getParameterMatrix
#' @name getParameterMatrix
#' @aliases getParameterMatrix
#' @title getParameterMatrix
#' @detail getParameterMatrix
#' @param x x
#' @param ... ... 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
getParameterMatrix = function(x,...) UseMethod("getParameterMatrix", x)

#' difference of log normal density between x and y 
#' @name dnorm.diff
#' @aliases dnorm.diff
#' @title dnorm.diff
#' @detail dnorm.diff
#' @param x x
#' @param y y
#' @param log logical. log density?
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
dnorm.diff = function(x, y, log=TRUE){
  dnorm(x,log=log) - dnorm(y,log=log)
}

addStar = function(x){
  x$star 
}

#' Check whether folder name from 1 to n exists
#' @name checkFolder
#' @aliases checkFolder
#' @title checkFolder
#' @detail checkFolder
#' @param n n 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
checkFolder = function(n=100){
  folder_exist = rep(F,n)
  for (i in 1:n){
    folder_exist[i] = file.exists(as.character(i)) # %+% "/last.rData")
  }
  if (!any(folder_exist)){
    return(0)
  }
  which(folder_exist)
}

#' What is the largest value of folder name. (Folder name must be integer)
#' @name checkMaxFolder
#' @aliases checkMaxFolder
#' @title checkMaxFolder
#' @detail checkMaxFolder
#' @param n max integer for folder name
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
checkMaxFolder = function(n=20){
  folder_exist = rep(F,n)
  for (i in 1:n){
    folder_exist[i] = file.exists(as.character(i))
  }
  if (!any(folder_exist)){
    return(0)
  }
  max(which(folder_exist))
}


#' commandArgsParser
#' @name commandArgsParser
#' @aliases commandArgsParser
#' @title commandArgsParser
#' @detail This function will check whether the commandArgs contain "-ID". If it does, the argument will be return, if not, it will return the max folder number +1. that's mean if there is not argument, it will use a new folder used in beginning of runme.r to determine where to store the mcmc reuslt
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
commandArgsParser = function(){
  command_args = commandArgs(F)
  command_args = gsub("-ID","",command_args[grep("-ID", command_args)])

  if ( length(command_args) ==0 ) {
    command_args = as.character(checkMaxFolder() +1)
  }
  command_args
}

#' gen an unique ID
#' @name genUniqueID
#' @aliases genUniqueID
#' @title genUniqueID
#' @detail it depends on the time and a 8 digits random integer
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
genUniqueID = function(){
  format(Sys.time(), "%y%m%d%H%M%S") %+% "_" %+% as.character(round(runif(1)*10000000))
}

#' betterTable
#' @name betterTable
#' @aliases betterTable
#' @title betterTable
#' @detail betterTable
#' @param x x
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
betterTable = function(x){
  name = rownames(x)
  out = signif(x$estimates,2) %+% "(" %+% signif(x$sd,2) %+% ")" %+% x$significance
  names(out) = name
  out
}


#' mergeTable
#' @name mergeTable
#' @aliases mergeTable
#' @title mergeTable
#' @detail mergeTable
#' @param x x
#' @param y y
#' @param seperator seperator
#' @param title title
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
mergeTable = function(x,y, seperator=c("friends_","studymates_"), title=c("Friends","Studymates")){

  name_x = names(x)
  name_y = names(y)

  same_name = intersect(name_x,name_y)

  simpleMergeTable1(x[same_name],y[same_name])

  same_table = simpleMergeTable1(x[same_name],y[same_name])

  x2 = x[! name_x %in% same_name]
  y2 = y[! name_y %in% same_name]

  names(x2) = gsub(seperator[1],"", names(x2))
  names(y2) = gsub(seperator[2],"",names(y2))

  same_name2 = intersect(names(x2),names(y2))


  same_table2 = simpleMergeTable1(x2[same_name2],y2[same_name2])

  rownames(same_table2) = seperator[1]%+% same_name2

  out = rbind(same_table,same_table2)
  colnames(out) = title
  out
}

#' SplitTable
#' @name SplitTable
#' @aliases SplitTable
#' @title SplitTable
#' @detail SplitTable
#' @param x x
#' @param seperator seperator
#' @param title title
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
SplitTable = function(x, seperator = c("friends_","studymates_"), title=c("Friends", "Studymates")){

  name = names(x)

  index1 = grep(seperator[1], name)
  index2 = grep(seperator[2], name)

  x1 = x[-index2]
  x2 = x[-index1]

  mergeTable(x1,x2)
}

#' simpleMergeTable1
#' @name simpleMergeTable1
#' @aliases simpleMergeTable1
#' @title simpleMergeTable1
#' @detail simpleMergeTable1
#' @param x x
#' @param y y 
#' @param colname colname
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
simpleMergeTable1 = function(x,y, colname=c("x","y")){
  name_x = names(x)
  name_y = names(y)

  id1 = data.frame(name=name_x,id=seq(length(x)))
  id2 = data.frame(name=name_y,id=seq(length(y)))

  xx = data.frame(name=name_x, x )
  yy = data.frame(name=name_y, y)

  out = join(xx,yy,by="name",type="full")
  rownames(out) = out$name
  out = out[c("x","y")]
  colnames(out) = colname
  out
}


#' Generate random variable from truncated normal
#' @name my.rtnorm
#' @aliases my.rtnorm
#' @title my.rtnorm
#' @detail truncated normal
#' @param mean mean
#' @param sd sd
#' @param a a 
#' @param b b 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
my.rtnorm=function(mean=0,sd=1 , a, b){
  n = length(mean)
  axb= pnorm((a-mean)/sd) 
  bxb= pnorm((b-mean)/sd)

  mean + qnorm(runif(n) * (bxb-axb) + axb)*sd
}

# test

# q1=replicate(10000, my.rtnorm(1:5, sd=2, a=2,b=4))
# q2=replicate(10000, rtruncnorm(1,a= 2,b=4,mean=as.numeric(1:5),sd=2))

# rowMeans(q1)- rowMeans(q2)
# apply(q1,1,sd)- apply(q2,1,sd)


# splitBy = cmpfun(function(x,by){
#   start = c(0,head(cumsum(by),-1)) +1
#   end = cumsum(by)
#   out = list()
#   for (i in 1:length(by)){
#     if (is.matrix(x))
#       out[[i]] = x[start[i]:end[i],]
#     else 
#       out[[i]] = x[start[i]:end[i]]
#   }
#   out
# })


#' transform a vector of error into absolute difference form
#' @name transform_e_by_each
#' @aliases transform_e_by_each
#' @title transform_e_by_each
#' @detail abs(ei-ej)
#' @param data data
#' @param e e
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
transform_e_by_each = function(data,e){
  e_i = e[data$group_index[,1]]
  e_j = e[data$group_index[,2]]

  abs(e_i - e_j )    #- mean(abs(e_i - e_j ) )
  # (e_i - e_j )^2
}

#' Apply transform_e_by_each to each network
#' @name transform_e
#' @aliases transform_e
#' @title transform_e
#' @detail transform_e
#' @param data data
#' @param e e
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
transform_e =function(data,e){
  out = vector("list",length(data))
  for (i in 1:length(data)){
    out[[i]] = transform_e_by_each(data[[i]],e[[i]])
  }
  unlist(out)
}

#' Apply transform_e2_by_each to each network
#' @name transform_e2
#' @aliases transform_e2
#' @title transform_e2
#' @detail transform_e2
#' @param data data
#' @param e e
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
transform_e2 =function(data,e){
  out = vector("list",length(data))
  for (i in 1:length(data)){
    out[[i]] = apply(e[[i]], 2, transform_e_by_each, data=data[[i]])
  }
  (out)
}




#' genPoly
#' @name genPoly
#' @aliases genPoly
#' @title genPoly
#' @detail genPoly
#' @param e_i e_i
#' @param e_j e_j
#' @param delta delta
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
genPoly = function(e_i,e_j,delta){
  e_ij_diff = (e_i - e_j )^2 * delta
  cbind(e_ij_diff + e_i , e_ij_diff + e_j )
}


#' Return e_i e_j and |e_i - e_j|
#' @name genFulleByNetwork
#' @aliases genFulleByNetwork
#' @title genFulleByNetwork
#' @detail genFulleByNetwork
#' @param data data
#' @param e e 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
genFulleByNetwork = function(data,e){
  e_i=NULL
  e_j=NULL
  for (i in 1:length(data)){
    e_i = c(e_i, e[[i]][data[[i]]$group_index[,1]] )
    e_j = c(e_j, e[[i]][data[[i]]$group_index[,2]] )
  }
  list(e_i=e_i, e_j=e_j, e_diff=abs(e_i-e_j))
}

#' genFullGroupIndex
#' @name genFullGroupIndex
#' @aliases genFullGroupIndex
#' @title genFullGroupIndex
#' @detail genFullGroupIndex
#' @param data data
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
genFullGroupIndex = function(data){
  n = sapply(data, "[[", i ="n")

  adjust_factor =  c(0,cumsum(head(n,-1)))

  out=NULL
  for(i in 1:length(data)){
    tmp = data[[i]]$group_index + adjust_factor[i]
    out = rbind(out,tmp)
  }
  out
}

#' genFullPositionIndex
#' @name genFullPositionIndex
#' @aliases genFullPositionIndex
#' @title genFullPositionIndex
#' @detail genFullPositionIndex
#' @param  data data
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
genFullPositionIndex = function(data){
  n = sapply(data,"[[","n")
  full_group_index = genFullGroupIndex(data)
  out = vector("list",sum(n))
  for (i in 1:sum(n)){
    out[[i]] = which(full_group_index==i)
  }
  out
}

#' genFullPositionMatrix
#' @name genFullPositionMatrix
#' @aliases genFullPositionMatrix
#' @title genFullPositionMatrix
#' @detail genFullPositionMatrix
#' @param data data
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
genFullPositionMatrix = function(data){
  n = sapply(data,"[[","n")
  full_group_index = genFullGroupIndex(data)
  out = NULL
  for (i in 1:sum(n)){
    out = rbind(out, cbind(i,which(full_group_index==i)))
  }
  sparseMatrix(i=out[,1], j =out[,2])
}


# e = lapply(sapply(data,"[[","n"), rnorm)
# full_e1 = genFulleByNetwork(data,e)
# full_e2 = genFulle(data,e)


# all.equal(full_e1[[1]],full_e2[[1]],check.attributes = FALSE)
# all.equal(full_e1[[2]],full_e2[[2]],check.attributes = FALSE)
# all.equal(full_e1[[3]],full_e2[[3]],check.attributes = FALSE)

#' Update the value of tau (scale of the sampling distribution)
#' @name updateTau
#' @aliases updateTau
#' @title updateTau
#' @detail updateTau
#' @param tau previous tau
#' @param update_rate update rate of last estimation
#' @param lower_bound update if the update_rate < lower_bound
#' @param upper_bound update if the update_rate > upper_bound
#' @param optim_rate targeted updated rate
#' @param min_rate minimum value of tau
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
updateTau = function(tau, update_rate, lower_bound=0.3, upper_bound=0.4, optim_rate=0.35, min_rate =0.001){
  for ( j in 1:length(tau)){
    if ( is.list(tau[[j]])  ){
      tau[[j]] = updateTau(tau[[j]], update_rate[[j]], lower_bound=lower_bound, upper_bound=upper_bound, optim_rate=optim_rate, min_rate =min_rate)
    } else if (any(update_rate[[j]] > upper_bound | any(update_rate[[j]]< lower_bound) )){
      cat("update tau-", names(tau)[[j]] , "\n")
      tau[[j]] = tau[[j]] * update_rate[[j]] / optim_rate
      tau[[j]] = ifelse(tau[[j]]==0, min_rate, tau[[j]])
    }
  }
  tau
}

# mean and var of i given j, a is the demean value of j

#' Compute the conditional normal distribution 
#' @name find_normal_conditional_dist
#' @aliases find_normal_conditional_dist
#' @title find_normal_conditional_dist
#' @detail Conditional distribution of x1 given x2, where x1 and x2 are joint normal
#' @param a mean 
#' @param i index or location of x1
#' @param j index or location of x2
#' @param Sigma Covariance matrix
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
find_normal_conditional_dist = function(a, i, j, Sigma){
  if (NROW(Sigma)==1){
    return(list(mean=0,var=matrix(1,1,1)))
  }
  Sigma_12 = Sigma[i,j] 
  S = Sigma_12 %*% solve(Sigma[j,j])
  Sigma_new = Sigma[i,i] - S %*% t(Sigma_12)
  if (!missing(a)){
    return(  list(mean=as.vector(tcrossprod(S,a[,j]) ), var = Sigma_new))
  } else {
    return(  list(var = Sigma_new))
  }
}

#' Check whether there is folder named as \code{command_args}
#' @name checkNewStart
#' @aliases checkNewStart
#' @title checkNewStart
#' @detail checkNewStart
#' @param command_args folder name
#' @param create_dir logical. If the folder doesn't exist, create a new directory with name \code{command_args}
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
checkNewStart = function(command_args, create_dir=TRUE){
  if (file.exists(command_args)) {
    return(FALSE)
  } else {
    if (!file.exists(command_args)){
      dir.create(command_args)
    }
    return(TRUE)
  }
}
