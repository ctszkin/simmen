#' @import Jmisc SESHK2011 MASS foreach mvtnorm parallel maxLik plyr truncnorm Matrix MCMCpack Formula
NULL

extractTail = function(x, tail=0){
  if (tail==0){
    return(x)
  }
  if( abs(tail) < 1){
    tail= round( tail*nrow(x) )
  }
  tail(x,tail)
}

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



getParameterMatrix = function(x,...) UseMethod("getParameterMatrix", x)

dnorm.diff = function(x, y, log=TRUE){
  dnorm(x,log=log) - dnorm(y,log=log)
}

addStar = function(x){
  x$star 
}

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

## This function will check whether the commandArgs contain "-ID". If it does, the argument will be return, if not, it will return the max folder number +1
## that's mean if there is not argument, it will use a new folder
## used in beginning of runme.r to determine where to store the mcmc reuslt
commandArgsParser = function(){
  command_args = commandArgs(F)
  command_args = gsub("-ID","",command_args[grep("-ID", command_args)])

  if ( length(command_args) ==0 ) {
    command_args = as.character(checkMaxFolder() +1)
  }
  command_args
}

genUniqueID = function(){
  format(Sys.time(), "%y%m%d%H%M%S") %+% "_" %+% as.character(round(runif(1)*10000000))
}

betterTable = function(x){
  name = rownames(x)
  out = signif(x$estimates,2) %+% "(" %+% signif(x$sd,2) %+% ")" %+% x$significance
  names(out) = name
  out
}

# betterTable(single_1_outcome[[2]])

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

SplitTable = function(x, seperator = c("friends_","studymates_"), title=c("Friends", "Studymates")){

  name = names(x)

  index1 = grep(seperator[1], name)
  index2 = grep(seperator[2], name)

  x1 = x[-index2]
  x2 = x[-index1]

  mergeTable(x1,x2)
}

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


transform_e_by_each = function(data,e){
  e_i = e[data$group_index[,1]]
  e_j = e[data$group_index[,2]]

  abs(e_i - e_j )    #- mean(abs(e_i - e_j ) )
  # (e_i - e_j )^2
}

transform_e =function(data,e){
  out = vector("list",length(data))
  for (i in 1:length(data)){
    out[[i]] = transform_e_by_each(data[[i]],e[[i]])
  }
  unlist(out)
}

transform_e2 =function(data,e){
  out = vector("list",length(data))
  for (i in 1:length(data)){
    out[[i]] = apply(e[[i]], 2, transform_e_by_each, data=data[[i]])
  }
  (out)
}





genPoly = function(e_i,e_j,delta){
  e_ij_diff = (e_i - e_j )^2 * delta
  cbind(e_ij_diff + e_i , e_ij_diff + e_j )
}


genFulleByNetwork = function(data,e){
  e_i=NULL
  e_j=NULL
  for (i in 1:length(data)){
    e_i = c(e_i, e[[i]][data[[i]]$group_index[,1]] )
    e_j = c(e_j, e[[i]][data[[i]]$group_index[,2]] )
  }
  list(e_i=e_i, e_j=e_j, e_diff=(e_i-e_j)^2)
}


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

genFulle = function(e, full_group_index){
  # full_group_index = genFullGroupIndex(data)
  if (is.list(e)){
    e= unlist(e)
  }

  e_i = e[full_group_index[,1]] 
  e_j = e[full_group_index[,2]]
  e_diff=(e_i-e_j)^2
  out = list(e_i=e_i, e_j=e_j, e_diff=(e_i-e_j)^2)
}

genFullPositionIndex = function(data){
  n = sapply(data,"[[","n")
  full_group_index = genFullGroupIndex(data)
  out = vector("list",sum(n))
  for (i in 1:sum(n)){
    out[[i]] = which(full_group_index==i)
  }
  out
}

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
