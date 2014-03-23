# #' genPairwiseIndex 
# #' @name genPairwiseIndex
# #' @aliases genPairwiseIndex
# #' @title genPairwiseIndex
# #' @param n
# #' @return data.frame
# #' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
# #' @keywords internal 
# #' @export
# genPairwiseIndex = cmpfun(function(n){
#   # i<- foreach( i=2:n,.combine=c ) %do% seq(from=i,to=n)
#   i = unlist(lapply(2:n, seq, to=n))
#   j = rep(1:(n-1),times=(n-1):1)
#   cbind(i,j)
# })


# #' getPairwiseFriendshipData 
# #' @name getPairwiseFriendshipData
# #' @aliases getPairwiseFriendshipData
# #' @title getPairwiseFriendshipData
# #' @param network_data network_data
# #' @param network_formation_formula network_formation_formula
# #' @return A data.frame consist of the estimates.
# #' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
# #' @keywords internal 
# #' @export
# getPairwiseFriendshipData<-cmpfun(function(network_data,network_formation_formula){
#   data<-network_data$data
#   n <- nrow(data)
#   require_variable <-unique(sub("friends_","",all.vars(network_formation_formula)))

#   require_variable <- setdiff(require_variable,".")

#   data<-data[require_variable]
  
#   pairwise_index = genPairwiseIndex(n)

#   i = pairwise_index[,1]
#   j = pairwise_index[,2]
  
#   name_of_self_data <- names(data)
#   name_of_friends_data <- paste("friends_",name_of_self_data,sep="")
  
#   self_data <- data.frame(data[i,])
  
#   friends_data <- data.frame(data[j,])
#   names(self_data)<-name_of_self_data
#   names(friends_data)<-name_of_friends_data
  
#   self_data_matrix <- model.matrix(network_formation_formula,cbind(self_data,friends_data))

#   names(friends_data)<-name_of_self_data
#   names(self_data)<-name_of_friends_data
  
#   friends_data_matrix <- model.matrix(network_formation_formula,cbind(self_data,friends_data))

#   D<-network_data$network_matrix_list
#   response_self = sapply(D, function(x) x[lower.tri(x)] )
#   response_friends = sapply(D, function(x) t(x)[lower.tri(t(x))] )

#   response_self=!!response_self
#   response_friends=!!response_friends

#   list(response_self=response_self, response_friends=response_friends,self_data_matrix=self_data_matrix,friends_data_matrix=friends_data_matrix)
# })



# #' genPairwiseHobbyData 
# #' @name genPairwiseHobbyData
# #' @aliases genPairwiseHobbyData
# #' @title genPairwiseHobbyData
# #' @param H H
# #' @return a vector
# #' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
# #' @keywords internal 
# #' @export
# genPairwiseHobbyData<-cmpfun(function(H){
#   tHH<-tcrossprod(H)
#   # stopifnot(nrow(tHH)==ncol(tHH))
#   n<-nrow(tHH)

#   index<-genPairwiseIndex(n)
#   index2= (index[,1]-1) * n + index[,2]

#   tHH[index2]
# })


# #' generateWeighting 
# #' @name generateWeighting
# #' @aliases generateWeighting
# #' @title generateWeighting
# #' @param D a network matrix
# #' @return weighting matrix
# #' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
# #' @keywords internal 
# #' @export
# generateWeighting<-cmpfun(function(D){
#   W<-D/rowSums(D)
#   W[is.nan(W)]<-0
#   W
# })

# prepareData <- function (.raw_data, .spec, .school ){

#     # .hobby <- match.arg(.hobby, several.ok = TRUE)
#     # .school <- match.arg(.school, several.ok = TRUE)
#     if ( missing(.school) ) {
#       out <- foreach(i = .spec$school_name, .combine = c) %do% {
#         prepareData(.raw_data=.raw_data, .school = i, .spec=.spec)
#       }
#       return(out)
#     }

#     data_wide <- getDataWide(.raw_data, .school)

#     network_matrix_list = 
#     lapply(.spec$network_info_list, function(x){
#       network_matrix1 = getNetwork(.raw_data, .school, x$definition[1])
#       network_matrix2 = getNetwork(.raw_data, .school, x$definition[2])
#       out = x$process_network(network_matrix1, network_matrix2, data_wide)
#       out
#     })

#     drop_case_id = .spec$findDropCaseID(data_wide, network_matrix_list)


#     data_wide <- getDataWide(.raw_data, .school, drop_case_id)

#     network_matrix_list = 
#     lapply(.spec$network_info_list, function(x){
#       network_matrix1 = getNetwork(.raw_data, .school, x$definition[1], .drop_by_case_id=drop_case_id)
#       network_matrix2 = getNetwork(.raw_data, .school, x$definition[2], .drop_by_case_id=drop_case_id)
#       # out = process_network_function_example(network_matrix1, network_matrix2, data_wide)
#       out = x$process_network(network_matrix1, network_matrix2, data_wide)
#       out
#     })



#     H <- foreach(i = .spec$hobby, .combine = c) %do% {
#         out <- list(getHobby(.raw_data, .school, i, drop_case_id))
#         names(out) <- i
#         out
#     }

#     network_name = sapply(.spec$network_info_list, "[[", "name")

#     group_index = list(id=which( .school == .spec$school_name ), all=length(.spec$school_name))

#     # degree = sapply(network_matrix_list, rowSums)

#     # colnames(degree) = paste0(    network_name = sapply(.spec$network_info_list, "[[", "name") , "_degree")

#     # data_wide= cbind(data_wide, degree)


#     data_wide = cbind(data_wide, .spec$genNetworkStatistics(network_matrix_list) )


#     out <- list(network_matrix_list = network_matrix_list, data = data_wide, H = H, network_name=network_name,group_index=group_index)
#     out <- list(out)
#     names(out) <- .school
#     return(out)
# }







# extractData <- cmpfun(function(spec, data){
    
#   # source("model_spec.r")
#   formula = spec$formula
#   network_formation_formula = Formula(spec$network_formation_formula)

  
#   f1 = formula(network_formation_formula,rhs=1)

#   ## if there is fixed effect, take out the intercept
#   if ( !is.null(spec$network_formation_fixed_effect) && spec$network_formation_fixed_effect ){
#     f1 = update(f1, ~.+-1)
#   }

#   other_network_variables<-NULL
#   use_network_variable<-FALSE
#   if (length(network_formation_formula)[2]==2){
#     use_network_variable<-TRUE
#     f2<-formula(network_formation_formula,rhs=2)
#     other_network_variables<-attr(terms(f2),"term.labels")
#   }
#   H_name <-  attr(terms(f2),"term.labels")
  
#   out<-getPairwiseFriendshipData(data,f1)
#   out$formula = formula
#   out$network_formation_formula = network_formation_formula


#   if (length(network_formation_formula)[2]==2){
#     H_pair<-sapply(data$H[H_name],genPairwiseHobbyData)
#     out$self_data_matrix<-cbind(out$self_data_matrix,H_pair)
#     out$friends_data_matrix<-cbind(out$friends_data_matrix,H_pair)
#   }
  
#   ## generate dummy matrix
#   if ( !is.null(spec$network_formation_fixed_effect) && spec$network_formation_fixed_effect ){
#     nn = NROW(out$self_data_matrix)
#     dummy_matrix = matrix(0, nrow=nn, ncol = data$group_index$all)
#     dummy_matrix[,data$group_index$id] = 1

#     colnames(dummy_matrix) = spec$school_names

#     out$self_data_matrix = cbind(dummy_matrix,out$self_data_matrix)
#     out$friends_data_matrix = cbind(dummy_matrix,out$friends_data_matrix)
#   }

#   # TODO: add seat assignment here
#   # out$response1 <- as.logical(out$response)
#   # out$response <- list(response1=out$response1)
  
#   # if(!single_network){
#   #   out$response2 <- as.logical(out2$response)
#   #   out$response$response2 = out$response2
#   # }

# # single network
#   out$D_list = data$network_matrix_list
#   out$W_list = lapply(out$D_list, generateWeighting)

#   y_x_wx = getXandWX(formula,data)
#   out$y = y_x_wx$y
#   out$x = y_x_wx$x
#   out$wx = y_x_wx$wx
#   out$wy = y_x_wx$wy
#   out$x_wx = y_x_wx$X

#   out$n = length(out$y)
#   out$k_x_wx = ncol(out$x_wx)
#   out$k_gamma = ncol(out$self_data_matrix)

#   out$group_index = genPairwiseIndex(length(out$y))
#   out$network_name = data$network_name 

#   return(out)
# })



# tic <- function(){
#   assign(".time", Sys.time(), envir= .GlobalEnv)
# }

# toc <- function(){
#   # as.numeric(Sys.time() - .time, units = "secs")
#   print(Sys.time() - .time)
# }

## should put this in Jmisc
exportAllFunction <- function(cl, envir=.GlobalEnv){
  all_stuff = ls(envir=envir )
  is_fun = sapply(lapply(all_stuff,get), is.function)
  clusterExport(cl,all_stuff[is_fun])
}

dnorm.diff = function(x, y, log=TRUE){
  dnorm(x,log=log) - dnorm(y,log=log)
}

# genDataMatrix = cmpfun(function(data, any_wx = TRUE){    
#   Y = unlist( sapply(data, function(z) z$y ) )
#   WY = do.call( rbind, lapply(data, "[[", "wy") )

#   if (any_wx){
#     X = do.call(rbind, lapply(data, function(z) z$x_wx ))
#   } else {
#     X = do.call(rbind, lapply(data, function(z) z$x ))
#   }

#   n_vector = sapply(data, function(z) z$n )
#   X = demean(X)
#   Y = demean(Y)
#   WY = demean(WY)

#   W_list = lapply(data, "[[", "W_list")

#   list(
#     Y = Y,
#     X = X,
#     WY = WY,
#     W_list = W_list,
#     n_vector = n_vector
#   )
# })

# genNetworkData = cmpfun(function(data){
#   lapply(data, function(x){
#     location_index1 = lapply(1:x$n, function(z){ which(x$group_index[,1] == z) })
#     location_index2 = lapply(1:x$n, function(z){ which(x$group_index[,2] == z) })
#     location_index_all = lapply(1:x$n, function(z){ c(location_index1[[z]], location_index2[[z]])} )
#     list(location_index1=location_index1,location_index2=location_index2,location_index_all=location_index_all)
#   })
# })



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

# genModelData = cmpfun(function(.spec, save=FALSE){
#   if (class(.spec)!="SESHK_Spec")
#     stop("spec must be an SESHK_Spec object")

#    raw_data<-readSeshkNetwork(.version = "../../" %+% .spec$data_version %+% "/")

#    data = prepareData(.raw_data=raw_data, .spec=.spec )

#    out = lapply(data, extractData, spec=.spec)
   
#   if (save){
#     cat("Saving data to model_data.rData\n")
#     save(out,file="model_data.rData")
#   }
#   return(out)
# })

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


# updateTau = cmpfun(function(tau,update_rate, less_update=FALSE){
#   for ( j in 1:length(tau)){
#     adjust_rate = update_rate[[j]] / 0.3
#     if ( less_update){
#       to_update = (tau[[j]] < 0.25)  | (tau[[j]] > 0.4)
#       tau[[j]][to_update] =  tau[[j]][to_update] * adjust_rate[to_update]
#     } else {
#       tau[[j]] = tau[[j]] * adjust_rate
#       tau[[j]] = ifelse(tau[[j]]==0, 0.001, tau[[j]])
#     }
#   }
#   tau
# })

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
