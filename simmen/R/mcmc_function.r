#' priorFunction
#' @name priorFunction
#' @aliases priorFunction
#' @title priorFunction
#' @details priorFunction
#' @param x x
#' @param type "normal","invgamma","uniform","trunc_normal"
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
priorFunction = function(x, type=c("normal","invgamma","uniform","trunc_normal")) { 
  type = match.arg(type)
  if (type=="normal")
    return(sum(dnorm(x,sd=2,log=TRUE)))
  if (type=="invgamma")
    return(dinvgamma(x, shape=2,scale=5))
  if (type=="uniform")
    return(ifelse(abs(x) <= 1, log(1/2), -Inf ) )
  if (type=="trunc_normal")
    return(sum(log(dtruncnorm(x,sd=2, a=0))))
}

#' samplingFunction
#' @name samplingFunction
#' @aliases samplingFunction
#' @title samplingFunction
#' @details samplingFunction
#' @param beta beta
#' @param tau tau 
#' @param type "normal","invgamma","uniform","trunc_normal"
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
samplingFunction = function(beta, tau, type=c("normal","invgamma","uniform","trunc_normal")) {
  type = match.arg(type)
  if (type=="normal")
    return(rnorm(length(beta))*tau + beta)
  if (type=="absnormal")
    return(abs(rnorm(length(beta),sd=tau) + beta))
  if (type=="uniform")
    return(runif(1,min=-tau,max=tau) +beta )
  if (type=="trunc_normal"){
    return(rtruncnorm(1, mean=beta, sd=tau, a=0))
  }
}

#' metropolis
#' @name metropolis
#' @aliases metropolis
#' @title metropolis
#' @details metropolis
#' @param beta_previous beta_previous
#' @param tau scale of sampling distribution
#' @param likelihoodFunction likelihoodFunction
#' @param prior_type Default is "normal" , see \link{priorFunction}
#' @param sampling_type Default is "normal", see \link{samplingFunction}
#' @param ... ...
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
metropolis = function(beta_previous, tau, likelihoodFunction, prior_type="normal", sampling_type="normal", ...){

  beta_new = samplingFunction(beta_previous, tau, sampling_type)

  lik_previous = likelihoodFunction(beta_previous, ...)
  lik_new = likelihoodFunction(beta_new, ...)

  # lik_previous = likelihoodFunction(beta_previous, data=data,sigma2=1)
  # lik_new = likelihoodFunction(beta_new, data=data,sigma2=1)
  # beta_new
  # exp(lik_new- lik_previous)


  prob = lik_new + priorFunction(beta_new, type=prior_type) - lik_previous - priorFunction(beta_previous, type=prior_type)

  prob = pmin(exp(prob),1)

  prob[!is.finite(prob)] = 0

  if ( prob > runif(1) ){
    return(list(beta_new=beta_new, beta=beta_new,update=1,prob=prob))
  } else{
    return(list(beta_new=beta_new, beta=beta_previous,update=0,prob=prob))
  }
}

#' metropolis2
#' @name metropolis2
#' @aliases metropolis2
#' @title metropolis2
#' @details metropolis2
#' @param beta_previous beta_previous
#' @param tau scale of sampling distribution
#' @param likelihoodFunction likelihoodFunction
#' @param prior_type Default is "normal" , see \link{priorFunction}
#' @param sampling_type Default is "normal", see \link{samplingFunction}
#' @param ... ... 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
metropolis2 = function(beta_previous, tau, likelihoodFunction, prior_type="normal", sampling_type="normal", ...){

  k = length(beta_previous)
  update = rep(0,k)

  beta_new = samplingFunction(beta_previous, tau, sampling_type)

  lik_previous = likelihoodFunction(beta_previous, ...)

  for (i in 1:k){
    beta_star = beta_previous
    beta_star[i] = beta_new[i]
    lik_star =  likelihoodFunction(beta_star, ...)
    if (lik_star > -1e+20) { 
      prob = lik_star + priorFunction(beta_star[i], type=prior_type) - lik_previous - priorFunction(beta_previous[i], type=prior_type)
      prob = exp(prob)
      if (prob > runif(1)){
        lik_previous = lik_star
        beta_previous = beta_star
        update[i] = 1
      }
    }
  }

  return(list(beta = beta_previous, update = update ))

}


#' findplotIndex
#' @name findplotIndex
#' @aliases findplotIndex
#' @title findplotIndex
#' @details findplotIndex
#' @param x x
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
findplotIndex = function(x){
  n=x/2
  index = 1:n* 2:(n+1) < x

  out = max(max(which(index)))+1

  return(c(out,out+1))
}


# transform.mcmc = function(x, remove=0){
#   i=j=NULL
#   all_name = 
#     foreach(i = x, j = names(x),.combine=c ) %do% {
#       if (is.matrix(i))
#         paste0(j,seq(ncol(i)))
#       else
#         j
#     }
#   x_matrix = do.call(cbind, x)
#   colnames(x_matrix) = all_name
#   if (remove>0)
#     return(tail(x_matrix,-remove))
#   else
#     return(x_matrix)  
# }

#' Inverse gamma Distribution
#' @name dinvgamma
#' @rdname invgamma
#' @aliases dinvgamma rinvgamma
#' @title invgamma
#' @details Inverse gamma Distribution
#' @param x x
#' @param n n
#' @param shape shape
#' @param scale scale
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
dinvgamma = function (x, shape, scale = 1) {
  if (shape <= 0 | scale <= 0) {
      stop("Shape or scale parameter negative in dinvgamma().\n")
  }
  alpha <- shape
  beta <- scale
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
      1) * log(x) - (beta/x)
  return(exp(log.density))
}

#' @rdname invgamma
#' @export
rinvgamma = function (n, shape, scale = 1) 
{
  if (scale==0 ) return(1)
  return(1/rgamma(n = n, shape = shape, rate = scale))
}

# dim2 = function(x, type=c("row","col")){
#   type = match.arg(type)
#   if (is.vector(x))
#     return(length(x))
#   else if (type=="row")
#     return(nrow(x))
#   else if (type=="col")
#     return(ncol(x))
# }

#' uniqueRow
#' @name uniqueRow
#' @aliases uniqueRow
#' @title uniqueRow
#' @details uniqueRow
#' @param x x
#' @param number_of_col number_of_col
#' @param percentage percentage
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
uniqueRow = function(x, number_of_col = "all", percentage=TRUE){
  if (is.numeric(number_of_col)){
    out = NROW(unique(x[,1:number_of_col]))
  } else {
    out = NROW(unique(x))
  }
  if (percentage){
    return(out/NROW(x))
  } else{
    return(out)
  }
}

#' load first/last object
#' @name loadLastObject
#' @rdname loadLastObject
#' @aliases loadLastObject loadFirstObject readAllLastObject
#' @title loadLastObject
#' @details load first/last/all last object. \code{readAllLastObject} will read all last object from every folder (with name from 1 to 20)
#' @param foldername foldername
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
loadLastObject = function(foldername){
  file_name = sort( setdiff(dir(foldername), "last.rData") , decreasing=TRUE)[1]
  obj_name = load(foldername %+% file_name)
  get(obj_name)
}

#' @rdname loadLastObject
#' @export
loadFirstObject = function(foldername){
  file_name = sort( setdiff(dir(foldername), "last.rData") , decreasing=FALSE)[1]
  obj_name = load(foldername %+% file_name)
  get(obj_name)
}

#' @rdname loadLastObject
#' @export
readAllLastObject = function(){
  i=NULL
  foreach (i = checkFolder()) %do% {
    loadLastObject(i %+% "/")
  }
}

#' readAllObject
#' @name readAllObject
#' @aliases readAllObject
#' @title readAllObject
#' @details readAllObject
#' @param path path
#' @param except vector of filenames for exclusion
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
readAllObject = function(path, except){
  env1 = new.env()
  file_name = dir(path)
  if (!missing(except)){
    file_name = setdiff(file_name,except)
  }

  sapply(path %+% file_name, load, envir = env1 )
  as.list(env1)
}

# readResult = function(path=".", start=6, data){
#   e1 = new.env()
#   sapply(path %+% "/" %+% dir(path), load, envir = e1)

#   name = createVariableName(data)
#   summary_list = 
#     lapply(start:length(as.list(e1)), function(i) {
#       out_obj = get("out"%+%i, envir = e1)
#       summary.mcmc(out_obj, name=name)
#     })

#   m = sum(
#     sapply(start:length(as.list(e1)), function(i) {
#       out_obj = get("out"%+%i, envir = e1)
#       nrow(out_obj$lambda_matrix)
#     })
#   )

#   mean_var = rowMeans( sapply(summary_list, function(i) i[,2] ) )
#   var_mean = apply(sapply(summary_list, function(i) i[,1 ]) , 1 , function(j)  mean((j-mean(j))^2) )

#   all_mean = rowMeans(sapply(summary_list, function(i) i[,1 ]))
#   all_var = mean_var + var_mean
#   all_sd = sqrt(all_var)
#   all_t = all_mean/all_sd

#   star <- ifelse(abs(all_t)>1.64,"*"," ") 
#   star <- star %+% ifelse(abs(all_t)>1.96,"*"," ") 
#   star <- star %+% ifelse(abs(all_t)>2.34,"*"," ")  

#   summary_table = data.frame(mean= all_mean, sd= sqrt(all_var),t=all_t, star=star)
#   row.names(summary_table) = name

#   out = list( summary_table=summary_table,   m=m)

#   # all_list = 
#   #   lapply(1:length(e1), function(i) {
#   #     out_obj = get("out"%+%i, envir = e1)
#   #     merge.mcmc(out_obj, name=name)
#   #   })
#   # q = do.call(rbind,all_list)
#   # colnames(q) = name

#   # all_mean2 = colMeans(q)
#   # all_var2 = apply(q,2,function(j)  mean((j-mean(j))^2) )


#   # sum(abs(all_mean-all_mean2) )
#   # sum(abs(all_var-all_var2))

#   return(out)
# }


# var2 = function(x){
#   if (is.matrix(x))
#     return(apply(x,2,var2))
#   if (is.vector(x))
#     return(mean(x^2)-mean(x)^2 )
# }

# colVar=function(x){
#   var2(x)
# }

# computeMeanVAR = function(x){
#   name = colnames(x)
#   out = cbind(mean=apply(x,2,mean), var=apply(x,2, function(z) mean((z-mean(z))^2) ))
#   rownames(out) = name
#   out
# }

# computeMeanSD = function(x){
#   name = colnames(x)
#   out = cbind(mean=apply(x,2,mean), sd=apply(x,2, sd ))
#   rownames(out) = name
#   out
# }


# summary.mcmc = function(x, name){
#   obj_name = c("phi_matrix","lambda_matrix","delta1_matrix","delta2_matrix")
#   out = 
#     do.call(rbind, 
#     do.call(rbind, 
#       lapply(obj_name, function(z) {
#         computeMeanVAR(x[[z]])  
#       })
#     )
#   rownames(out) = name
#   out
# }

#' plotmcmc_byid
#' @name plotmcmc_byid
#' @aliases plotmcmc_byid
#' @title plotmcmc_byid
#' @details plotmcmc_byid
#' @param x x
#' @param id id
#' @param mean mean
#' @param min min
#' @param max max
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
plotmcmc_byid = function(x, id, mean,min,max){
  m = length(x)
  name = colnames(x[[1]])[id]
  plot.new()
  par(mfrow=c(ceiling(m/2) ,2))
  for (i in 1:m){
    y = x[[i]][,id]
    plot(y,type='l', main = name %+% "="%+% round(mean(y),4), ylim=c(min,max))
    abline(h=mean(y),lwd=2,col="red")
    abline(h=mean,lwd=2,col="blue"  )
  }
}

#' plotmcmc_ma
#' @name plotmcmc_ma
#' @aliases plotmcmc_ma
#' @title plotmcmc_ma
#' @details plotmcmc_ma
#' @param x x
#' @param tail tail
#' @param ... ... 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export

plotmcmc_ma = function(x, tail,...){
  parameter_matrix_list = lapply(x, getParameterMatrix, tail=tail)
  j = 50
  par(mfrow=c(4,2))
  for (i in 1:8) {
    y = parameter_matrix_list[[i]][,j]
    plot( y , type="l")
    lines( cumsum(y) / 1:length(y), type = "l" , col="red")
  }

  par(mfrow=c(4,2))
  for (i in 1:8) {
    y = parameter_matrix_list[[i]][,j]
    plot( filter(y,rep(1/50000,50000) ,sides=1) , type="l")
    lines( cumsum(y) / 1:length(y), type = "l" , col="red")
  }

  par(mfrow=c(4,2))
  for (i in 1:8) {
    y = parameter_matrix_list[[i]][,j]
    plot( y , type="l")
    lines( filter(y,rep(1/50000,50000) ,sides=1) , type="l", col="red")
    # lines( cumsum(y) / 1:length(y), type = "l" , col="red")
  }



}

#' compareMCMC
#' @name compareMCMC
#' @aliases compareMCMC
#' @title compareMCMC
#' @details compareMCMC
#' @param x x
#' @param tail tail
#' @param plot Logical. Whether to have plot
#' @param pdf Logical
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
compareMCMC = function(x, tail, plot=FALSE, pdf=FALSE){
  if (missing(tail)){
    tail = NA
  }

  parameter_matrix_list = lapply(x, getParameterMatrix, tail=tail)

  summary_table_list = lapply(x, computeSummaryTable, tail=tail)


  parameter_matrix_all = do.call(rbind, parameter_matrix_list)

  summary_table_all = computeSummaryTable(parameter_matrix_all)


  all_mean = sapply(summary_table_list, function(z) z[,1] )
  all_var = sapply(summary_table_list, function(z) z[,2]^2 )

  between_mean_min = apply(all_mean,1,min)
  between_mean_max = apply(all_mean,1,max)
  between_mean_var = apply(all_mean,1,var)
  mean_within_var = rowMeans(all_var)

  ratio_var = between_mean_var/mean_within_var
  ratio_range_mean = (between_mean_max-between_mean_min) /2   / rowMeans(all_mean)

  converge1 = ratio_var < 0.2
  converge3 = abs(ratio_range_mean) < 0.2
  
  var_con  = ifelse(ratio_var < 0.1, "!!!", ifelse(ratio_var < 0.2, "!!", ifelse(ratio_var < 0.3, "!", "")))
  range_con  = ifelse(ratio_range_mean < 0.1, "###", ifelse(ratio_range_mean < 0.2, "##", ifelse(ratio_range_mean < 0.3, "#", "")))

  summary_table_simple = cbind(summary_table_all, 
    min=between_mean_min,
    max=between_mean_max,
    ratio_var=ratio_var,
    var_con=var_con
  )

  # summary_table_all = cbind(summary_table_all, 
  #   min=between_mean_min, 
  #   max=between_mean_max, 
  #   # between_mean_var=between_mean_var,
  #   # mean_within_var=mean_within_var,
  #   ratio_var=ratio_var,
  #   ratio_range_mean= ratio_range_mean, 
  #   converge1=converge1, 
  #   converge3=converge3 
  # )

  # summary_table_all = round(summary_table_all,4)
  non_converge = lapply(summary_table_list, function(z) z[!converge1, ] )
  non_converge_matrix_list = lapply(parameter_matrix_list, function(x) x[,!converge1] )

  non_converge_summary_table_simple = summary_table_simple[which(!converge1),]
  
  if (plot){
    if (file.exists("non_converge.pdf")){
      file.remove("non_converge.pdf")
    }
    if (pdf){
      pdf(file="non_converge.pdf",height=12,width=10)
    }
    for (j in which(!converge1)){
      mean = summary_table_simple[j,]$estimates
      min = min(parameter_matrix_all[,j])
      max = max(parameter_matrix_all[,j])
      plot.new()
      plotmcmc_byid(parameter_matrix_list,j,mean,min,max)
    }
    if (pdf){
      dev.off() 
    }
    # plotmcmc4(x,tail)
  }

  list(summary_table_simple = summary_table_simple, non_converge=non_converge, non_converge_summary_table_simple=non_converge_summary_table_simple)
}


#' plotmcmc.default
#' @name plotmcmc.default
#' @aliases plotmcmc.default
#' @title plotmcmc.default
#' @details plotmcmc.default
#' @param x x
#' @param tail tail
#' @param ... ...
#' @method plotmcmc default
#' @export plotmcmc
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
plotmcmc.default = function(x,tail=0,...){
  if (!is.matrix(x)){
    x = getParameterMatrix(x)
  }

  x = extractTail(x, tail)

  name=colnames(x)
  n = ncol(x)
  par(mfrow=findplotIndex(n))
  for (i in 1:ncol(x)){
    t = mean(x[,i])/sd(x[,i])
    t[!is.finite(t)] =0
    col = "black"
    if (abs(t) > 1.50) col = "yellow"
    if (abs(t) > 1.64) col = "blue"
    if (abs(t) > 1.96) col = "green"

    plot(x[,i], type="l", col=col)
    abline(h=0,col="purple", lty=2)
    abline(h=mean(x[,i]), col="red" , lwd=2 )
    if (!missing(name))
      title(name[i], sub=round(t,2))
  }
}

#' plotmcmc
#' @name plotmcmc
#' @rdname plotmcmc
#' @aliases plotmcmc plotmcmc2
#' @title plotmcmc
#' @details plotmcmc
#' @param x x
#' @param ... ... 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
plotmcmc = function(x,...) UseMethod("plotmcmc", x)




# mergeParameter = function(x,name,colnames){
#   out = do.call(cbind,x[name])
#   if (!missing(colnames))
#     colnames(out) = colnames
#   out
# }

# merge.mcmc = function(x, obj_name){
#   obj_name = c("phi_matrix","lambda_matrix","delta1_matrix","delta2_matrix","sigma2e_matrix")
#   out = 
#     do.call(cbind, 
#       lapply(obj_name, function(z) {
#         x[[z]]
#       })
#     )
#   colnames(out)
#   out
# }



#' @rdname plotmcmc 
#' @export
plotmcmc2 = function(x,...) UseMethod("plotmcmc2", x)



#' plotmcmc2.default
#' @name plotmcmc2.default
#' @aliases plotmcmc2.default
#' @title plotmcmc2.default
#' @details plotmcmc2.default
#' @param x x
#' @param name name
#' @param file file
#' @param tail tail
#' @param ... ... 
#' @method plotmcmc2 default
#' @export plotmcmc2 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
plotmcmc2.default = function(x, name, file="mcmc_plot.pdf", tail,...){
  if (!missing(tail)){
    x = lapply(x, function(z) tail(z,tail))
  }
  number_of_chain = length(x)
  number_of_variable = ncol(x[[1]])
  if (missing(name)){
    name = colnames(x[[1]])
  }
  if (number_of_chain<=4){
    mfrow_dim = c(1,number_of_chain)
    paper_dim = list(width=4*number_of_chain*2, height=4)
  } else if (number_of_chain<=6){
      mfrow_dim = c(3,3)
      paper_dim = list(width=4*4*3, height=3*4)
  }else if (number_of_chain<=8){
    mfrow_dim = c(2,4)
    paper_dim = list(width=4*4*2, height=2*4)
  } else {
    warning("too many chians, only the first 8 chains is plotted")
    mfrow_dim = c(2,4)
    number_of_chain=8
    paper_dim = c(width=4*4*2, height=2*4)
  }

  pdf(file=file, width=paper_dim$width, height=paper_dim$height)
  for (i in 1:number_of_variable){
    par(mfrow=mfrow_dim)
    min_x = min(unlist( lapply(x, function(z) z[,i ])))
    max_x = max(unlist( lapply(x, function(z) z[,i ])))
    for (j in 1:number_of_chain){
      xx = x[[j]][,i]
      t = mean(xx)/sd(xx)
      col = "black"
      if (abs(t) > 1.50) col = "yellow"
      if (abs(t) > 1.64) col = "blue"
      if (abs(t) > 1.96) col = "green"
      plot(xx, type='l' ,col=col , ylim=c(min_x,max_x))
      abline(h=0,col="dark red")
      abline(h=mean(xx),col="red", lwd=2)
      abline(h=mean(xx)+sd(xx)*c(1,-1),col="red",lty=2)
      title(main=name[i] , sub=paste0(round(mean(xx),2),"(" , round(t,2),")" ) )
    }
  }
  dev.off()

}


#' save.mcmc
#' @name save.mcmc
#' @aliases save.mcmc
#' @title save.mcmc
#' @details save.mcmc
#' @param x x
#' @param i i 
#' @param foldername foldername
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
save.mcmc = function(x, i , foldername){
  if (missing(i)){
    i = x$index 
  }
  if (!is.null(x$index)){
    i = x$index
  }

  obj_name = "out"%+%i
  assign(obj_name, x)
  # save(list=obj_name , file=foldername %+% "last.rData")
  save(list=obj_name , file=foldername %+% format(Sys.time(), "%y%m%d%H%M%S")%+%".rData")
  rm(list = "obj_name")
}

#' loadMergeMCMC
#' @name loadMergeMCMC
#' @aliases loadMergeMCMC
#' @title loadMergeMCMC
#' @details loadMergeMCMC
#' @param max max number of merge
#' @param folder_exist vector of names of folder 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
loadMergeMCMC = function(max=40,folder_exist){
  if (missing(folder_exist))
    folder_exist = checkFolder()

  number_of_chain = length(folder_exist)

  chain_list = vector("list", number_of_chain)

  for (i in 1:number_of_chain) {
    chain_list[[i]] = readAllObject(folder_exist[i] %+% "/", except="last.rData")
    ordering = unlist(sapply(chain_list[[i]], "[[", "index"))
    if (is.null(unlist(ordering))){
      ordering = as.numeric(gsub("out","",names(ordering)))
    }

    chain_list[[i]] = chain_list[[i]][order(ordering)]
    if (length(chain_list[[i]])>40)
      chain_list[[i]]= tail(chain_list[[i]] ,max)

    print(length(chain_list[[i]]))

    if (length(chain_list[[i]])>1){
      chain_list[[i]] = do.call(merge, chain_list[[i]])
    } else {
      chain_list[[i]] = chain_list[[i]][[1]]
    }


  }

  chain_list
}

#' Merge and delete all object in all folder
#' @name mergeAndDelete
#' @aliases mergeAndDelete
#' @title mergeAndDelete
#' @details mergeAndDelete
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
mergeAndDelete = function(){
  folder_exist = checkFolder()

  number_of_chain = length(folder_exist)

  for (i in 1:number_of_chain) {
    cat(i," ... \n")
    path = folder_exist[i] %+% "/"

    env1 = new.env()
    file_name = setdiff(dir(path),"last.rData")

    sapply(path %+% file_name, load, envir = env1 )

    chain_list = as.list(env1)

    ordering = unlist(sapply(chain_list, "[[", "index"))
    if (is.null(unlist(ordering))){
      ordering = as.numeric(gsub("out","",names(ordering)))
    }

    chain_list = chain_list[order(ordering)]

    if (length(chain_list)>1){
      chain_list = do.call(merge, chain_list)
    } else {
      chain_list = chain_list[[1]]
    }
    obj_name = "out"%+%max(ordering)
    assign(obj_name, chain_list)

    save(list=obj_name , file=path %+% format(Sys.time(), "%y%m%d%H%M%S")%+%".rData")
    file.remove(path %+% file_name)
  }
}

#' Remover file with size zero(corrupted files)
#' @name removeFileSizeZero
#' @aliases removeFileSizeZero
#' @title removeFileSizeZero
#' @details removeFileSizeZero
#' @param folder_name folder_name
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
removeFileSizeZero = function (folder_name){
  folder_name=as.character(folder_name)
  file_name = as.character(folder_name) %+% "/" %+%  dir(folder_name)
  file_with_zero_size =file_name[ which( file.info(file_name)$size == 0 ) ] 
  cat("\nRemoving:\n")
  cat(paste(file_with_zero_size,collapse='\n'))
  file.remove(file_with_zero_size)
}

#' computeSummaryTableFromDrive
#' @name computeSummaryTableFromDrive
#' @aliases computeSummaryTableFromDrive
#' @title computeSummaryTableFromDrive
#' @details computeSummaryTableFromDrive
#' @param path path
#' @param tail tail
#' @param ... ...
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
computeSummaryTableFromDrive = function(path, tail,...){
  path =  as.character(path)
  file_name = path %+% "/" %+% dir(path)
  out_list_mean = vector("list",length(file_name))
  out_list_mean_square = vector("list",length(file_name))
  m_vector = rep(NA, length(file_name))
  for (i in 1:length(file_name)){
    env1 = new.env()
    load(file_name[i], envir = env1)
    env1 = as.list(env1)
    tmp = getParameterMatrix( env1[[1]] , tail,...)
    out_list_mean[[i]] = colMeans(tmp)
    out_list_mean_square[[i]] = colMeans(tmp^2)
    m_vector[i] = env1[[1]]$m
  }
  mean_matrix = do.call(rbind, out_list_mean)
  mean_square_matrix = do.call(rbind, out_list_mean_square)

  name = colnames(mean_matrix)
  sum_m = sum(m_vector)
  mean_parameter = as.vector(m_vector %*% mean_matrix )  / sum_m
  square_mean_parameter = as.vector(m_vector %*% mean_square_matrix)   / sum_m
  names(mean_parameter) = name
  names(square_mean_parameter) = name

  var_parameter = square_mean_parameter - mean_parameter^2

  data.frame(mean=mean_parameter, var=var_parameter, sd= sqrt(var_parameter), m = sum_m)
} 

#' compareMCMC2
#' @name compareMCMC2
#' @aliases compareMCMC2
#' @title compareMCMC2
#' @details compareMCMC2
#' @param x x
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
compareMCMC2 = function(x){
  name = rownames(x[[1]])
  m = foreach(i=x, .combine=rbind) %do% mean(i$m)
  w = m / sum(m)
  between_var = diag( var(t(foreach(i=x, .combine=cbind) %do% i$mean)) )
  within_var = as.vector(foreach(i=x, .combine=cbind) %do% i$var %*% w )
  total_var = within_var + between_var
  total_mean = as.vector(foreach(i=x, .combine=cbind) %do% i$mean %*% w )
  var_ratio = between_var/within_var
  between_mean_min = apply(foreach(i=x, .combine=cbind) %do% i$mean ,1,min)
  between_mean_max = apply(foreach(i=x, .combine=cbind) %do% i$mean ,1,max)

  converge1 = var_ratio < 0.2
  
  var_con  = ifelse(var_ratio < 0.1, "!!!", ifelse(var_ratio < 0.2, "!!", ifelse(var_ratio < 0.3, "!", "")))


  summary_table_simple = cbind(generateSignificance(cbind(total_mean,sqrt(total_var)), name), 
    min=between_mean_min,
    max=between_mean_max,
    var_ratio=var_ratio,
    var_con=var_con
  )
  i=NULL
  not_converge = foreach(i = x) %do% generateSignificance(i[,c("mean","sd")], name)[!converge1,]

  not_converge_summary_table_simple = summary_table_simple[!converge1,]
  list(summary_table_simple=summary_table_simple, not_converge=not_converge,not_converge_summary_table_simple=not_converge_summary_table_simple)
}


# out = lapply(1:8, computeSummaryTableFromDrive )
# tab1 = compareMCMC2(out)

# ## between group
# loaded_files = loadMergeMCMC(100)
# tab2 = compareMCMC(loaded_files, 500000)

# ## within group


