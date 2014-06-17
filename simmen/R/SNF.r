#' loglikelihood_SNF
#' @name loglikelihood_SNF
#' @aliases loglikelihood_SNF
#' @title loglikelihood_SNF
#' @param y indicator of whether i and j is connected
#' @param x1 variables of i 
#' @param x2 variables of j
#' @param delta parameters
#' @param y_not complement of y
#' @return value of log likelihood
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
loglikelihood_SNF = function(y, x1, x2, delta, y_not){
  if (missing(y_not))
    y_not = !y
  p = pnorm(x1 %*% delta) * pnorm(x2 %*% delta)
  out = sum(log(p^y*(1-p)^y_not))

  if (!is.finite(out)){
    return(-1e+20)
  }
  return(out)
}

#' lik_grad_single_SNF
#' @name lik_grad_single_SNF
#' @aliases lik_grad_single_SNF
#' @title lik_grad_single_SNF
#' @param y indicator of whether i and j is connected
#' @param x1 variables of i 
#' @param x2 variables of j
#' @param delta parameters
#' @param y_not complement of y
#' @return value of gradient of log likelihood
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
lik_grad_single_SNF = function(y, x1, x2, delta, y_not){
  R1 = x1 %*% delta
  R2 = x2 %*% delta

  pi = pnorm(R1)
  pj = pnorm(R2)
  di = dnorm(R1)
  dj = dnorm(R2)
  p = pi * pj

  f =  (y-p) / (p* (1-p)) 
  f[is.nan(f)] = 0

  out =  as.vector( f  ) * (as.vector(pj * di) * x1  + as.vector(pi *dj) * x2 )
  as.vector(colSums(out))
}

#' drawYstar
#' @name drawYstar
#' @aliases drawYstar
#' @title drawYstar
#' @param y indicator of whether i and j is connected
#' @param ystar_other latent value of j
#' @param mean x*delta
#' @param y_not complement of y
#' @param sd sd of the error
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
drawYstar = function(y, ystar_other, mean, y_not=!y, sd=1){

  ystar_other_positive = ystar_other>=0  

  index_case1 = y
  index_case2 = as.logical(y_not * ystar_other_positive)
  index_case3 = as.logical(y_not * !ystar_other_positive)

  n1=sum(index_case1)
  n2=sum(index_case2)
  n3=sum(index_case3)
  
  ystar_new = rep(NA, length(y))

  if (n1>0)
    ystar_new[index_case1] = rtruncnorm(1,a=0,b=Inf, mean=mean[index_case1],sd=sd)  
  if (n2>0)
    ystar_new[index_case2] =rtruncnorm(1,a=-Inf,b=0,mean=mean[index_case2],sd=sd) 
  if (n3>0)
    ystar_new[index_case3] = mean[index_case3] +rnorm(n3,sd=sd)

  ystar_new
}




#' Strategy Network Formation
#' @name SNF
#' @rdname SNF
#' @aliases SNF
#' @aliases SNF.static.maxLik
#' @aliases SNF.static.mcmc
#' @aliases SNF.dynamic.mcmc
#' @title SNF
#' @param data data
#' @param method Estimation method, either "static.maxLik","static.mcmc","dynamic.mcmc". Default is "static.maxLik"
#' @param m m
#' @param last_estimation last_estimation
#' @param update_tau update_tau
#' @param tau tau 
#' @param ... others argument.
#' @return SNF object
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
SNF = function(data , method=c("static.maxLik","static.mcmc","dynamic.mcmc"), ...){
  method = match.arg(method)
  switch(method,
    static.maxLik = SNF.static.maxLik(data,...),
    static.mcmc = SNF.static.mcmc(data,...),
    dynamic.mcmc = SNF.dynamic.mcmc(data,...),
  )
}

#' @rdname SNF
#' @export
SNF.static.maxLik = function(data,...){
  tic()

  self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
  friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))

  y = do.call(rbind, lapply(data, function(z) z$response_self))

  number_of_network = ncol(y)
  network_name = data[[1]]$network_name

  out = vector("list",number_of_network)
  summary_table = vector("list",number_of_network)

  for (i in 1:number_of_network){
    yy = y[,i]
    yy_not = !yy
    start = glm(yy~self_data_matrix-1, family=binomial(link="probit"))$coef

    out[[i]] = maxLik(function(z, ...) loglikelihood_SNF(delta=z, ...) , start=start ,  x1=self_data_matrix, x2=friends_data_matrix, y=yy, y_not=yy_not , grad= lik_grad_single_SNF, method="BFGS")

      summary_table[[i]] = generateSignificance(summary(out[[i]])$estimate[,1:2])
      rownames(summary_table[[i]]) = network_name[i] %+% "_" %+%colnames(self_data_matrix)
  }

  summary_table = do.call(rbind,summary_table)
  toc()
  out2 = list(out=out, summary_table=summary_table)
  class(out2) = "SNF.static.maxLik"
  out2
}



#' @rdname SNF
#' @export
SNF.static.mcmc = function(data, m=1000, last_estimation,...){

  self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
  friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))
  
  y = do.call(rbind, lapply(data,"[[","response_self"))
  y_not = !y
  number_of_network = ncol(y)


  name = colnames(self_data_matrix)
  k =  ncol(self_data_matrix)
  n = NROW(y)*2
  delta_matrix = rep(list(matrix(0, nrow=m, ncol= k )), number_of_network)

  if (number_of_network>1){
    number_col_Sigma_matrix = number_of_network*(number_of_network-1)/2
    Sigma_matrix = matrix(0, nrow=m, ncol = number_col_Sigma_matrix )

    sigma_name = genPairwiseIndex(number_of_network)
    sigma_name = sigma_name[,1] %+% sigma_name[,2]
    
    colnames(Sigma_matrix) = "Sigma_" %+%  sigma_name
  } else{
    number_col_Sigma_matrix=1
    Sigma_matrix = matrix(0, nrow=m, ncol = number_col_Sigma_matrix )
    colnames(Sigma_matrix) = "Sigma_11" 
  }

  network_name = data[[1]]$network_name

  for (i in 1:number_of_network){
    colnames(delta_matrix[[i]]) = network_name[[i]] %+% "_" %+% name
  }

  ystar1 = matrix(0, nrow=nrow(y), ncol=ncol(y))
  ystar2 = matrix(0, nrow=nrow(y), ncol=ncol(y))
  Sigma = matrix(0.5,number_of_network,number_of_network)
  diag(Sigma) = 1
  delta = matrix(0, nrow=k, ncol=number_of_network)


  if (!missing(last_estimation)){
    ystar1 = last_estimation$ystar1
    ystar2 = last_estimation$ystar2
    delta = last_estimation$delta
    Sigma = last_estimation$Sigma
  }


  X = rbind(self_data_matrix, friends_data_matrix)
  XX_inv = solve(crossprod(X))


  xb1 = self_data_matrix %*% delta 
  xb2 = friends_data_matrix %*% delta 


  tic()
  for (i in 1:m){
    if (i %% 1000 == 0 )
      cat(i, ">\n")

    ## update ystar
      for( j in 1:number_of_network){
        ystar1_demean = ystar1 - xb1
        ystar2_demean = ystar2 - xb2

        temp = find_normal_conditional_dist(a=ystar1_demean, i=j, j=-j, Sigma=Sigma)
        ystar1[,j] = drawYstar(y=y[,j] , ystar_other=ystar2[,j], mean=xb1[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )

        temp = find_normal_conditional_dist(a= ystar2_demean, i=j, j=-j, Sigma=Sigma)
        ystar2[,j] = drawYstar(y=y[,j] , ystar_other=ystar1[,j], mean=xb2[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )
      }
      ystar1_demean = ystar1 - xb1
      ystar2_demean = ystar2 - xb2

      ystar_demean = rbind(ystar1_demean,ystar2_demean)
      ystar = rbind(ystar1,ystar2)

      for ( j in 1:number_of_network){
        temp = find_normal_conditional_dist(a=ystar_demean, i=j, j=-j, Sigma=Sigma)

        beta_coef = XX_inv %*% crossprod(X, (ystar[,j]-temp$mean ) )

        delta[,j] = mvrnorm(n=1, mu=beta_coef, XX_inv * as.vector(temp$var) )

        ystar_demean = ystar[,j] - X %*% delta[,j]

        delta_matrix[[j]][i,] = delta[,j]
      }
      xb1 = self_data_matrix %*% delta 
      xb2 = friends_data_matrix %*% delta 

      ystar1_demean = ystar1 - xb1
      ystar2_demean = ystar2 - xb2
      ystar_demean = rbind(ystar1_demean,ystar2_demean)


      ## Sigma
      if (number_of_network > 1 ){
        Sigma = solve( rwish(n , solve( crossprod(ystar_demean )) ) )

        normalization = diag(1/sqrt(diag(Sigma)))

        Sigma = normalization %*%  Sigma %*% t(normalization) 
        Sigma_matrix[i,] = Sigma[lower.tri(Sigma)]
      } else {
        Sigma = as.matrix(1)
        Sigma_matrix[i,] = 1
      }
  }   
  toc()

  out = list(delta_matrix=delta_matrix, ystar1=ystar1,ystar2=ystar2, Sigma=Sigma,Sigma_matrix=Sigma_matrix, delta=delta)

  class(out) = "network_formation.mcmc" 
  out
}

#' merge.SNF.static.mcmc
#' @name merge.SNF.static.mcmc
#' @aliases merge.SNF.static.mcmc
#' @title merge.SNF.static.mcmc
#' @param x First object to merge with
#' @param y Second object to merge with 
#' @param ... not used
#' @return A new SNF.static.mcmc object
#' @method merge SNF.static.mcmc
#' @export merge SNF.static.mcmc
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
merge.SNF.static.mcmc = function(x,y,...){
  out = y
  if (is.list(x$delta_matrix)){
    for (i in 1:length(x$delta_matrix)){
      out$delta_matrix[[i]] = rbind(x$delta_matrix[[i]], y$delta_matrix[[i]] )
    }
    out$Sigma_matrix = rbind(x$Sigma_matrix, y$Sigma_matrix)
  } else{
    out$delta_matrix = rbind(x$delta_matrix , y$delta_matrix) 
  }
  out
}

#' Get a matrix of parameter 
#' @name getParameterMatrix.SNF.static.mcmc
#' @aliases getParameterMatrix.SNF.static.mcmc
#' @title getParameterMatrix.SNF.static.mcmc
#' @param x SNF.static.mcmc
#' @param tail iteration to be used. Negative value: Removing the first \code{tail} iterations. Positive value: keep the last \code{tail} iterations. If -1< code{tail}< 1, it represent the percentage of iterations.
#'' @param ... not used
#' @return A matrix
#' @method getParameterMatrix SNF.static.mcmc
#' @export getParameterMatrix SNF.static.mcmc
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
getParameterMatrix.SNF.static.mcmc = function(x, tail, ...){
  if (is.list(x$delta_matrix)){
    out = do.call(cbind, x$delta_matrix)
    out = cbind(out, x$Sigma_matrix )
  } else{
    out = x$delta_matrix
  }
  if (!missing(tail)) {
    out = extractTail(out, tail)
  }
  out
}


#' Create a summary table
#' @name summary.SNF.static.mcmc
#' @aliases summary.SNF.static.mcmc
#' @title summary.SNF.static.mcmc
#' @param object SNF.static.mcmc object
#' @param ... tail:  iteration to be used. Negative value: Removing the first \code{tail} iterations. Positive value: keep the last \code{tail} iterations. If -1< code{tail}< 1, it represent the percentage of iterations.
#' @return A summary table
#' @method summary SNF.static.mcmc
#' @export summary SNF.static.mcmc
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
summary.SNF.static.mcmc = function(object,...){
  computeSummaryTable(object,...)
}

#' Create a summary table
#' @name summary.SNF.static.maxLik
#' @aliases summary.SNF.static.maxLik
#' @title summary.SNF.static.maxLik
#' @param object SNF.static.maxLik object
#' @param ... not used
#' @return A summary table
#' @method summary SNF.static.maxLik
#' @export summary SNF.static.maxLik
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
summary.SNF.static.maxLik = function(object,...){
  object$summary_table
}


# single_network_formation_mcmc = function(data,  m=1000, last_estimation){
#   self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
#   friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))

#   response = do.call(rbind, lapply(data, function(z) z$response_self))
#   response_not = !response

#   name = colnames(self_data_matrix)
#   k =  ncol(self_data_matrix)
#   delta_matrix = matrix(0, nrow=m+1, ncol=k)
#   ystar1 = rep(0, length(response))
#   ystar2 = rep(0, length(response))
#   network_name = data[[1]]$network_name

#   if (!missing(last_estimation)){
#     delta_matrix[1,] = tail(last_estimation$delta_matrix,1)
#     ystar1 = last_estimation$ystar1
#     ystar2 = last_estimation$ystar2
#   }

#   colnames(delta_matrix) = network_name %+% "_" %+% colnames(self_data_matrix)

#   X = rbind(self_data_matrix, friends_data_matrix)
#   XX_inv = solve(crossprod(X))

#   tic()
#   for (i in 1:m){
#     if (i %% 1000 == 0 ){
#       cat(i ,">\n")
#     }

#     xb1 = self_data_matrix %*% delta_matrix[i, ]
#     xb2 = friends_data_matrix %*% delta_matrix[i, ]

#     ystar1 = drawYstar(y=response , ystar_other=ystar2, mean=xb1, y_not= response_not)
#     ystar2 = drawYstar(y=response, ystar_other=ystar1, mean=xb2, y_not= response_not)

#     delta_matrix[i+1, ] =   mvrnorm(n=1, mu=XX_inv %*% crossprod(X,c(ystar1,ystar2)), XX_inv)

#   }
#   toc()
#   delta_matrix = tail(delta_matrix,-1)

#   out = list(delta_matrix=delta_matrix, ystar1=ystar1, ystar2=ystar2)
#   class(out) = "network_formation" 
#   out
# }





# lik_single_network_formation_parser = function(data, delta, network_id=1){
#       loglikelihood_network_formation(
#         x1=data$self_data_matrix, 
#         x2= data$friends_data_matrix, 
#         y=data$response1, 
#         delta

#   if (network_id==1){
#     return( 
#       )
#     )
#   } else if (network_id==2){
# return( 
#       loglikelihood_network_formation(
#         x1=data$self_data_matrix, 
#         x2= data$friends_data_matrix, 
#         y=data$response2, 
#         delta
#       )
#     )  }
# })

# lik_single_network_formation_par = function(cl, delta, network_id=1, G=5){
#   sum(
#     parSapply(cl, 1:G, 
#       function(i,delta,network_id) lik_single_network_formation_parser(
#           data[[i]], 
#           delta=delta, 
#           network_id=network_id
#         ),
#       delta=delta,
#       network_id=network_id
#     ) 
#   )
# })


# lik_grad_single_network_formation_parser = function(data, delta, network_id=1){

#   if (network_id==1){
#     return( 
#       lik_grad_single_network_formation(
#         x1=data$self_data_matrix, 
#         x2= data$friends_data_matrix, 
#         y=data$response1, 
#         delta
#       )
#     )
#   } else if (network_id==2){
#     return( 
#       lik_grad_single_network_formation(
#         x1=data$self_data_matrix, 
#         x2= data$friends_data_matrix, 
#         y=data$response2, 
#         delta
#       )
#     )  }
# })

# lik_grad_single_network_formation_par = function(cl, delta, network_id=1, G=5){
#   rowSums(
#     parSapply(cl, 1:G, 
#       function(i,delta,network_id) lik_grad_single_network_formation_parser(
#           data[[i]], 
#           delta=delta, 
#           network_id=network_id
#         ), 
#       delta=delta,
#       network_id=network_id
#     ) 
#   )
# })




# single_network_formation = function(data, network_id=1){
#   tic()
#   require("maxLik")

#   self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
#   friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))
#   if (network_id==1){
#     response = unlist(lapply(data, function(z) z$response1))
#   } else if (network_id==2){
#     response = unlist(lapply(data, function(z) z$response2))
#   }

#   start = rep(0,ncol(self_data_matrix))

#   system.time({
#   out= maxLik(function(z, ...) loglikelihood_network_formation(delta=z, ...) , start=start ,  self_data_matrix=self_data_matrix, friends_data_matrix=friends_data_matrix, response=response , grad= lik_grad_single_network_formation, method="BFGS")
#   })

#   summary_table = generateSignificance(summary(out)$estimate[,1:2])
#   rownames(summary_table) = colnames(self_data_matrix)

#   toc()
#   list(out, summary_table)
# })



# single_network_formation_parallel = function(data, cl, network_id=1){
#   tic()
#   name = colnames(data[[1]]$self_data_matrix)

#   start = rep(0,length(name))

#   out= maxLik(
#     lik_single_network_formation_par, 
#     start=start ,  
#     cl=cl, 
#     G=length(data),
#     network_id=network_id , 
#     grad= lik_grad_single_network_formation_par, 
#     method="BFGS"
#   )

#   summary_table = summary(out)$estimate
#   rownames(summary_table) = name
#   toc()
#   list(maxLik_object=out, summary_table=generateSignificance(summary_table[,1:2]))
# })


# single_network_formation_mcmc_v1 = function(start, tau, m, cl, network_id, G){
#   k = length(start)
#   delta_matrix = matrix(0, nrow=m+1, ncol=k)
#   delta_matrix[1, ] = start
#   update_rate=0

#   for (i in 1:m){
#     metro_obj = 
#       metropolis2(
#         beta_previous=delta_matrix[i,],
#         tau=tau,
#         likelihoodFunction=lik_single_network_formation_par,
#         cl=cl,
#         network_id=network_id,
#         G=G
#       )
#     delta_matrix[i+1,] = metro_obj$beta
#     update_rate = update_rate + metro_obj$update
#   }

#   delta_matrix = tail(delta_matrix,-1)
#   update_rate =update_rate / m
#   next_tau = tau * ifelse(update_rate==0,0.1,update_rate) / 0.27 

#   return(list(delta_matrix = delta_matrix, update_rate =update_rate, next_parameter = tail(delta_matrix,1), tau=tau, next_tau = next_tau))
# })

# single_network_formation_mcmc_v2 = function(start, tau, m, cl, network_id, G){
#   k = length(start)
#   delta_matrix = matrix(0, nrow=m+1, ncol=k)
#   delta_matrix[1, ] = start
#   update_rate=0

#   for (i in 1:m){
#     metro_obj = 
#       metropolis(
#         beta_previous=delta_matrix[i,],
#         tau=tau,
#         likelihoodFunction=lik_single_network_formation_par,
#         cl=cl,
#         network_id=network_id,
#         G=G
#       )
#     delta_matrix[i+1,] = metro_obj$beta
#     update_rate = update_rate + metro_obj$update
#   }

#   delta_matrix = tail(delta_matrix,-1)
#   update_rate =update_rate / m
#   next_tau = tau * ifelse(update_rate==0,0.1,update_rate) / 0.27 

#   return(list(delta_matrix = delta_matrix, update_rate =update_rate, next_parameter = tail(delta_matrix,1), tau=tau, next_tau = next_tau))
# })

## method 1 : update delta as vector
## method 2 : udpate delta one by one. Method 2 is more efficient, because it reduces the call to the likelihood function by half.


# single_network_formation_mcmc_RE = function(data, network_id, m=1000, last_estimation){
#   self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
#   friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))

#   response = as.logical(unlist(lapply(data, function(z) z$response[[network_id]])))
#   response_not = !response
  
#   n = sapply(data, function(z) length(z$y))
#   n2 = sapply(data, function(z) length(z$response[[network_id]]))

#   name = colnames(self_data_matrix)
#   k =  ncol(self_data_matrix)
#   delta_matrix = matrix(0, nrow=m+1, ncol=k)
#   sigma2e_matrix = matrix(1, nrow=m+1, ncol=1)
#   ystar1 = rep(0, length(response))
#   ystar2 = rep(0, length(response))
#   e = rep(0,sum(n)) #rnorm(sum(n))

#   if (!missing(last_estimation)){
#     delta_matrix[1,] = tail(last_estimation$delta_matrix,1)
#     sigma2e_matrix[1,] = tail(last_estimation$sigma2e_matrix,1)
#     ystar1 = last_estimation$ystar1
#     ystar2 = last_estimation$ystar2
#     e = last_estimation$e
#   }


#   colnames(delta_matrix) = colnames(self_data_matrix)

#   X = rbind(self_data_matrix, friends_data_matrix)
#   XX_inv = solve(crossprod(X))

#   full_group_index = genFullGroupIndex(data)
#   full_position_index = genFullPositionIndex(data)
#   full_position_matrix = genFullPositionMatrix(data)
#   row_sums_full_position_matrix = rowSums(full_position_matrix)

#   tic()
#   for (i in 1:m){
#     if (i %% 1000 == 0 ){
#       cat(i ,">\n")
#     }
#     full_e = genFulle(e,full_group_index )

#     xb1 = self_data_matrix %*% delta_matrix[i, ] + full_e$e_i
#     xb2 = friends_data_matrix %*% delta_matrix[i, ] + full_e$e_j

#     # update ystar
#     ystar1 = drawYstar(y=response , ystar_other=ystar2, mean=xb1, y_not= response_not)
#     ystar2 = drawYstar(y=response, ystar_other=ystar1, mean=xb2, y_not= response_not)
    
#     # update delta
#     mu_delta = XX_inv %*% crossprod(X,c(ystar1-full_e$e_i,ystar2-full_e$e_j))
#     delta_matrix[i+1, ] =   mvrnorm(n=1, mu=mu_delta, XX_inv)

#     # update e
#     # actually i dont need previous e, just need ystar1-xb1, ystar2-xb2 and the correct position.
#     #  
#     residual = c(ystar1, ystar2) - X %*% delta_matrix[i+1,]
#     # mean of residual by individual 

#     var_e = 1 / (row_sums_full_position_matrix + sigma2e_matrix[i,])
#     mean_e = as.numeric( full_position_matrix %*% residual  ) * var_e
#     e<-rnorm(sum(n),mean_e,sqrt(var_e))

#     sigma2e_matrix[i+1,]<-1/rgamma(1,length(e)/2,crossprod(e,e)/2) 
#   }
#   toc()

#   delta_matrix = tail(delta_matrix,-1)
#   sigma2e_matrix= tail(sigma2e_matrix,-1)

#   out = list(ystar1=ystar1, ystar2=ystar2,e=e,delta_matrix=delta_matrix, sigma2e_matrix=sigma2e_matrix)

#   class(out) = "single_network_formation_RE" 
#   out
# })


# single_network_formation_mcmc_parallel = function(data,cl, network_id, m=1000, last_estimation){

#   name = colnames(data[[1]]$self_data_matrix)
#   k =  ncol(data[[1]]$self_data_matrix)
#   n2 = sapply(data,function(z) length(z$response1)) 
#   G = length(data)
#   delta_matrix = matrix(0, nrow=m+1, ncol=k)
#   ystar1 = lapply(n2, rep,x=0 )

#   ystar2 = lapply(n2, rep,x=0 )

#   if (!missing(last_estimation)){
#     delta_matrix[1,] = tail(last_estimation$delta_matrix,1)
#     ystar1 = tail(last_estimation$ystar1)
#     ystar2 = tail(last_estimation$ystar2)
#   }

#   colnames(delta_matrix) = name

#   X = rbind(
#     do.call(rbind,lapply(data,"[[",i="self_data_matrix")), 
#     do.call(rbind,lapply(data,"[[",i="friends_data_matrix"))
#   )

#   XX_inv = solve(crossprod(X))

#   tic()
#   for (i in 1:m){
#     ystar1=
#     parLapply(cl, 1:G, function(z, ystar2, network_id,delta) {
#       drawYstar(
#         y= data[[z]]$response[[network_id]],
#         ystar_other = ystar2[[z]],
#         mean = data[[z]]$self_data_matrix %*% delta
#       )},
#       delta = delta_matrix[i,],
#       network_id = network_id,
#       ystar2=ystar2
#     )

#     ystar2=
#     parLapply(cl, 1:G, function(z, ystar1, network_id, delta) {
#       drawYstar(
#         y= data[[z]]$response[[network_id]],
#         ystar_other = ystar1[[z]],
#         mean = data[[z]]$friends_data_matrix %*% delta
#       )},
#       delta = delta_matrix[i,],
#       network_id = network_id,
#       ystar1=ystar1
#     )

#     delta_matrix[i+1, ] =   mvrnorm(n=1, mu=XX_inv %*% crossprod(X,c(unlist(ystar1), unlist(ystar2))), XX_inv)

#   }
#   toc()
#   delta_matrix = tail(delta_matrix,-1)

#   plotmcmc(delta_matrix,remove=remove)
#   print(computeSummaryTable(delta_matrix, remove=remove))

#   list(delta_matrix=delta_matrix, ystar1=ystar1, ystar2=ystar2)
# })




# drawYstar_multi = function(y, ystar_other, demean_ystar_corr, mean , y_not=!y, rho){

#   mean = mean + rho * (demean_ystar_corr) 
#   sd = sqrt(1-rho^2)

#   ystar_other_positive = ystar_other>=0  

#   index_case1 = y
#   index_case2 = as.logical(y_not * ystar_other_positive)
#   index_case3 = as.logical(y_not * !ystar_other_positive)

#   n = length(y)
#   n1=sum(index_case1)
#   n2=sum(index_case2)
#   n3=sum(index_case3)
  
#   stopifnot(n==n1+n2+n3)

#   ystar_new = rep(NA, length(y))

#   if (n1>0)
#     ystar_new[index_case1] = rtruncnorm(1,a=0,b=Inf, mean=mean[index_case1], sd=sd)  
#   if (n2>0)
#     ystar_new[index_case2] =rtruncnorm(1,a=-Inf,b=0,mean=mean[index_case2], sd=sd) 
#   if (n3>0)
#     ystar_new[index_case3] = mean[index_case3] +rnorm(n3, sd=sd)

#     stopifnot(!any(is.nan(ystar_new)))
#   ystar_new
# })



# multi_network_formation_mcmc_RE = function(data, m=1000, last_estimation){

#   self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
#   friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))
#   y1 = as.logical(unlist(lapply(data, function(z) z$response1)))
#   y2 = as.logical(unlist(lapply(data, function(z) z$response2)))
#   y1_not = !y1
#   y2_not = !y2

#   name = colnames(self_data_matrix)
#   k =  ncol(self_data_matrix)
#   n = length(y1)
#   delta1_matrix = matrix(0, nrow=m+1, ncol=k)
#   delta2_matrix = matrix(0, nrow=m+1, ncol=k)
#   rho_matrix = matrix(0, nrow=m+1,ncol=1)
#   ystar11 = rep(0, length(y1))
#   ystar12 = rep(0, length(y1))
#   ystar21 = rep(0, length(y2))
#   ystar22 = rep(0, length(y2))
#   e1 = rep(0,sum(sapply(data, "[[", "n")))
#   e2 = rep(0,sum(sapply(data, "[[", "n")))
#   sigma2e1_matrix = matrix(1, nrow=m+1,ncol=1)
#   sigma2e2_matrix = matrix(1, nrow=m+1,ncol=1)

#   colnames(delta1_matrix) = name
#   colnames(delta2_matrix) = name

#   if (!missing(last_estimation)){
#     delta1_matrix[1,] = tail(last_estimation$delta1,1)
#     delta2_matrix[1,] = tail(last_estimation$delta2,1)
#     rho_matrix[1,] = tail(last_estimation$rho,1)
#     ystar11 = last_estimation$ystar11
#     ystar12 = last_estimation$ystar12
#     ystar21 = last_estimation$ystar21
#     ystar22 = last_estimation$ystar22
#     e1 = last_estimation$e1
#     e2 = last_estimation$e2
#   }



#   X = rbind(self_data_matrix, friends_data_matrix)
#   XX_inv = solve(crossprod(X))

#   # xb11 = self_data_matrix %*% delta1_matrix[1, ]
#   # xb12 = friends_data_matrix %*% delta1_matrix[1, ]
#   # xb21 = self_data_matrix %*% delta2_matrix[1, ]
#   # xb22 = friends_data_matrix %*% delta2_matrix[1, ]

#   full_group_index = genFullGroupIndex(data)
#   full_position_index = genFullPositionIndex(data)
#   full_position_matrix = genFullPositionMatrix(data)
#   row_sums_full_position_matrix = rowSums(full_position_matrix)

#   tic()
#   for (i in 1:m){
#     if (i %% 1000 == 0 )
#       cat(i, ">\n")
#     rho = rho_matrix[i, 1]

#     full_e1 = genFulle(e1, full_group_index )
#     full_e2 = genFulle(e2, full_group_index )

#     xb11 = self_data_matrix %*% delta1_matrix[i, ] + full_e1$e_i
#     xb12 = friends_data_matrix %*% delta1_matrix[i, ] + full_e1$e_j
#     xb21 = self_data_matrix %*% delta2_matrix[i, ] + full_e2$e_i
#     xb22 = friends_data_matrix %*% delta2_matrix[i, ] + full_e2$e_j

#     ystar11 = drawYstar_multi(y=y1 , ystar_other=ystar12, demean_ystar_corr=ystar21-xb21 ,mean=xb11,  y_not= y1_not, rho=rho)
#     ystar12 = drawYstar_multi(y=y1, ystar_other=ystar11, demean_ystar_corr=ystar22-xb22, mean=xb12, y_not= y1_not, rho=rho)
#     ystar21 = drawYstar_multi(y=y2 , ystar_other=ystar22, demean_ystar_corr=ystar11-xb11, mean=xb21, y_not= y2_not, rho=rho)
#     ystar22 = drawYstar_multi(y=y2, ystar_other=ystar21, demean_ystar_corr=ystar12-xb12, mean=xb22, y_not= y2_not, rho=rho)

#     ystar1 = c(ystar11,ystar12)
#     ystar2 = c(ystar21,ystar22)

#     ystar2_demean = ystar2 - c(xb21,xb22)  - c(full_e2$e_i,full_e2$e_j)
#     new_y1 = ystar1 - rho*ystar2_demean - c(full_e1$e_i,full_e1$e_j)
#     lm1 = myFastLm(X, new_y1)
#     delta1_matrix[i+1, ] =   mvrnorm(n=1, mu=lm1$coef, (lm1$cov)/(lm1$s^2) * (1-rho^2) )
#     xb11 = self_data_matrix %*% delta1_matrix[i+1, ]
#     xb12 = friends_data_matrix %*% delta1_matrix[i+1, ]


#     ystar1_demean = ystar1 - c(xb11,xb12) - c(full_e1$e_i,full_e1$e_j)
#     new_y2 = ystar2 - rho*ystar1_demean - c(full_e2$e_i,full_e2$e_j)
#     lm2 = myFastLm(X, new_y2 )
#     delta2_matrix[i+1, ] =   mvrnorm(n=1, mu=lm2$coef, (lm2$cov)/(lm2$s^2) * (1-rho^2) )
#     xb21 = self_data_matrix %*% delta2_matrix[i+1, ]
#     xb22 = friends_data_matrix %*% delta2_matrix[i+1, ]

#     ystar2_demean= ystar2 - c(xb21,xb22)- c(full_e2$e_i,full_e2$e_j)
#     # update rho
#     # var_rho = (1-rho^2)/ sum((ystar1_demean)^2)
#     # mean_rho = var_rho / (1-rho^2) * crossprod(ystar1_demean, ystar2_demean  )

#     # mean_rho = cov(ystar1_demean,ystar2_demean)
#     # var_rho = (mean(ystar1_demean^2* ystar2_demean^2 ) - mean( ystar1_demean * ystar2_demean )^2) / 2/n

#     mean_rho = mean(ystar1_demean * ystar2_demean)
#     var_rho = var(ystar1_demean * ystar2_demean)/2/n

#     rho_matrix[i+1,1 ] = rtruncnorm(1,mean= mean_rho, sd= sqrt(var_rho),a=-.999,b=.999) 

#     # update e1
#     residual1 = ystar1 - rho*ystar2_demean - c(full_e1$e_i,full_e1$e_j) - c(xb11,xb12)
#     residual2 = ystar2 - rho*ystar1_demean - c(full_e2$e_i,full_e2$e_j) - c(xb21,xb22)


#     # mean of residual by individual 

#     var_e1 = 1 / (row_sums_full_position_matrix + sigma2e1_matrix[i,])
#     mean_e1 = as.numeric( full_position_matrix %*% residual1  ) * var_e1
#     e1<-rnorm(length(e1),mean_e1,sqrt(var_e1))

#     sigma2e1_matrix[i+1,]<-1/rgamma(1,length(e1)/2,crossprod(e1,e1)/2) 

#     var_e2 = 1 / (row_sums_full_position_matrix + sigma2e2_matrix[i,])
#     mean_e2 = as.numeric( full_position_matrix %*% residual2  ) * var_e2
#     e2<-rnorm(length(e2),mean_e2,sqrt(var_e2))

#     sigma2e2_matrix[i+1,]<-1/rgamma(1,length(e2)/2,crossprod(e2,e2)/2) 

#     # cat(i, "> ", rho_matrix[i+1,],"\n")
#   }   
#   toc()

#   delta1_matrix = tail(delta1_matrix,-1)
#   delta2_matrix = tail(delta2_matrix,-1)
#   rho_matrix = tail(rho_matrix,-1)
#   sigma2e1_matrix = tail(sigma2e1_matrix,-1)
#   sigma2e2_matrix = tail(sigma2e2_matrix,-1)
#   out = list(delta1_matrix=delta1_matrix, delta2_matrix=delta2_matrix, rho_matrix=rho_matrix, ystar11=ystar11, ystar12=ystar12, ystar21=ystar21, ystar22=ystar22, e1=e1, e2=e2, sigma2e1_matrix=sigma2e1_matrix, sigma2e2_matrix=sigma2e2_matrix)

#   class(out) = "multi_network_formation_RE" 
#   out
# })




# plotmcmc.single_network_formation_RE = function(x, tail=-0.2){
#   data_matrix = cbind(x$delta, x$sigma2e_matrix)
#   colnames(data_matrix) = c(colnames(x$delta), "Sigma2")
#   plotmcmc.default(data_matrix, tail=tail)
# })

# merge.single_network_formation_RE = function(x,y,...){
#   out = y
#   out$delta_matrix = rbind(x$delta_matrix , y$delta_matrix) 
#   out$sigma2e_matrix = rbind(x$sigma2e_matrix , y$sigma2e_matrix) 
#   out
# })


# getParameterMatrix.multi_network_formation = function(x){
#   out = do.call(cbind, x$delta_matrix)
#   out = cbind(out, x$Sigma_matrix )
#   out
# }

# merge.multi_network_formation = function(x,y){
#   out = y
#   for (i in 1:length(x$delta_matrix)){
#     out$delta_matrix[[i]] = rbind(x$delta_matrix[[i]], y$delta_matrix[[i]] )
#   }
#   out$Sigma_matrix = rbind(x$Sigma_matrix, y$Sigma_matrix)
#   out
# })

# plotmcmc.multi_network_formation_RE = function(x, tail=-0.2){
#   data_matrix = cbind(x$delta1_matrix, x$delta2_matrix, x$rho_matrix, x$sigma2e1_matrix, x$sigma2e2_matrix)
#   colnames(data_matrix) = c(colnames(x$delta1_matrix), colnames(x$delta2_matrix), "rho", "sigma2e1", "sigma2e2")

#   plotmcmc.default(data_matrix, tail=tail)
# })

# merge.multi_network_formation_RE = function(x,y){
#   out = y
#   out$delta1_matrix = rbind(x$delta1_matrix , y$delta1_matrix) 
#   out$delta2_matrix = rbind(x$delta2_matrix , y$delta2_matrix) 
#   out$rho_matrix = rbind(x$rho_matrix , y$rho_matrix) 
#   out$sigma2e1_matrix = rbind(x$sigma2e1_matrix , y$sigma2e1_matrix) 
#   out$sigma2e2_matrix = rbind(x$sigma2e2_matrix , y$sigma2e2_matrix) 
#   out
# })


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################






## Given a seq, beta, D, compute the likelihood
## Given a seq m  n*(n-1)/2 x 2 matrix
## D: n by n network matrix
## U_xb : an n by n utility matrix, i,j element is the utility of i to make friends with j. (Xbeta)
## delta1 delta2 

#' Draw random sample of meeting sequence
#' @name DrawSeqSample
#' @aliases DrawSeqSample
#' @title DrawSeqSample
#' @param x x
#' @param p p
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
DrawSeqSample = function(x , p =0.01){
	if (is.list(x)){
		out = vector("list",length(x))
		if (length(p)==1)
			p = rep(p,length(x))
		for (i in 1:length(x)){
			out[[i]] = DrawSeqSample(x[[i]], p[[i]])
		}
		return(out)
	}	else if (is.vector(x)){
		n = length(x)
		pn = pmax(2,ceiling(n * p)  )
		to_change = sample(n,pn)
		reorder = sample(to_change)

		x[to_change] = x[reorder]
		return(x)
	} else if (is.matrix(x)){
		n = nrow(x)
		pn = pmax(2,ceiling(n * p)  )
		to_change = sample(n,pn)
		reorder = sample(to_change)
		x[to_change,] = x[reorder,]
		return(x)	
	}
}

# repeat

# computeNetworkSummary_r = function(seq_m,D){
# 	n = nrow(D)
# 	D0 = matrix(0,ncol=n,nrow=n)
# 	degree = matrix(0,ncol=n,nrow=n)
# 	common_frds_1 = matrix(0,ncol=n,nrow=n)
# 	common_frds_2 = matrix(0,ncol=n,nrow=n)

#   nn = n*(n-1)/2
# 	for ( i in 1:nn){
# 		index = seq_m[i,]
# 		index1 = index[1]
# 		index2 = index[2]
# 		if (D[index1,index2]==1) {
# 			D0[index1,index2] = D0[index2,index1]= 1
# 		}
# 		d1 = D0[index1,]
# 		d2 = D0[index2,]
# 		degree[index1,index2] = sum(d1)
# 		degree[index2,index1] = sum(d2)
# 		common_frds_1[index1,index2] = sum(d1*d2)
# 	}

# 	lower_tri = lower.tri(degree)
#   degree1 = degree[ lower_tri ] 
# 	degree2 = t(degree)[ lower_tri ] 
# 	common_frds_1 = common_frds_1[ lower_tri ] 

# 	list(self=cbind(degree1,degree1^2,common_frds_1), friends=cbind(degree2,degree2^2,common_frds_1))
# })


# src=
# '
# arma::mat seq_m2 = Rcpp::as<arma::mat>(seq_m);
# arma::mat DD = Rcpp::as<arma::mat>(D);
# int nn = DD.n_rows;
# arma::mat D00 = arma::zeros(nn,nn);
# arma::mat degreee = arma::zeros(nn,nn);
# arma::mat common_frds_11 = arma::zeros(nn,nn);

# for ( int i=0 ; i<nn*(nn-1)/2; i++ ){
# 	int index1 = seq_m2(i,0) -1 ;
# 	int index2 = seq_m2(i,1) -1;
# 	if (DD(index1,index2)==1) {
# 		D00(index1,index2) = 1;
# 		D00(index2,index1) = 1;
# 	}

# 	degreee(index1,index2) = sum(D00.col(index1)) ; 
# 	degreee(index2,index1) = sum(D00.col(index2))  ;
# 	common_frds_11(index1,index2) = arma::as_scalar(D00.col(index1).t() * D00.col(index2)) ; 
# }

# return Rcpp::List::create(Rcpp::Named("degree")=degreee, Rcpp::Named("common_frds_1")=common_frds_11);
# '

# g <- cxxfunction(signature(seq_m="integer", D="integer"),
# 	plugin="RcppArmadillo", 
# 	body=src)

# computeNetworkSummary_cxx <- cxxfunction(
# 	signature(seq_m="integer", D="integer"),
# 	plugin="RcppArmadillo", 
# 	body=
# 		'
# 		arma::mat seq_m2 = Rcpp::as<arma::mat>(seq_m);
# 		arma::mat DD = Rcpp::as<arma::mat>(D);
# 		int nn = DD.n_rows;
# 		arma::mat D00 = arma::zeros(nn,nn);
# 		arma::mat degreee = arma::zeros(nn,nn);
# 		arma::mat common_frds_11 = arma::zeros(nn,nn);

# 		for ( int i=0 ; i<nn*(nn-1)/2; i++ ){
# 			int index1 = seq_m2(i,0) -1 ;
# 			int index2 = seq_m2(i,1) -1;
# 			if (DD(index1,index2)==1) {
# 				D00(index1,index2) = 1;
# 				D00(index2,index1) = 1;
# 			}

# 			degreee(index1,index2) = sum(D00.col(index1)) ; 
# 			degreee(index2,index1) = sum(D00.col(index2))  ;
# 			common_frds_11(index1,index2) = arma::as_scalar(D00.col(index1).t() * D00.col(index2)) ; 
# 		}
# 		arma::mat out1 = arma::zeros(nn*(nn-1)/2, 3);
# 		arma::mat out2 = arma::zeros(nn*(nn-1)/2, 3);


# 		int k = 0;
# 		for ( int j=0 ; j < nn ; j++){
# 			for ( int i=j+1  ; i< nn ; i++){
# 				out1(k,0) = arma::as_scalar( degreee(i,j) );
# 				out1(k,1) = arma::as_scalar( degreee(i,j)*degreee(i,j) );
# 				out1(k,2) = arma::as_scalar( common_frds_11(i,j) );

# 				out2(k,0) = arma::as_scalar( degreee(j,i) );
# 				out2(k,1) = arma::as_scalar( degreee(j,i)*degreee(j,i) );
# 				out2(k,2) = arma::as_scalar( common_frds_11(i,j) );

# 				k++;
# 			}
# 		}

# 		return Rcpp::List::create(Rcpp::Named("self")=out1, Rcpp::Named("friends")=out2);
# 		'
# 	)

# q1=computeNetworkSummary_r(seq_m[[1]],D[[1]])
# q2=computeNetworkSummary_cxx(seq_m[[1]],D[[1]])
# all.equal(q1$self,q2$self,check.attributes=F)
# all.equal(q1$friends,q2$friends,check.attributes=F)


# benchmark(
# b={
# q2 = computeNetworkSummary_cxx(seq_m[[1]],D[[1]])
#   }
# )


# all.equal(q1$self,q2$self,check.attributes=F)
# all.equal(q1$friends,q2$friends,check.attributes=F)


################################

# library(Matrix)
# load("model_data.rData")
# library(Matrix)
# library(rbenchmark)
# library(compiler)

# D = (data[[1]]$W!=0) + 0 

# n=nrow(D)
# index_table = data[[1]]$group_index
# seq_m = index_table[sample(nn,nn),]

# benchmark(
#   a={
# 		q1= f(seq_m=seq_m, D=D)
#   }
#   ,b={
# 		q2= g(seq_m=seq_m, D=D)
#   }
# )

# all(q1[[1]]==q2[[1]])
# all(q1[[2]]==q2[[2]])


#' computeNetworkSummary
#' @name computeNetworkSummary
#' @aliases computeNetworkSummary
#' @title computeNetworkSummary
#' @param seq_m meeting sequence 
#' @param D adjacency matrix 
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
computeNetworkSummary=function(seq_m, D){
	# out = list()
	# for (i in 1:length(D)){
	# 	temp = computeNetworkSummary_r(g(seq_m=seq_m[[i]], D=D[[i]]))
	# 	out$self = rbind(out$self, temp$self)
	# 	out$friends = rbind(out$friends, temp$friends)
	# }

	# out
  out = mapply(function(x,y) computeNetworkSummary_cxx(seq_m=x, D=y), x= seq_m, y=D  )

  self = do.call(rbind, out[1,] )
  friends = do.call(rbind, out[1,] )

  colnames(self) = c("degree","degree_seq","common_frds")
  colnames(friends) = c("degree","degree_seq","common_frds")

  list(self=self,friends=friends)
}

# q1 = computeNetworkSummary(seq_m,D)



#############################################




# SNF_single_mcmc = function(m, data,  last_estimation, update_tau=TRUE,tau=0.0005){

#   # if (network_id==1){
# 		# D = lapply(data, function(z) z$W!=0 )
# 		# y = do.call(c, lapply(data, "[[", "response1"))
# 		# y_not = !y
#   # } else{
# 		# D = lapply(data, function(z) z$W2!=0 )
# 		# y = do.call(c, lapply(data, "[[", "response2"))
# 		# y_not = !y
#   # }

#   D = lapply(data, function(z) z$D_list[[1]])
#   y = do.call(c, lapply(data,"[[","response_self") )
#   y_not = !y

# 	n= sapply(data,"[[","n")
# 	n2 = sapply(data, function(z) length(z$response_self))

# 	x1 = do.call(rbind, lapply(data, "[[" , "self_data_matrix") )
# 	x2 = do.call(rbind, lapply(data, "[[" , "friends_data_matrix") )
	
# 	number_of_network_variable = 3

# 	postition = mapply(seq, c(0,head(cumsum(n2),-1)) +1,cumsum(n2))
# 	ystar1 = rep(0,length(y))
# 	ystar2 = rep(0,length(y))
# 	seq_m = lapply(data,"[[","group_index")

# 	delta_matrix = matrix(0, nrow=m+1, ncol= ncol(x1) + number_of_network_variable )
# 	update_rate = 0

# 	## initialization 

# 	if (!missing(last_estimation) && !is.null(last_estimation) ){
# 		cat("Using last_estimation \n")
# 		ystar1 = last_estimation$ystar1
# 		ystar2 = last_estimation$ystar2
# 		seq_m = last_estimation$seq_m
# 		delta_matrix[1,] = as.vector(tail(last_estimation$delta,1))

#     index = last_estimation$index+1
#     # name = last_estimation$name
#     ID = last_estimation$ID

# 		if (update_tau){
#     	tau=updateTau(last_estimation$tau, last_estimation$update_rate, lower_bound=0.2, upper_bound=0.4,optim_rate=0.3,min_rate=0.00001)
# 		} else{
# 			tau=last_estimation$tau
# 		}
# 	} else {
#     index = 1
#     ID = genUniqueID()
#     cat("Start new instance with ID ", ID, "\n")
#   }

# 	network_summary = computeNetworkSummary(seq_m=seq_m, D=D)

# 	xx1 = cbind(x1,network_summary$self)
# 	xx2 = cbind(x2,network_summary$friends)

# 	xb1 = xx1 %*% delta_matrix[1,]
# 	xb2 = xx2 %*% delta_matrix[1,]

# 	colnames(delta_matrix) = colnames(xx1)
#   name = colnames(xx1) 

# 	delta_x_index = 1:ncol(x1)
# 	delta_network_index = 1:number_of_network_variable + ncol(x1)

# 	tic()
# 	## start the gibbs
# 	for (i in 1:m){
# 		## base on the seq, compute the network summary
# 		## draw ystar

# 		ystar1 = drawYstar(y=y , ystar_other=ystar2, mean=xb1, y_not= y_not)
# 		ystar2 = drawYstar(y=y, ystar_other=ystar1, mean=xb2, y_not= y_not)

# 		## draw delta
# 		lm_fit = myFastLm(X= rbind(xx1,xx2), y = c(ystar1,ystar2))
# 		delta_matrix[i+1, ] =   mvrnorm(n=1, mu=lm_fit$coef, lm_fit$cov/lm_fit$s^2)

#     R1 = x1 %*% delta_matrix[i+1, delta_x_index]
#     R2 = x2 %*% delta_matrix[i+1, delta_x_index]

#     xb1 = R1 + network_summary$self %*% delta_matrix[i+1,delta_network_index]
#     xb2 = R2 + network_summary$friends %*% delta_matrix[i+1,delta_network_index]


# 		## update sequence 
# 		seq_m_new = DrawSeqSample(seq_m,p=tau)

# 		# sapply(1:5, function(i) sum(seq_m_new[[i]]!=seq_m[[i]]) )

# 		network_summary_new = computeNetworkSummary(seq_m=seq_m_new, D=D)

#     xb1_new = R1 + network_summary_new$self %*% delta_matrix[i+1,delta_network_index]
#     xb2_new = R2 + network_summary_new$friends %*% delta_matrix[i+1,delta_network_index]


# 		p1 = splitBy(dnorm(ystar1 - xb1, log=TRUE),by=n2)
# 		p2 = splitBy(dnorm(ystar2 - xb2, log=TRUE),by=n2)
# 		p1_new = splitBy(dnorm(ystar1 - xb1_new, log=TRUE),by=n2)
# 		p2_new = splitBy(dnorm(ystar2 - xb2_new, log=TRUE),by=n2)

# 		p1 = sapply(p1, sum)
# 		p2 = sapply(p2, sum)
# 		p1_new = sapply(p1_new, sum)
# 		p2_new = sapply(p2_new, sum)

# 		alpha = exp( p1_new+ p2_new - p1- p2  )
# 		update_index = alpha  > runif(5)
# 		seq_m[update_index] = seq_m_new[update_index]
# 		update_rate = update_rate + update_index


# 		update_position = unlist(postition[update_index])
#     network_summary$self[update_position,] = network_summary_new$self[update_position,]
#     network_summary$friends[update_position,] = network_summary_new$friends[update_position,]
# 	  xb1[update_position] = xb1_new[update_position]
# 	  xb2[update_position] = xb2_new[update_position]

# 		xx1[update_position,delta_network_index] = network_summary$self[update_position,]
# 		xx2[update_position,delta_network_index] = network_summary$friends[update_position,]

# 		# test
# 		# xx1_q = cbind(x1,network_summary$self)
# 		# xx2_q = cbind(x2,network_summary$friends)

# 		# network_summary_q = computeNetworkSummary(seq_m=seq_m, D=D)
# 		# xx1_q = cbind(x1,network_summary_q$self)
# 		# xx2_q = cbind(x2,network_summary_q$friends)
# 		# xb1_q = xx1_q %*% delta_matrix[i+1,]
# 		# xb2_q = xx2_q %*% delta_matrix[i+1,]

# 		# identical(xx1,xx1_q)
# 		# identical(xx2,xx2_q)

# 		# identical(xb1,xb1_q)
# 		# identical(xb2,xb2_q)


# 	}
# 	toc()


# 	update_rate = update_rate/m
# 	cat("Update rate : \n")
# 	print(update_rate)

# 	out = list(delta=tail(delta_matrix,-1) , seq_m=seq_m,ystar1=ystar1,ystar2=ystar2, tau=tau, update_rate=update_rate, index=index,ID=ID, name=name)

# 	class(out) = "SNF_single"
# 	out
# }

# merge.SNF_single = function(x,y,...){
#   out = y 
#   out$delta_matrix = rbind(x$delta_matrix, y$delta_matrix)
#   out 
# }

# getParameterMatrix.SNF_single = function(x ){
#   x$delta
# }






# ## update by network
# ## ystar1_demean
# updateSequence = function(ystar1_demean, ystar2_demean, seq_m, tau, delta_network, D){

# 	network_summary = computeNetworkSummary(g(seq_m=seq_m, D=D))

# 	seq_m_new = DrawSeqSample(seq_m,p=tau)

# 	network_summary_new = computeNetworkSummary(g(seq_m=seq_m, D=D))

# 	lik_old = sum(dnorm(ystar1_demean - network_summary$self %*% delta_network, log=TRUE)) + sum(dnorm(ystar2_demean - network_summary$friends %*% delta_network, log=TRUE))
# 	lik_new = sum(dnorm(ystar1_demean - network_summary_new$self %*% delta_network, log=TRUE)) + sum(dnorm(ystar2_demean - network_summary_new$friends %*% delta_network, log=TRUE))
# 	alpha = exp(lik_new - lik_old )
# 	if (alpha>runif(1)){
# 		return(list(seq_m=seq_m_new, update=TRUE))
# 	} else{
# 		return(list(seq_m=seq_m, update=FALSE))
# 	}
# })





# ######parallel 
# library(Matrix)
# load("model_data.rData")
# library(Matrix)
# library(rbenchmark)
# library(compiler)
# library(parallel)


# D = lapply(data, function(z) z$W!=0 )


# n= sapply(data,"[[","n")

# x1 = do.call(rbind, lapply(data, "[[" , "self_data_matrix") )
# x2 = do.call(rbind, lapply(data, "[[" , "friends_data_matrix") )

# y = do.call(c, lapply(data, "[[", "response1"))
# y_not = !y

# n2 = sapply(data, function(z) length(z$response1))


# parameter=list()
# parameter$delta = rep(0, ncol(x1)+3)
# parameter$ystar1 = rep(0,length(y))
# parameter$ystar2 = rep(0,length(y))
# parameter$seq_m = lapply(data,"[[","group_index")
# parameter$tau = parameter$tau
# parameter$m = 100



# ystar1 = parameter$ystar1
# ystar2 = parameter$ystar2
# seq_m = parameter$seq_m
# tau = parameter$tau 
# m = parameter$m 


# delta_matrix = matrix(0,nrow=m+1,ncol=length(parameter$delta))
# delta_matrix[1,] = parameter$delta


# ## initialization 


# network_summary = computeNetworkSummary(seq_m=seq_m, D=D)

# xx1 = cbind(x1,network_summary$self)
# xx2 = cbind(x2,network_summary$friends)
# xb1 = xx1 %*% delta_matrix[1,]
# xb2 = xx2 %*% delta_matrix[1,]


# cl=makeCluster(6)
# exportAllFunction(cl)
# clusterExport(cl,c("D","src"))
# clusterEvalQ(cl,{library(inline);library(RcppArmadillo)})
# clusterEvalQ(cl,{g <- cxxfunction(signature(seq_m="integer", D="integer"),
# 	plugin="RcppArmadillo", 
# 	body=src)
# })


# tic()
# ## start the gibbs
# for (i in 1:m){
# 	## base on the seq, compute the network summary
# 	## draw ystar
# 	network_summary = computeNetworkSummary(seq_m=seq_m, D=D)

# 	xx1 = cbind(x1,network_summary$self)
# 	xx2 = cbind(x2,network_summary$friends)
# 	xb1 = xx1 %*% delta_matrix[1,]
# 	xb2 = xx2 %*% delta_matrix[1,]

# 	ystar1 = drawYstar(y=y , ystar_other=ystar2, mean=xb1, y_not= y_not)
# 	ystar2 = drawYstar(y=y, ystar_other=ystar1, mean=xb2, y_not= y_not)

# 	## draw delta
# 	lm_fit = myFastLm(XX= rbind(xx1,xx2), yy = c(ystar1,ystar2))
# 	delta_matrix[i+1, ] =   mvrnorm(n=1, mu=lm_fit$coef, lm_fit$cov/lm_fit$s^2)

# 	xb1 = xx1 %*% delta_matrix[i+1,]
# 	xb2 = xx2 %*% delta_matrix[i+1,]


# 	ystar1_demean = ystar1 - x1 %*% head(delta_matrix[i+1,],ncol(x1))
# 	ystar2_demean = ystar2 - x2 %*% head(delta_matrix[i+1,],ncol(x1))

# 	ystar1_demean_list = splitBy(ystar1_demean,n2)
# 	ystar2_demean_list = splitBy(ystar2_demean,n2)
# 	out = parLapply(cl, 1:length(D), 
# 		function(x,ystar1_demean_list,ystar2_demean_list, seq_m, tau, delta ) {
# 				updateSequence(
# 					ystar1_demean=ystar1_demean_list[[x]],
# 					ystar2_demean=ystar2_demean_list[[x]],
# 					seq_m= seq_m[[x]],
# 					tau = tau ,
# 					delta_network=delta,
# 					D=D[[x]]
# 				)
# 		},
# 		delta=tail(delta_matrix[i+1,],-ncol(x1)),
# 		ystar1_demean_list = ystar1_demean_list,
# 		ystar2_demean_list = ystar2_demean_list,
# 		tau=tau,
# 		seq_m=seq_m
# 	)
# 	seq_m = lapply(out,"[[","seq_m" )
# 	update = update + out$update
# }
# toc()



##############################################################################################################
## Given a seq, beta, D, compute the likelihood
## Given a seq m  n*(n-1)/2 x 2 matrix
## D: n by n network matrix
## U_xb : an n by n utility matrix, i,j element is the utility of i to make friends with j. (Xbeta)
## delta1 delta2 


# DrawSeqSample = function(x , p =0.01){
#   if (is.list(x)){
#     out = vector("list",length(x))
#     if (length(p)==1)
#       p = rep(p,length(x))
#     for (i in 1:length(x)){
#       out[[i]] = DrawSeqSample(x[[i]], p[[i]])
#     }
#     return(out)
#   } else if (is.vector(x)){
#     n = length(x)
#     pn = pmax(2,ceiling(n * p)  )
#     to_change = sample(n,pn)
#     reorder = sample(to_change)

#     x[to_change] = x[reorder]
#     return(x)
#   } else if (is.matrix(x)){
#     n = nrow(x)
#     pn = pmax(2,ceiling(n * p)  )
#     to_change = sample(n,pn)
#     reorder = sample(to_change)
#     x[to_change,] = x[reorder,]
#     return(x) 
#   }
# })

# # repeat

# computeNetworkSummary = function(seq_m,D){
#   n = nrow(D)
#   D0 = matrix(0,ncol=n,nrow=n)
#   degree = matrix(0,ncol=n,nrow=n)
#   common_frds_1 = matrix(0,ncol=n,nrow=n)
#   common_frds_2 = matrix(0,ncol=n,nrow=n)

#   nn = n*(n-1)/2
#   for ( i in 1:nn){
#     index = seq_m[i,]
#     index1 = index[1]
#     index2 = index[2]
#     if (D[index1,index2]==1) {
#       D0[index1,index2] = D0[index2,index1]= 1
#     }
#     d1 = D0[index1,]
#     d2 = D0[index2,]
#     degree[index1,index2] = sum(d1)
#     degree[index2,index1] = sum(d2)
#     common_frds_1[index1,index2] = sum(d1*d2)
#   }

#   lower_tri = lower.tri(degree)
#   degree1 = degree[ lower_tri ] 
#   degree2 = t(degree)[ lower_tri ] 
#   common_frds_1 = common_frds_1[ lower_tri ] 

#   list(self=cbind(degree1,degree1^2,common_frds_1), friends=cbind(degree2,degree2^2,common_frds_1))
# })


# src=
# '
# arma::mat seq_m2 = Rcpp::as<arma::mat>(seq_m);
# arma::mat DD = Rcpp::as<arma::mat>(D);
# int nn = DD.n_rows;
# arma::mat D00 = arma::zeros(nn,nn);
# arma::mat degreee = arma::zeros(nn,nn);
# arma::mat common_frds_11 = arma::zeros(nn,nn);

# for ( int i=0 ; i<nn*(nn-1)/2; i++ ){
#   int index1 = seq_m2(i,0) -1 ;
#   int index2 = seq_m2(i,1) -1;
#   if (DD(index1,index2)==1) {
#     D00(index1,index2) = 1;
#     D00(index2,index1) = 1;
#   }

#   degreee(index1,index2) = sum(D00.col(index1)) ; 
#   degreee(index2,index1) = sum(D00.col(index2))  ;
#   common_frds_11(index1,index2) = arma::as_scalar(D00.col(index1).t() * D00.col(index2)) ; 
# }

# return Rcpp::List::create(Rcpp::Named("degree")=degreee, Rcpp::Named("common_frds_1")=common_frds_11);
# '

# g <- cxxfunction(signature(seq_m="integer", D="integer"),
#   plugin="RcppArmadillo", 
#   body=src)

# computeNetworkSummary_cxx <- cxxfunction(
#   signature(seq_m="integer", D="integer"),
#   plugin="RcppArmadillo", 
#   body=
#     '
#     arma::mat seq_m2 = Rcpp::as<arma::mat>(seq_m);
#     arma::mat DD = Rcpp::as<arma::mat>(D);
#     int nn = DD.n_rows;
#     arma::mat D00 = arma::zeros(nn,nn);
#     arma::mat degreee = arma::zeros(nn,nn);
#     arma::mat common_frds_11 = arma::zeros(nn,nn);

#     for ( int i=0 ; i<nn*(nn-1)/2; i++ ){
#       int index1 = seq_m2(i,0) -1 ;
#       int index2 = seq_m2(i,1) -1;
#       if (DD(index1,index2)==1) {
#         D00(index1,index2) = 1;
#         D00(index2,index1) = 1;
#       }

#       degreee(index1,index2) = sum(D00.col(index1)) ; 
#       degreee(index2,index1) = sum(D00.col(index2))  ;
#       common_frds_11(index1,index2) = arma::as_scalar(D00.col(index1).t() * D00.col(index2)) ; 
#     }
#     arma::mat out1 = arma::zeros(nn*(nn-1)/2, 3);
#     arma::mat out2 = arma::zeros(nn*(nn-1)/2, 3);


#     int k = 0;
#     for ( int j=0 ; j < nn ; j++){
#       for ( int i=j+1  ; i< nn ; i++){
#         out1(k,0) = arma::as_scalar( degreee(i,j) );
#         out1(k,1) = arma::as_scalar( degreee(i,j)*degreee(i,j) );
#         out1(k,2) = arma::as_scalar( common_frds_11(i,j) );

#         out2(k,0) = arma::as_scalar( degreee(j,i) );
#         out2(k,1) = arma::as_scalar( degreee(j,i)*degreee(j,i) );
#         out2(k,2) = arma::as_scalar( common_frds_11(i,j) );

#         k++;
#       }
#     }

#     return Rcpp::List::create(Rcpp::Named("self")=out1, Rcpp::Named("friends")=out2);
#     '
#   )

# q1=computeNetworkSummary(seq_m[[1]],D[[1]])
# q2=computeNetworkSummary_cxx(seq_m[[1]],D[[1]])
# all.equal(q1$self,q2$self,check.attributes=F)
# all.equal(q1$friends,q2$friends,check.attributes=F)


# benchmark(
# b={
# q2 = computeNetworkSummary_cxx(seq_m[[1]],D[[1]])
#   }
# )


# all.equal(q1$self,q2$self,check.attributes=F)
# all.equal(q1$friends,q2$friends,check.attributes=F)


################################

# library(Matrix)
# load("model_data.rData")
# library(Matrix)
# library(rbenchmark)
# library(compiler)

# D = (data[[1]]$W!=0) + 0 

# n=nrow(D)
# index_table = data[[1]]$group_index
# seq_m = index_table[sample(nn,nn),]

# benchmark(
#   a={
#     q1= f(seq_m=seq_m, D=D)
#   }
#   ,b={
#     q2= g(seq_m=seq_m, D=D)
#   }
# )

# all(q1[[1]]==q2[[1]])
# all(q1[[2]]==q2[[2]])


# computeNetworkSummary=function(seq_m=seq_m_new, D=D){
#   # out = list()
#   # for (i in 1:length(D)){
#   #   temp = computeNetworkSummary(g(seq_m=seq_m[[i]], D=D[[i]]))
#   #   out$self = rbind(out$self, temp$self)
#   #   out$friends = rbind(out$friends, temp$friends)
#   # }

#   # out
#   out = mapply(function(x,y) computeNetworkSummary_cxx(seq_m=x, D=y), x= seq_m, y=D  )

#   self = do.call(rbind, out[1,] )
#   friends = do.call(rbind, out[1,] )

#   colnames(self) = c("degree","degree_seq","common_frds")
#   colnames(friends) = c("degree","degree_seq","common_frds")

#   list(self=self,friends=friends)
# })

# # q1 = computeNetworkSummary(seq_m,D)


# computeConditionalVariance = function(Sigma){
#   k = nrow(Sigma)
#   ols_coef = vector("list", k)
#   sd_new = vector("list",k)
#   for (i in 1:k){
#     ols_coef[[i]] = Sigma[i,-i] %*% solve(Sigma[-i,-i])
#     sd_new[[i]] = sqrt ( Sigma[i,i] - ols_coef[[i]] %*% Sigma[-i,i] )
#   }
#   list(sd=sd_new, ols_coef=ols_coef)
# })
# #############################################


# update_seq_multi = function(seq_m, D_list, xb1, xb2, x1_network, x2_network, delta_network_index, ystar1,ystar2, Sigma,n2,update_rate, tau){
#     seq_m_new = DrawSeqSample(seq_m,p=tau)

#     network_summary_new = lapply(D_list,computeNetworkSummary, seq_m=seq_m_new )

#     x1_network_new = do.call(cbind, lapply(network_summary_new, "[[", "self")) 
#     x2_network_new = do.call(cbind, lapply(network_summary_new, "[[", "friends")) 

#     xb1_new = R1 + x1_network_new %*% delta[delta_network_index,]
#     xb2_new = R2 + x2_network_new %*% delta[delta_network_index,]



#     p1 = splitBy( dmvnorm(ystar1 - xb1, sigma=Sigma, log=TRUE),by=n2)
#     p2 = splitBy( dmvnorm(ystar2 - xb2, sigma=Sigma, log=TRUE),by=n2)
#     p1_new = splitBy( dmvnorm(ystar1 - xb1_new, sigma=Sigma, log=TRUE),by=n2)
#     p2_new = splitBy( dmvnorm(ystar2 - xb2_new, sigma=Sigma, log=TRUE),by=n2)

#     p1 = sapply(p1, sum)
#     p2 = sapply(p2, sum)
#     p1_new = sapply(p1_new, sum)
#     p2_new = sapply(p2_new, sum)

#     alpha = exp( p1_new+ p2_new - p1- p2  )
#     update_index = alpha  > runif(5)
#     seq_m[update_index] = seq_m_new[update_index]
#     update_rate = update_rate + update_index


#     update_position = unlist(position [update_index])

#     x1_network[update_position,] = x1_network_new[update_position,]
#     x2_network[update_position,] = x2_network_new[update_position,]

#     xb1[update_position] = xb1_new[update_position]
#     xb2[update_position] = xb2_new[update_position]

#     list(seq_m, xb1, xb2, x1_network, x2_network,update_rate)
# })


# drawYstar_multi_SNF = function(y, y_not=!y, ystar, ystar_other, xb,  Sigma){
#   # update ystar given ystar_other
#   ystar_other_positive = ystar_other>=0  
#   number_of_network = ncol(y)
#   n = nrow(y)
#   for (i in 1:number_of_network){
#     ols_coef = Sigma[i,-i] %*% solve(Sigma[-i,-i])
#     sd_new = sqrt ( Sigma[i,i] - ols_coef %*% Sigma[-i,i] )
#     mean_new = xb1[,i] + ols_coef %*% (ystar[,-i] - xb[,-i])

#     index_case1 = y[,i]
#     index_case2 = as.logical(y_not[,i] * ystar_other_positive[,i])
#     index_case3 = as.logical(y_not[,i] * !ystar_other_positive[,i])

#     n1=sum(index_case1)
#     n2=sum(index_case2)
#     n3=sum(index_case3)

#     stopifnot(n==n1+n2+n3)

#     if (n1>0)
#       ystar[index_case1,i] = rtruncnorm(1,a=0,b=Inf, mean=mean_new[index_case1], sd=sd_new)  
#     if (n2>0)
#       ystar[index_case2,i] =rtruncnorm(1,a=-Inf,b=0,mean=mean_new[index_case2], sd=sd_new) 
#     if (n3>0)
#       ystar[index_case3,i] = mean_new[index_case3] +rnorm(n3, sd=sd_new)
#   }

#   ystar
# })




# SNF_single_mcmc = function(m, data, network_id, last_estimation, update_tau=TRUE,tau=0.005){

#   if (network_id==1){
#     D = lapply(data, function(z) z$W!=0 )
#     y = do.call(c, lapply(data, "[[", "response1"))
#     y_not = !y
#   } else{
#     D = lapply(data, function(z) z$W2!=0 )
#     y = do.call(c, lapply(data, "[[", "response2"))
#     y_not = !y
#   }
#   n= sapply(data,"[[","n")
#   n2 = sapply(data, function(z) length(z$response1))

#   x1 = do.call(rbind, lapply(data, "[[" , "self_data_matrix") )
#   x2 = do.call(rbind, lapply(data, "[[" , "friends_data_matrix") )
  
#   number_of_network_variable = 3

#   position  = mapply(seq, c(0,head(cumsum(n2),-1)) +1,cumsum(n2))
#   ystar1 = rep(0,length(y))
#   ystar2 = rep(0,length(y))
#   seq_m = lapply(data,"[[","group_index")

#   delta_matrix = matrix(0, nrow=m+1, ncol= ncol(x1) + number_of_network_variable )
#   update_rate = 0


#   if (!missing(last_estimation) && !is.null(last_estimation) ){
#     cat("Using last_estimation \n")
#     ystar1 = last_estimation$ystar1
#     ystar2 = last_estimation$ystar2
#     seq_m = last_estimation$seq_m
#     delta_matrix[1,] = as.vector(tail(last_estimation$delta,1))
    
#     if (update_tau){
#       tau = last_estimation$tau 
#       for ( j in 1:length(tau)){
#         if (any(last_estimation$update_rate[[j]] >0.5 | any(last_estimation$update_rate[[j]]< 0.2) )){
#           cat("update tau-", j , "\n")
#           tau[[j]] = tau[[j]] * last_estimation$update_rate[[j]] / 0.4
#           tau[[j]] = ifelse(tau[[j]]==0, 0.0001, tau[[j]])
#         }
#       }
#     }
#   }

#   ## initialization 
#   network_summary = computeNetworkSummary(seq_m=seq_m, D=D)

#   xx1 = cbind(x1,network_summary$self)
#   xx2 = cbind(x2,network_summary$friends)
#   xb1 = xx1 %*% delta_matrix[1,]
#   xb2 = xx2 %*% delta_matrix[1,]

#   colnames(delta_matrix) = colnames(xx1)

#   delta_x_index = 1:ncol(x1)
#   delta_network_index = 1:number_of_network_variable + ncol(x1)

#   tic()
#   ## start the gibbs
#   for (i in 1:m){
#     ## base on the seq, compute the network summary
#     ## draw ystar

#     ystar1 = drawYstar(y=y , ystar_other=ystar2, mean=xb1, y_not= y_not)
#     ystar2 = drawYstar(y=y, ystar_other=ystar1, mean=xb2, y_not= y_not)

#     ## draw delta
#     lm_fit = my.fastLm(XX= rbind(xx1,xx2), yy = c(ystar1,ystar2))
#     delta_matrix[i+1, ] =   mvrnorm(n=1, mu=lm_fit$coef, lm_fit$cov/lm_fit$s^2)

#     R1 = x1 %*% delta_matrix[i+1, delta_x_index]
#     R2 = x2 %*% delta_matrix[i+1, delta_x_index]

#     xb1 = R1 + network_summary$self %*% delta_matrix[i+1,delta_network_index]
#     xb2 = R2 + network_summary$friends %*% delta_matrix[i+1,delta_network_index]


#     ## update sequence 
#     seq_m_new = DrawSeqSample(seq_m,p=tau)
#     network_summary_new = computeNetworkSummary(seq_m=seq_m_new, D=D)

#     xb1_new = R1 + network_summary_new$self %*% delta_matrix[i+1,delta_network_index]
#     xb2_new = R2 + network_summary_new$friends %*% delta_matrix[i+1,delta_network_index]


#     p1 = splitBy(dnorm(ystar1 - xb1, log=TRUE),by=n2)
#     p2 = splitBy(dnorm(ystar2 - xb2, log=TRUE),by=n2)
#     p1_new = splitBy(dnorm(ystar1 - xb1_new, log=TRUE),by=n2)
#     p2_new = splitBy(dnorm(ystar2 - xb2_new, log=TRUE),by=n2)

#     p1 = sapply(p1, sum)
#     p2 = sapply(p2, sum)
#     p1_new = sapply(p1_new, sum)
#     p2_new = sapply(p2_new, sum)

#     alpha = exp( p1_new+ p2_new - p1- p2  )
#     update_index = alpha  > runif(5)
#     seq_m[update_index] = seq_m_new[update_index]
#     update_rate = update_rate + update_index


#     update_position = unlist(position [update_index])
#     network_summary$self[update_position,] = network_summary_new$self[update_position,]
#     network_summary$friends[update_position,] = network_summary_new$friends[update_position,]
#     xb1[update_position] = xb1_new[update_position]
#     xb2[update_position] = xb2_new[update_position]

#     xx1[update_position,delta_network_index] = network_summary$self[update_position,]
#     xx2[update_position,delta_network_index] = network_summary$friends[update_position,]


#   }
#   toc()


#   update_rate = update_rate/m
#   cat("Update rate : \n")
#   print(update_rate)

#   out = list(delta=tail(delta_matrix,-1) , seq_m=seq_m,ystar1=ystar1,ystar2=ystar2, tau=tau, update_rate=update_rate)

#   class(out) = "SNF_single"
#   out
# })

# merge.SNF_single = function(x,y,...){
#   if (length(list(...))>0){
#     list_args = list(...)
#     out = list_args[[1]]
#     for (i in 2:length(list_args)){
#       out = merge(out, list_args[[i]])
#     }
#     return(out)
#   }

#   out =y 
#   out$delta = rbind(x$delta, y$delta)
#   out
# })

# plotmcmc.SNF_single = function(x,remove=0.2,...){
#   plotmcmc(x$delta, remove=nrow(x$delta)*remove )
# })



# getParameterMatrix.SNF_single = function(x, tail){
#   out = cbind(x$delta)

#   if (tail!=0) {
#     out = tail(out,tail)
#   }
#   out
# })




















#' @rdname SNF
#' @export
SNF.dynamic.mcmc = function(m, data, last_estimation, update_tau=TRUE,tau=0.005){

  G = length(data)
  number_of_network = length(data[[1]]$D_list)

  D_list = vector("list", number_of_network)
  for (i in 1:number_of_network){
    D_list[[i]] = lapply(data, function(z) z$D_list[[i]])
  }

  # D_list = lapply(data, "[[", "D_list")

  n= sapply(data,"[[","n")
  n2 = sapply(data, function(z) NROW(z$response_self))

  nn = sum(n2) * 2

  y = do.call(rbind, lapply(data, "[[", "response_self") )

  y_not = !y

  x1 = do.call(rbind, lapply(data, "[[" , "self_data_matrix") )
  x2 = do.call(rbind, lapply(data, "[[" , "friends_data_matrix") )
  
  number_of_network_variable = 3
  number_of_variable = ncol(x1) + number_of_network_variable*number_of_network
  ## store all the network matrix of different group into one vector. position is the location of them.
  position  = mapply(seq, c(0,head(cumsum(n2),-1)) +1,cumsum(n2))

  ystar1 = array(0, dim=dim(y))
  ystar2 = array(0, dim=dim(y))

  seq_m = lapply(data,"[[","group_index")

  delta_matrix = rep(list(matrix(0, nrow=m, ncol= number_of_variable )), number_of_network)
  delta= matrix(0, nrow=number_of_variable , ncol=number_of_network) 
  for (i in 1:number_of_network){
    delta[,i] = delta_matrix[[i]][1,]
  }

  update_rate = 0
  Sigma = diag(number_of_network)

  if (!missing(last_estimation) && !is.null(last_estimation) ){
    cat("Using last_estimation \n")
    ystar1 = last_estimation$ystar1
    ystar2 = last_estimation$ystar2
    seq_m = last_estimation$seq_m
    delta = last_estimation$delta
    Sigma = last_estimation$Sigma

    if (update_tau){
      tau = last_estimation$tau 
      for ( j in 1:length(tau)){
        if (any(last_estimation$update_rate[[j]] >0.5 | any(last_estimation$update_rate[[j]]< 0.2) )){
          cat("update tau-", j , "\n")
          tau[[j]] = tau[[j]] * last_estimation$update_rate[[j]] / 0.4
          tau[[j]] = ifelse(tau[[j]]==0, 0.0001, tau[[j]])
        }
      }
    }
  }

  ## initialization 

  network_summary = lapply(D_list,computeNetworkSummary, seq_m=seq_m )

  ## repeat that for serial D.

  # network_summary reduce to 1 variable ( # of common frds, then include all the network , or we can have all, it would be 3 x number of network )

  x1_network = do.call(cbind, lapply(network_summary, "[[", "self")) 
  x2_network = do.call(cbind, lapply(network_summary, "[[", "friends"))               
  delta_x_index = 1:ncol(x1)
  delta_network_index = ncol(x1)+ seq(ncol(x1_network))


  R1 = x1 %*% delta[delta_x_index,]
  R2 = x2 %*% delta[delta_x_index,]

  xb1 = R1 + x1_network %*% delta[delta_network_index,]
  xb2 = R2 + x2_network %*% delta[delta_network_index,]

  network_name = data[[1]]$network_name

  rownames(delta) = c( colnames(x1) , colnames(x1_network) )

  colname_network = colnames(x1_network)[1:3]

  colname_network = unlist( lapply(network_name, function(z) z %+% "_" %+% colname_network) )


  for (i in 1:number_of_network){
    colnames(delta_matrix[[i]]) = c(network_name[[i]] %+% "_" %+% colnames(x1) , network_name[[i]] %+% "_" %+%  colname_network )  
  }
  X= rbind(x1,x2)

  if (number_of_network>1){
    number_col_Sigma_matrix = number_of_network*(number_of_network-1)/2
    sigma_name = genPairwiseIndex(number_of_network)
    sigma_name = sigma_name[,1] %+% sigma_name[,2]
  } else {
    number_col_Sigma_matrix =1
    sigma_name = 11
  }
  Sigma_matrix = matrix(0, nrow=m, ncol = number_col_Sigma_matrix )
  colnames(Sigma_matrix) = "Sigma_" %+%  sigma_name

  tic()
  ## start the gibbs
  for (i in 1:m){

    ##### base on the seq, compute the network summary
    ##### drawing ystar
    for (j in 1:number_of_network){
      ystar1_demean = ystar1 - xb1
      ystar2_demean = ystar2 - xb2

      temp = find_normal_conditional_dist(a= ystar1_demean, i=j, j=-j, Sigma=Sigma)

      ystar1[,j] = drawYstar(y=y[,j] , ystar_other=ystar2[,j], mean=xb1[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )

      temp = find_normal_conditional_dist(a= ystar2_demean, i=j, j=-j, Sigma=Sigma)

      ystar2[,j] = drawYstar(y=y[,j] , ystar_other=ystar1[,j], mean=xb2[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )
    }

    ystar1_demean = ystar1 - xb1
    ystar2_demean = ystar2 - xb2
    ystar_demean = rbind(ystar1_demean,ystar2_demean)
    ystar = rbind(ystar1,ystar2)


    ##### draw delta
    XX = cbind(X, rbind(x1_network,x2_network) )
    YY = rbind(ystar1,ystar2)

    for ( j in 1:number_of_network){
      temp = find_normal_conditional_dist(a=ystar_demean, i=j, j=-j, Sigma=Sigma)
      lm_fit = myFastLm(X=XX, y =YY[,j]-temp$mean)

      delta[,j] = mvrnorm(n=1, mu=lm_fit$coef, lm_fit$cov/lm_fit$s^2 * as.vector(temp$var) )

      delta_matrix[[j]][i,] = delta[,j]
    }

    R1 = x1 %*% delta[delta_x_index,]
    R2 = x2 %*% delta[delta_x_index,]
    xb1 = R1 + x1_network %*% delta[delta_network_index,]
    xb2 = R2 + x2_network %*% delta[delta_network_index,]

    ##### update sequence , only need to do it once. but may need to modify the likelihood 
    seq_m_new = DrawSeqSample(seq_m,p=tau)

    network_summary_new = lapply(D_list,computeNetworkSummary, seq_m=seq_m_new )

    x1_network_new = do.call(cbind, lapply(network_summary_new, "[[", "self")) 
    x2_network_new = do.call(cbind, lapply(network_summary_new, "[[", "friends")) 

    xb1_new = R1 + x1_network_new %*% delta[delta_network_index,]
    xb2_new = R2 + x2_network_new %*% delta[delta_network_index,]



    p1 = splitBy( dmvnorm(ystar1 - xb1, sigma=Sigma, log=TRUE),by=n2)
    p2 = splitBy( dmvnorm(ystar2 - xb2, sigma=Sigma, log=TRUE),by=n2)
    p1_new = splitBy( dmvnorm(ystar1 - xb1_new, sigma=Sigma, log=TRUE),by=n2)
    p2_new = splitBy( dmvnorm(ystar2 - xb2_new, sigma=Sigma, log=TRUE),by=n2)

    p1 = sapply(p1, sum)
    p2 = sapply(p2, sum)
    p1_new = sapply(p1_new, sum)
    p2_new = sapply(p2_new, sum)

    alpha = exp( p1_new+ p2_new - p1- p2  )
    update_index = alpha  > runif(5)

    update_rate = update_rate + update_index

    seq_m[update_index] = seq_m_new[update_index]

    update_position = unlist(position[update_index])

    x1_network[update_position,] = x1_network_new[update_position,]
    x2_network[update_position,] = x2_network_new[update_position,]

    xb1[update_position] = xb1_new[update_position]
    xb2[update_position] = xb2_new[update_position]

    ##### compute Sigma
    ystar1_demean = ystar1 - xb1
    ystar2_demean = ystar2 - xb2
    ystar_demean = rbind(ystar1_demean,ystar2_demean)
    

    if (number_of_network > 1 ){
      Sigma = solve( rwish(nn , solve( crossprod(ystar_demean )) ) )
      normalization = diag(1/sqrt(diag(Sigma)))

      Sigma = normalization %*%  Sigma %*% t(normalization) 
      Sigma_matrix[i,] = Sigma[lower.tri(Sigma)]
    } else {
      Sigma = as.matrix(1)
      Sigma_matrix[i,] = 1
    }

  }
  toc()

  update_rate = update_rate/m
  cat("Update rate : \n")
  print(update_rate)

  out = list(delta_matrix=delta_matrix, delta=delta, seq_m=seq_m, ystar1=ystar1,ystar2=ystar2, tau=tau, update_rate=update_rate, Sigma=Sigma,Sigma_matrix=Sigma_matrix)

  class(out) = "SNF.dynamic.mcmc"
  out
}

#' Get a matrix of parameter 
#' @name getParameterMatrix.SNF.dynamic.mcmc
#' @aliases getParameterMatrix.SNF.dynamic.mcmc
#' @title getParameterMatrix.SNF.dynamic.mcmc
#' @param x SNF.dynamic.mcmc object
#' @param tail iteration to be used. Negative value: Removing the first \code{tail} iterations. Positive value: keep the last \code{tail} iterations. If -1< code{tail}< 1, it represent the percentage of iterations.
#'' @param ... not used
#' @return A matrix
#' @method getParameterMatrix SNF.dynamic.mcmc
#' @export getParameterMatrix SNF.dynamic.mcmc
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
getParameterMatrix.SNF.dynamic.mcmc = function(x, tail, ...){
  if (is.list(x$delta_matrix)){
    out = do.call(cbind, x$delta_matrix)
    out = cbind(out, x$Sigma_matrix )
  } else{
    out = x$delta_matrix
  }
  if (!missing(tail)) {
    out = extractTail(out, tail)
  }
  out
}

#' merge.SNF.dynamic.mcmc
#' @name merge.SNF.dynamic.mcmc
#' @aliases merge.SNF.dynamic.mcmc
#' @title merge.SNF.dynamic.mcmc
#' @param x First object to merge with
#' @param y Second object to merge with 
#' @param ... not used
#' @return A new SNF.dynamic.mcmc object
#' @method merge SNF.dynamic.mcmc
#' @export merge SNF.dynamic.mcmc
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export

merge.SNF.dynamic.mcmc = function(x,y,...){
  out = y 

  for (i in 1:length(y$delta_matrix)){
    out$delta_matrix[[i]] = rbind(x$delta_matrix[[i]], y$delta_matrix[[i]])
  }

  out$Sigma_matrix = rbind(x$Sigma_matrix, y$Sigma_matrix)

  out
}

#' Create a summary table
#' @name summary.SNF.dynamic.mcmc
#' @aliases summary.SNF.dynamic.mcmc
#' @title summary.SNF.dynamic.mcmc
#' @param object SNF.dynamic.mcmc object
#' @param ... tail:  iteration to be used. Negative value: Removing the first \code{tail} iterations. Positive value: keep the last \code{tail} iterations. If -1< code{tail}< 1, it represent the percentage of iterations.
#' @return A summary table
#' @method summary SNF.dynamic.mcmc
#' @export summary SNF.dynamic.mcmc
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export

summary.SNF.dynamic.mcmc = function(object,...){
  computeSummaryTable(object,...)
}


