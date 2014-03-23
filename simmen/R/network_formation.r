loglikelihood_network_formation = function(y, x1, x2, delta, y_not){
  if (missing(y_not))
    y_not = !y
  p = pnorm(x1 %*% delta) * pnorm(x2 %*% delta)
  out = sum(log(p^y*(1-p)^y_not))

  if (!is.finite(out)){
    return(-1e+20)
  }
  return(out)
}

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

lik_grad_single_network_formation = function(x1,x2,y, delta,...){
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


network_formation = function(data){
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

    out[[i]] = maxLik(function(z, ...) loglikelihood_network_formation(delta=z, ...) , start=start ,  x1=self_data_matrix, x2=friends_data_matrix, y=yy, y_not=yy_not , grad= lik_grad_single_network_formation, method="BFGS")

      summary_table[[i]] = generateSignificance(summary(out[[i]])$estimate[,1:2])
      rownames(summary_table[[i]]) = network_name[i] %+% "_" %+%colnames(self_data_matrix)
  }

  summary_table = do.call(rbind,summary_table)
  toc()
  list(out=out, summary_table=summary_table)
}

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

drawYstar_single = function(y, ystar_other, mean, y_not=!y, sd=1){

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


single_network_formation_mcmc = function(data,  m=1000, last_out){
  self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
  friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))

  response = do.call(rbind, lapply(data, function(z) z$response_self))
  response_not = !response

  name = colnames(self_data_matrix)
  k =  ncol(self_data_matrix)
  delta_matrix = matrix(0, nrow=m+1, ncol=k)
  ystar1 = rep(0, length(response))
  ystar2 = rep(0, length(response))
  network_name = data[[1]]$network_name

  if (!missing(last_out)){
    delta_matrix[1,] = tail(last_out$delta_matrix,1)
    ystar1 = last_out$ystar1
    ystar2 = last_out$ystar2
  }

  colnames(delta_matrix) = network_name %+% "_" %+% colnames(self_data_matrix)

  X = rbind(self_data_matrix, friends_data_matrix)
  XX_inv = solve(crossprod(X))

  tic()
  for (i in 1:m){
    if (i %% 1000 == 0 ){
      cat(i ,">\n")
    }

    xb1 = self_data_matrix %*% delta_matrix[i, ]
    xb2 = friends_data_matrix %*% delta_matrix[i, ]

    ystar1 = drawYstar_single(y=response , ystar_other=ystar2, mean=xb1, y_not= response_not)
    ystar2 = drawYstar_single(y=response, ystar_other=ystar1, mean=xb2, y_not= response_not)

    delta_matrix[i+1, ] =   mvrnorm(n=1, mu=XX_inv %*% crossprod(X,c(ystar1,ystar2)), XX_inv)

  }
  toc()
  delta_matrix = tail(delta_matrix,-1)

  out = list(delta_matrix=delta_matrix, ystar1=ystar1, ystar2=ystar2)
  class(out) = "network_formation" 
  out
}


# single_network_formation_mcmc_RE = function(data, network_id, m=1000, last_out){
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

#   if (!missing(last_out)){
#     delta_matrix[1,] = tail(last_out$delta_matrix,1)
#     sigma2e_matrix[1,] = tail(last_out$sigma2e_matrix,1)
#     ystar1 = last_out$ystar1
#     ystar2 = last_out$ystar2
#     e = last_out$e
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
#     ystar1 = drawYstar_single(y=response , ystar_other=ystar2, mean=xb1, y_not= response_not)
#     ystar2 = drawYstar_single(y=response, ystar_other=ystar1, mean=xb2, y_not= response_not)
    
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


# single_network_formation_mcmc_parallel = function(data,cl, network_id, m=1000, last_out){

#   name = colnames(data[[1]]$self_data_matrix)
#   k =  ncol(data[[1]]$self_data_matrix)
#   n2 = sapply(data,function(z) length(z$response1)) 
#   G = length(data)
#   delta_matrix = matrix(0, nrow=m+1, ncol=k)
#   ystar1 = lapply(n2, rep,x=0 )

#   ystar2 = lapply(n2, rep,x=0 )

#   if (!missing(last_out)){
#     delta_matrix[1,] = tail(last_out$delta_matrix,1)
#     ystar1 = tail(last_out$ystar1)
#     ystar2 = tail(last_out$ystar2)
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
#       drawYstar_single(
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
#       drawYstar_single(
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

multi_network_formation_mcmc = function(data, m=1000, last_out){

  self_data_matrix = do.call(rbind, lapply(data , function(z) z$self_data_matrix))
  friends_data_matrix = do.call(rbind, lapply(data , function(z) z$friends_data_matrix))
  
  y = do.call(rbind, lapply(data,"[[","response_self"))
  y_not = !y
  number_of_network = ncol(y)

  if (number_of_network==1){
    return( single_network_formation_mcmc(data,m,last_out) )
  }

  name = colnames(self_data_matrix)
  k =  ncol(self_data_matrix)
  n = NROW(y)*2
  delta_matrix = rep(list(matrix(0, nrow=m, ncol= k )), number_of_network)

  number_col_Sigma_matrix = number_of_network*(number_of_network-1)/2
  Sigma_matrix = matrix(0, nrow=m, ncol = number_col_Sigma_matrix )

  sigma_name = genPairwiseIndex(number_of_network)
  sigma_name = sigma_name[,1] %+% sigma_name[,2]
  
  colnames(Sigma_matrix) = "Sigma_" %+%  sigma_name

  network_name = data[[1]]$network_name

  for (i in 1:number_of_network){
    colnames(delta_matrix[[i]]) = network_name[[i]] %+% "_" %+% name
  }

  ystar1 = matrix(0, nrow=nrow(y), ncol=ncol(y))
  ystar2 = matrix(0, nrow=nrow(y), ncol=ncol(y))
  Sigma = matrix(0.5,number_of_network,number_of_network)
  diag(Sigma) = 1
  delta = matrix(0, nrow=k, ncol=number_of_network)


  if (!missing(last_out)){
    ystar1 = last_out$ystar1
    ystar2 = last_out$ystar2
    delta = last_out$delta
    Sigma = last_out$Sigma
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
        ystar1[,j] = drawYstar_single(y=y[,j] , ystar_other=ystar2[,j], mean=xb1[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )

        temp = find_normal_conditional_dist(a= ystar2_demean, i=j, j=-j, Sigma=Sigma)
        ystar2[,j] = drawYstar_single(y=y[,j] , ystar_other=ystar1[,j], mean=xb2[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )
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
      Sigma = solve( rwish(n , solve( crossprod(ystar_demean )) ) )

      normalization = diag(1/sqrt(diag(Sigma)))

      Sigma = normalization %*%  Sigma %*% t(normalization) 
      Sigma_matrix[i,] = Sigma[lower.tri(Sigma)]
  }   
  toc()

  out = list(delta_matrix=delta_matrix, ystar1=ystar1,ystar2=ystar2, Sigma=Sigma,Sigma_matrix=Sigma_matrix, delta=delta)

  class(out) = "network_formation" 
  out
}





# multi_network_formation_mcmc_RE = function(data, m=1000, last_out){

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

#   if (!missing(last_out)){
#     delta1_matrix[1,] = tail(last_out$delta1,1)
#     delta2_matrix[1,] = tail(last_out$delta2,1)
#     rho_matrix[1,] = tail(last_out$rho,1)
#     ystar11 = last_out$ystar11
#     ystar12 = last_out$ystar12
#     ystar21 = last_out$ystar21
#     ystar22 = last_out$ystar22
#     e1 = last_out$e1
#     e2 = last_out$e2
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



merge.network_formation = function(x,y){
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

getParameterMatrix.network_formation = function(x){
  if (is.list(x$delta_matrix)){
    out = do.call(cbind, x$delta_matrix)
    out = cbind(out, x$Sigma_matrix )
  } else{
    out = x$delta_matrix
  }
  out
}




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







