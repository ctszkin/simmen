# drawingPhi_multi = function(lambda, rho_e1, rho_e2 , data_matrix){
#   new_y = data_matrix$Y - lambda[1] * data_matrix$WY1 - lambda[2] * data_matrix$WY2 - rho_e1 - rho_e2
#   new_x = data_matrix$X

#   lm_fit = myFastLm(X=new_x, y=new_y)

#   mvrnorm(n=1, mu=lm_fit$coefficients, Sigma=lm_fit$cov)
# })

#' drawingPhi
#' @name drawingPhi
#' @aliases drawingPhi
#' @title drawingPhi
#' @detail drawingPhi
#' @param lambda lambda
#' @param rho_e rho_e
#' @param data_matrix data_matrix
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
drawingPhi= function(lambda, rho_e, data_matrix){
  new_y = data_matrix$Y - rho_e - data_matrix$WY %*% lambda
  # for (i in 1:length(lambda))
  #   new_y = new_y - lambda[i] * data_matrix$WY[,i]

  lm_fit = myFastLm(X=data_matrix$X, y=new_y)
  mvrnorm(n=1, mu=lm_fit$coefficients, Sigma=lm_fit$cov)
}

# data_matrix = genDataMatrix(data)
# benchmark(replications=1000,
# drawingPhi(parameter$lambda,parameter$e1,parameter$e2,data_matrix)
# )

#' logLikeliMultiNormalSigma
#' @name logLikeliMultiNormalSigma
#' @aliases logLikeliMultiNormalSigma
#' @title logLikeliMultiNormalSigma
#' @detail logLikeliMultiNormalSigma
#' @param x x
#' @param sigma sigma
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
logLikeliMultiNormalSigma = function(x, sigma){
  k = ncol(x)
  Sigma = diag(k)
  Sigma[lower.tri(Sigma)] = sigma
  Sigma[upper.tri(Sigma)] =  t(Sigma)[upper.tri(Sigma)]


  sum ( dmvnorm(x, sigma= Sigma ,log=T) )
}


# Likelihood function of the outcome equation of all different models.

#' loglikelihood_outcome_vector_multi
#' @name loglikelihood_outcome_vector_multi
#' @aliases loglikelihood_outcome_vector_multi
#' @title loglikelihood_outcome_vector_multi
#' @detail Likelihood function of the outcome equation in vector form
#' @param e_matrix e_matrix
#' @param new_e new_e
#' @param lambda lambda
#' @param phi phi
#' @param Y Y
#' @param X X
#' @param WY WY
#' @param network_id network_id
#' @param rho rho
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
loglikelihood_outcome_vector_multi = function(e_matrix, new_e, lambda, phi, Y, X, WY, network_id, rho){
  n = length(Y)
  if (!is.matrix(e_matrix)){
    e_matrix = as.matrix(e_matrix)
  }

  new_y = Y - WY %*% lambda

  e_matrix_new = e_matrix
  e_matrix_new[,network_id] = new_e

  e_matrix = scale(e_matrix,scale=FALSE)
  e_matrix_new = scale(e_matrix_new,scale=FALSE)

  x_phi = X%*% phi

  varepsilon_old = new_y - x_phi -  e_matrix %*% rho
  varepsilon_new = new_y - x_phi -  e_matrix_new %*% rho

  sigma2 = sum(varepsilon_old^2)/(n-1)

  out = (- varepsilon_new^2 + varepsilon_old^2) / 2 / sigma2
  out[!is.finite(out)] = -1e+20

  out
}

#' parser of likelihood function to rho 
#' @name lik_multi_exogenous_parser_rho
#' @aliases lik_multi_exogenous_parser_rho
#' @title lik_multi_exogenous_parser_rho
#' @detail lik_multi_exogenous_parser_rho
#' @param lambda lambda
#' @param phi phi
#' @param X X
#' @param Y Y
#' @param WY WY
#' @param W_list W_list
#' @param n_vector n_vector
#' @param e e
#' @param rho rho
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
lik_multi_exogenous_parser_rho = function(lambda, phi, X, Y, WY, W_list, n_vector, e, rho){

  lik_multi_exogenous(lambda, phi, X, Y-  e %*% rho  , WY,  W_list, n_vector)

}

#' update the value of e
#' @name update_e_multi
#' @aliases update_e_multi
#' @title update_e_multi
#' @detail update_e_multi
#' @param data data
#' @param e_matrix e_matrix
#' @param delta delta
#' @param network_data network_data
#' @param tau tau
#' @param lambda lambda
#' @param rho rho
#' @param phi phi
#' @param network_id network_id
#' @param ystar1 ystar1
#' @param ystar2 ystar2
#' @param Sigma Sigma
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
update_e_multi = function(data, e_matrix, delta, network_data, tau, lambda, rho, phi, network_id, ystar1, ystar2, Sigma){

  # e_dist = match.arg(e_dist)

  Y = data$y
  X = data$x_wx
  WY = data$wy

  n = data$n

  e = e_matrix[,network_id]



  e_i = e[data$group_index[,1]]
  e_j = e[data$group_index[,2]]

  diff_e = transform_e_by_each(data, e)


  xb1 = data$self_data_matrix %*% head(delta,-1)
  xb2 = data$friends_data_matrix %*% head(delta,-1)

  # this is not rho
  # it should be find_normal_conditional_dist(i=network_id,j=-network_id,Sigma=Sigma)$var
  sigma2= find_normal_conditional_dist(i=network_id,j=-network_id,Sigma=Sigma)$var
  prob_old1 = dnorm(ystar1 - xb1 - diff_e * tail(delta,1), log=TRUE, sd = sqrt(sigma2)) 
  prob_old2 = dnorm(ystar2 - xb2 - diff_e * tail(delta,1), log=TRUE, sd = sqrt(sigma2)) 
  # prob_old1 = dnorm(ystar1 - xb1 - diff_e * tail(delta,1), log=TRUE, sd = sqrt(1-rho)) 
  # prob_old2 = dnorm(ystar2 - xb2 - diff_e * tail(delta,1), log=TRUE, sd = sqrt(1-rho)) 


  # poly_e = genPoly(e_i,e_j, tail(delta,1))
  # p_R = pnorm(R_ij+poly_e[,1]) * pnorm(R_ji+poly_e[,2])
  # prob_old = log( p_R^response * (1-p_R)^(1-response) )

  # if (e_dist=="normal"){
    # normal dist
    new_e = e + rnorm(data$n, sd=tau )
    prob_diff = dnorm.diff(new_e,e) + loglikelihood_outcome_vector_multi(e_matrix=e_matrix, new_e=new_e, lambda=lambda, phi=phi,  Y=Y, X=X, WY=WY, network_id=network_id, rho=rho)
  # } else if (e_dist=="uniform"){
  #   # uniform
  #   new_e = e + runif(data$n, min=-tau,max=tau )

  #   prob_diff = ifelse(new_e>1 | new_e<0, -Inf,  0) + loglikelihood_outcome_vector_multi(e_matrix=e_matrix-mean(e_matrix), new_e=new_e, lambda=lambda, phi=phi,  Y=Y, X=X, WY=WY, network_id=network_id, rho=rho)
  # }


  uniform_rv = runif(n);


  update_e_internal_single(
    location_index_all=network_data$location_index_all,
    location_index1=network_data$location_index1,
    location_index2=network_data$location_index2,
    e=e,
    new_e=new_e,
    xb1=xb1,
    xb2=xb2,
    delta=tail(delta,1),
    prob_old1=prob_old1,
    prob_old2=prob_old2,
    ystar1=ystar1 ,
    ystar2=ystar2,
    e_i=e_i,
    e_j=e_j,
    uniform_rv=uniform_rv,
    prob_diff=prob_diff,
    sigma2=sigma2
  )

  # for (i in 1:n){
  #   location_index_all = network_data[[i]]$location_index_all
  #   location_index1 = network_data[[i]]$location_index1
  #   location_index2 = network_data[[i]]$location_index2

  #   e_i_new = e_i
  #   e_j_new = e_j

  #   e_i_new[location_index1] = new_e[i] 
  #   e_j_new[location_index2] = new_e[i] 

  #   # response_new = response[location_index_all]

  #   # poly_e_new  = genPoly(e_i_new[location_index_all],e_j_new[location_index_all], tail(delta,1))
  #   # p_R_new = pnorm(R_ij[location_index_all,]+poly_e_new[,1]) * pnorm(R_ji[location_index_all,]+poly_e_new[,2])
  #   # prob_new = log( p_R_new^response_new * (1-p_R_new)^(1-response_new) )

  #   diff_e_new = (e_i_new[location_index_all]- e_j_new[location_index_all])^2
  #   Ri_new = xb1[location_index_all,]+diff_e_new*tail(delta,1)
  #   Rj_new = xb2[location_index_all,]+diff_e_new*tail(delta,1)

  #   ystar1_new = ystar1[location_index_all]
  #   ystar2_new = ystar2[location_index_all]

  #   prob_new1 = dnorm(ystar1_new-Ri_new, log=TRUE, sd = sqrt(1-rho^2))
  #   prob_new2 = dnorm(ystar2_new-Rj_new, log=TRUE, sd = sqrt(1-rho^2))

  #   prob_new = sum(prob_new1,prob_new2)

  #   prob_old = sum(prob_old1[location_index_all]) + sum(prob_old2[location_index_all])

  #   alpha = exp( prob_new - prob_old + prob_diff[i] )

  #   if (!is.finite(alpha))
  #     alpha=0

  #   if ( alpha > runif(1)){
  #     e_i = e_i_new
  #     e_j = e_j_new 
  #     prob_old1[location_index_all] = prob_new1
  #     prob_old2[location_index_all] = prob_new2
  #     update_count[i] = TRUE
  #   }

  # }

  # e[update_count] = new_e[update_count]
  # return(list(e=e, update_rate= mean(update_count)))
}


# update_e_multi_parallel = function(cl, e1, e2, delta, network_data, tau, lambda, phi, network_id,G, rho, ystar1, ystar2,phi_rho){
#   parLapply(cl, 1:G, 
#     function(i, e1, e2, delta, lambda, phi, network_id, tau,rho,ystar1,ystar2,phi_rho){
#       update_e_multi(data=data[[i]], e1=e1[[i]], e2=e2[[i]], delta=delta, network_data=network_data[[i]], tau=tau[i], lambda=lambda, phi=phi, network_id=network_id,rho=rho, ystar1=ystar1[[i]], ystar2=ystar2[[i]] ,phi_rho=phi_rho)
#     },
#     delta=delta, 
#     tau, 
#     e1=e1,
#     e2=e2,
#     lambda=lambda, 
#     phi=phi, 
#     network_id=network_id,
#     rho=rho,
#     phi_rho=phi_rho,
#     ystar1=ystar1,
#     ystar2=ystar2
#   )

# })

# update_e_multi_parser = function(data, network_data, e_list, delta, tau, lambda, rho, phi, network_id,ystar1_demean_list, ystar2_demean_list, Sigma){

#   out = vector("list",length(data))
#   for (i in 1:length(out)){
#     out[[i]] = update_e_multi(data=data[[i]], e_matrix = e_list[[i]], delta=delta, network_data=network_data[[i]], tau = tau[[i]], lambda = lambda, rho= rho, phi= phi, network_id=network_id, ystar1=ystar1_demean_list[[i]], ystar2= ystar2_demean_list[[i]],Sigma=Sigma)
#   }
#   out 
# })

# countUpdateRate = function(x){
#   # c(
#   # lambda=uniqueRow(x$lambda_matrix),
#   # delta1=uniqueRow(x$delta1_matrix),
#   # delta2=uniqueRow(x$delta2_matrix),

#   # e=x$update_rate$e / nrow(x$lambda_matrix)
#   # # ui=x$update_rate$ui / nrow(x$lambda_matrix),
#   # # uj=x$update_rate$uj / nrow(x$lambda_matrix)
#   # )

#   c(
#   lambda=uniqueRow(x$lambda_matrix),
#   delta=mean(uniqueRow(x$delta1_matrix),uniqueRow(x$delta2_matrix)),

#   e= mean(x$update_rate$e1 , x$update_rate$e2)

#   )
# }




# initialization_multi = function(data){

#   n2 = sapply(data, function(z) length(z$response[[1]]))

#   data_matrix = genDataMatrix(data,"all")
#   # some initialization
#   k_phi = ncol( data_matrix$X )
#   k_lambda = 2

#   # coefficient for (ei-ej)^2, no coef for ei
#   k_delta1 = ( ncol(data[[1]]$self_data_matrix) + 1) 
#   k_delta2 = ( ncol(data[[1]]$self_data_matrix) + 1) 
#   k_phi_rho = 2
#   k_rho = 1

#   k_e1 = data_matrix$n_vector
#   k_e2 = data_matrix$n_vector


#   # define the initial value
#   parameter=list()
#   parameter$phi = rep(0,k_phi)
#   parameter$lambda = rep(0,k_lambda)
#   parameter$delta1 = rep(0,k_delta1)
#   parameter$delta2 = rep(0,k_delta2)
#   parameter$phi_rho_matrix = rep(0,k_phi_rho)
#   parameter$rho_matrix = rep(0,k_rho)


#   parameter$e1 =lapply(k_e1, function(z) rnorm(z))
#   parameter$e2 =lapply(k_e2, function(z) rnorm(z))

#   parameter$ystar11 = rep(0,sum(n2))
#   parameter$ystar12 = rep(0,sum(n2))
#   parameter$ystar21 = rep(0,sum(n2))
#   parameter$ystar22 = rep(0,sum(n2))

#   parameter$tau = list(lambda=0.3, rho=0.3, e1=rep(0.5,length(data)), e2=rep(0.5,length(data)))
#   parameter
# })


# createVariableName_multi = function(data,...){
#   c("phi_" %+% colnames(data[[1]]$x_wx)
#     ,c("phi_rho1","phi_rho2")
#     , paste0("lambda_",1:2)
#     ,"delta1_" %+%  c(colnames(data[[1]]$self_data_matrix),  "e_diff")
#     , "delta2_" %+%  c(colnames(data[[1]]$self_data_matrix),  "e_diff" )
#     , "rho"
#   )
# }


# getLastParameter_multi = function(x){
#   parameter=list()
#   parameter$phi = as.vector(tail(x$phi,1))
#   parameter$lambda = as.vector(tail(x$lambda,1))
#   parameter$delta1 = as.vector(tail(x$delta1,1))
#   parameter$delta2 = as.vector(tail(x$delta2,1))
#   parameter$rho = as.vector(tail(x$rho,1))
#   parameter$phi_rho = as.vector(tail(x$phi_rho,1))

#   parameter$e1 = x$e1
#   parameter$e2 = x$e2
#   parameter$ystar11 = x$ystar11
#   parameter$ystar12 = x$ystar12
#   parameter$ystar21 = x$ystar21
#   parameter$ystar22 = x$ystar22

#   tau=x$tau
#   for ( j in 1:length(tau)){
#     if (any(x$update_rate[[j]] >0.4 | any(x$update_rate[[j]]< 0.3) )){
#       cat("update tau-", names(tau)[[j]] , "\n")
#       tau[[j]] = tau[[j]] * x$update_rate[[j]] / 0.35
#       tau[[j]] = ifelse(tau[[j]]==0, 0.001, tau[[j]])
#     }

#   }
#   parameter$tau = tau

#   parameter


# }



#' Update the sign of e
#' @name update_sign_e_multi
#' @aliases update_sign_e_multi
#' @title update_sign_e_multi
#' @detail update_sign_e_multi
#' @param lambda lambda
#' @param phi phi
#' @param rho rho
#' @param data_matrix data_matrix
#' @param e_list e_list
#' @param network_id network_id
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
update_sign_e_multi = function(lambda,phi,rho, data_matrix, e_list, network_id){
  n = data_matrix$n_vector


  new_y = data_matrix$Y - data_matrix$WY %*% lambda  - data_matrix$X %*% phi 

  # unlist_e = unlist(e)
  demean_e_list = lapply(e_list, scale, scale=F)
  unlist_e = do.call(rbind, demean_e_list)

  unlist_e_new = unlist_e
  unlist_e_new[,network_id] = -unlist_e_new[,network_id]

  ssr_old = (new_y - (unlist_e %*% rho))^2
  ssr_new = (new_y - (unlist_e_new %*% rho))^2

  ssr_old_list = splitBy(ssr_old,n)
  ssr_new_list = splitBy(ssr_new,n)

  sum_ssr_old = sapply(ssr_old_list , sum )
  sum_ssr_new = sapply(ssr_new_list , sum)

  flip_indicator = which( sum_ssr_old > sum_ssr_new )

  if (length(flip_indicator)>0){
    # cat("flip sign :")
    for (i in flip_indicator ){
      # cat("Group=", i,";network=",network_id,";diff in SSR=", sum_ssr_old[i]-sum_ssr_new[i]," | ")
      e_list[[i]][,network_id] = -e_list[[i]][,network_id]
    }
    # cat("\n")
  }
  list(e_list=e_list, update=sum_ssr_old > sum_ssr_new)
}


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

#' Social Interaction Model with Multiple Endogenous Networks
#' @name simmen
#' @aliases simmen
#' @title Social Interaction Model with Multiple Endogenous Networks
#' @detail Social Interaction Model with Multiple Endogenous Networks
#' @param m Number of iteration
#' @param data data
#' @param last_out previous estimation
#' @param allow_correlation Logical. Allow correlation between network
#' @param update_sign Logical. Whether to update sign of e
#' @param initial_seed seed for the first round
#' @param seed_for_chain generate seed for each chain
#' @param start_value_delta_diff_e start value of the coefficient of abs(ei-ej)
#' @param tau_e sd of sampling distribution
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
simmen = function(m, data, last_out, allow_correlation=TRUE, update_sign=TRUE, initial_seed = 100, seed_for_chain = round(runif(1)*10000),   start_value_delta_diff_e  = 0 , tau_e=1){
  # e_dist = match.arg(e_dist)

  number_of_network = length(data[[1]]$D_list)

  if (number_of_network==1){
    # return( single_endogenous_mcmc(m=m,data=data,last_out=last_out) )
    allow_correlation=FALSE
  }

  network_data = genNetworkData(data)
  data_matrix = genDataMatrix(data)

  n = sapply(data, "[[", "n")
  n2 = sapply(data, function(z) NROW(z$response_self ))
  nn = sum(n2)
  G = length(data)

  self_data_matrix = do.call(rbind, lapply(data, "[[", i="self_data_matrix" ))
  friends_data_matrix = do.call(rbind, lapply(data , "[[", i="friends_data_matrix"))

  y = do.call(rbind, lapply(data, "[[", "response_self"))
  y_not = !y

  delta = rbind( matrix(0, nrow=NCOL(self_data_matrix), ncol=number_of_network), start_value_delta_diff_e)
  Sigma = diag(number_of_network)
  phi = rep(0, ncol(data_matrix$X))
  lambda = rep(0, number_of_network)
  rho = rep(0, number_of_network)

  network_name = data[[1]]$network_name

  number_col_Sigma_matrix = number_of_network*(number_of_network-1)/2
  sigma_index = genPairwiseIndex(number_of_network)
  sigma_name = network_name[sigma_index[,1]] %+% "_" %+% network_name[sigma_index[,2]]


  ystar1 = array(0, dim= dim(y))
  ystar2 = array(0, dim= dim(y))

  # tmp_tau_e = rep(list(rep(0.5,G)) , number_of_network)
  tmp_tau_e = lapply(n, function(z)  matrix(tau_e, nrow=z, ncol=number_of_network) )
  tmp_tau_e = lapply(tmp_tau_e, function(z) {colnames(z) = network_name; z} )

  tau = list(e=tmp_tau_e,lambda=0.3, rho=0.3, sigma_matrix=rep(0.1,number_col_Sigma_matrix))
  flip_sign_rate = matrix(0, nrow=number_of_network, ncol=G)

  colnames(flip_sign_rate) =names(data)
  rownames(flip_sign_rate) = network_name


  index = 1
  iteration_index  = c(start=1, end=m)
  ID = genUniqueID()
  # cat("Start new instance with ID ", ID, "\n")

  if (!missing(last_out)){
    if (missing(allow_correlation)){
      allow_correlation = last_out$allow_correlation
    }
    if (missing(update_sign)){
      update_sign = last_out$update_sign
    }

    delta = sapply(last_out$delta, tail, n=1)
    Sigma = last_out$Sigma
    phi = as.vector(tail(last_out$phi,1))
    lambda = as.vector(tail(last_out$lambda,1))
    rho = as.vector(tail(last_out$rho,1))

    e_list = last_out$e_list

    ## remove later: converge e_list to the new one
    # e_list = lapply(e_list, function(z) {colnames(z) = network_name; z} )


    ystar1=last_out$ystar1
    ystar2=last_out$ystar2

    ## remove later: converge tau to the new one
    tau=updateTau(last_out$tau, last_out$update_rate)
    # names(tau) = c("lambda","rho","e_" %+% network_name)

    index = last_out$index+1
    iteration_index = last_out$iteration_index[2] + c(start=1,end=m)


    initial_seed = last_out$initial_seed
    seed_for_chain = last_out$seed_for_chain
    # name = last_out$name
    ID = last_out$ID
  } else{
    set.seed(initial_seed)
    e_list = lapply(n, function(z)  matrix(rnorm(z*number_of_network), nrow=z, ncol=number_of_network) )
    e_list = lapply(e_list, function(z) {colnames(z) = network_name; z} )

  }
  set.seed(seed_for_chain[index])




  demean_e_list = lapply(e_list, scale, scale=F)
  demean_e = do.call(rbind,demean_e_list) 

  phi_matrix = matrix(NA, nrow=m, ncol=length(phi) )
  lambda_matrix = matrix(NA, nrow=m, ncol=length(lambda))
  rho_matrix = matrix(NA, nrow=m, ncol=length(rho))

  Sigma_matrix = matrix(NA, nrow=m, ncol = number_col_Sigma_matrix )

  colnames(phi_matrix) = c("phi_" %+% colnames(data[[1]]$x_wx))
  colnames(rho_matrix) = "rho_" %+% network_name
  colnames(lambda_matrix) = "lambda_" %+% network_name

  if (allow_correlation){
    colnames(Sigma_matrix) = "Sigma_" %+% sigma_name
  }

  delta_matrix = vector("list",number_of_network)

  tmp_update_rate = lapply(n, function(z)  matrix(0, nrow=z, ncol=number_of_network) )
  tmp_update_rate = lapply(tmp_update_rate, function(z) {colnames(z) = network_name; z} )



  update_rate =list(e=tmp_update_rate,lambda=0,rho=0, sigma_matrix=rep(0,number_col_Sigma_matrix))
  for (j in 1:number_of_network){
    delta_matrix[[j]] = matrix(NA, nrow=m, ncol=nrow(delta))
    colnames(delta_matrix[[j]]) = "delta_" %+% network_name[j] %+% "_" %+% c(colnames(self_data_matrix), "diff_e")
    # update_rate["e" %+% "_" %+% network_name[j]] = 0
  }



  diff_e = do.call(rbind, transform_e2(data,e_list) )

  R1 = self_data_matrix %*% head(delta,-1)
  R2 = friends_data_matrix %*% head(delta,-1)
  xb1 = R1 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1))) 
  xb2 = R2 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1))) 

  ystar1_demean = ystar1 - xb1
  ystar2_demean = ystar2 - xb2
  ystar = rbind(ystar1,ystar2)


  Z = rbind(self_data_matrix,friends_data_matrix)

  tic()
  for (i in 1:m ){
    ##### update ystar
    for (j in 1:number_of_network){
      temp = find_normal_conditional_dist(a= ystar1_demean, i=j, j=-j, Sigma=Sigma)

      ystar1[,j] = drawYstar(y=y[,j] , ystar_other=ystar2[,j], mean=xb1[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )

      temp = find_normal_conditional_dist(a= ystar2_demean, i=j, j=-j, Sigma=Sigma)

      ystar2[,j] = drawYstar(y=y[,j] , ystar_other=ystar1[,j], mean=xb2[,j] + temp$mean, y_not= y_not[,j], sd= sqrt(temp$var) )

      ystar1_demean = ystar1 - xb1
      ystar2_demean = ystar2 - xb2
      ystar_demean = rbind(ystar1_demean,ystar2_demean)
      ystar = rbind(ystar1,ystar2)

    }
    ### update e
    for (j in 1:number_of_network){
      temp = find_normal_conditional_dist(a= ystar1_demean, i=j, j=-j, Sigma=Sigma)

      ystar1_demean_list = splitBy(ystar1[,j]-temp$mean, n2)

      temp = find_normal_conditional_dist(a= ystar2_demean, i=j, j=-j, Sigma=Sigma)

      ystar2_demean_list = splitBy(ystar2[,j]-temp$mean, n2)
    
      metro_e1 = vector("list",G)
      for (ii in 1:G){
        metro_e1[[ii]] = update_e_multi(data=data[[ii]], e_matrix = e_list[[ii]], delta=delta[,j], network_data=network_data[[ii]], tau = tau$e[[ii]][,j], lambda = lambda, rho= rho, phi= phi, network_id=j, ystar1=ystar1_demean_list[[ii]], ystar2= ystar2_demean_list[[ii]],Sigma=Sigma)

        e_list[[ii]][,j] = metro_e1[[ii]]$e
        update_rate$e[[ii]][,j] = update_rate$e[[ii]][,j] + metro_e1[[ii]]$update_rate 
      }

      # update_rate[["e_" %+% network_name[j] ]] = update_rate[["e_" %+% network_name[j]]] + unlist(lapply(metro_e1, function(z) z$update_rate))

      if (update_sign){
        update_e1 = update_sign_e_multi(lambda=lambda,phi=phi,rho=rho,data_matrix=data_matrix, e_list=e_list, network_id=j)
        e_list = update_e1$e
        flip_sign_rate[j,] = update_e1$update + flip_sign_rate[j,] 
      }

      demean_e_list = lapply(e_list, scale, scale=F)
      demean_e = do.call(rbind,demean_e_list)

      diff_e = do.call(rbind, transform_e2(data,e_list) )

      xb1 = R1 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1))) 
      xb2 = R2 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1))) 
      ystar1_demean = ystar1 - xb1
      ystar2_demean = ystar2 - xb2
      ystar_demean = rbind(ystar1_demean,ystar2_demean)
    }



    
    ##### update delta
    # XX = cbind(Z, rbind(diff_e,diff_e))

    for ( j in 1:number_of_network){
      XX = cbind(Z, rep(diff_e[,j],2))
      # XX = Z
      temp = find_normal_conditional_dist(a=ystar_demean, i=j, j=-j, Sigma=Sigma)
      # lm_fit = myFastLm(X=XX, y =ystar[,j]-temp$mean + rep(diff_e[,j],2)*1/sqrt(2) )
      lm_fit = myFastLm(X=XX, y =ystar[,j]-temp$mean )

      delta[,j] = mvrnorm(n=1, mu=lm_fit$coef, lm_fit$cov/lm_fit$s^2 * as.vector(temp$var) )

      delta_matrix[[j]][i,] = delta[,j]

      R1 = self_data_matrix %*% head(delta,-1)
      R2 = friends_data_matrix %*% head(delta,-1)
      xb1 = R1 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1))) 
      xb2 = R2 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1))) 

      ystar1_demean = ystar1 - xb1
      ystar2_demean = ystar2 - xb2
      ystar_demean = rbind(ystar1_demean,ystar2_demean)
      ystar = rbind(ystar1,ystar2)
    }






    ##### update phi
    phi = drawingPhi(
      lambda=lambda,
      rho_e= demean_e %*% rho,
      data_matrix=data_matrix
    )
    phi_matrix[i,] = phi

    ##### update lambda
    metro_obj = metropolis2(
      beta_previous=lambda, 
      tau=tau$lambda, 
      likelihoodFunction=lik_multi_exogenous, 
      prior_type="normal", 
      sampling_type="normal", 
      phi = phi,
      W_list = data_matrix$W_list,
      X = data_matrix$X,
      Y = data_matrix$Y-    demean_e %*% rho ,
      WY = data_matrix$WY,
      n_vector = n
    )
    lambda = metro_obj$beta
    update_rate$lambda = update_rate$lambda + metro_obj$update
    lambda_matrix[i,] = lambda

    ##### update rho_matrix
    # if (e_dist == "normal"){
      metro_obj = metropolis2(
        beta_previous=rho, 
        tau=tau$rho, 
        likelihoodFunction=lik_multi_exogenous_parser_rho, 
        prior_type="trunc_normal", 
        sampling_type="trunc_normal",
        lambda = lambda, 
        phi = phi,
        W_list = data_matrix$W_list,
        X = data_matrix$X,
        Y = data_matrix$Y,
        WY = data_matrix$WY,
        n_vector = data_matrix$n_vector,
        e=demean_e
      )
      rho = metro_obj$beta
      update_rate$rho = update_rate$rho + metro_obj$update
      rho_matrix[i,] = rho
    # } else{
    #   ## better to combine this to phi
    #   new_y = data_matrix$Y - rowSums(data_matrix$WY %*% diag(lambda))  - data_matrix$X %*% phi 
    #   lm_fit = myFastLm(X=demean_e, y = new_y )
    #   rho = mvrnorm(n=1, mu=lm_fit$coefficients, Sigma=lm_fit$cov)
    #   rho_matrix[i,]  = rho
    # }

    ##### update correlation
    if (allow_correlation){

      Sigma_metro = metropolis2(
        beta_previous=Sigma[lower.tri(Sigma)], 
        tau=tau$sigma_matrix, 
        likelihoodFunction=logLikeliMultiNormalSigma, 
        prior_type = "uniform",
        x=ystar_demean
      )

      Sigma[lower.tri(Sigma)] =  Sigma_metro$beta

      Sigma[upper.tri(Sigma)] =  t(Sigma)[upper.tri(Sigma)]


      update_rate$sigma_matrix = Sigma_metro$update + update_rate$sigma_matrix
      Sigma_matrix[i,] = Sigma[lower.tri(Sigma)]

      # Sigma = solve( rwish(nn*2 , solve( crossprod(ystar_demean -colMeans(ystar_demean))) ) )
      # normalization = diag(1/sqrt(diag(Sigma)))

      # Sigma = normalization %*%  Sigma %*% t(normalization) 
      # Sigma_matrix[i,] = Sigma[lower.tri(Sigma)]
    }
    # cat(i, " ",cor(ystar_demean)[2], "\n")
    # cat(i, "sigma=", Sigma_matrix[i,], " | delta=", tail(delta,1), " | rho=", rho, "\n")
  }
  toc()
  # plot(Sigma_matrix,type='l')
  # Sigma



  out=list()
  out$allow_correlation=allow_correlation
  out$update_sign=update_sign
  out$seed = seed_for_chain
  out$initial_seed = initial_seed
  out$m = m 

  out$phi=phi_matrix 
  out$lambda=lambda_matrix
  out$delta=delta_matrix
  out$rho = rho_matrix
  out$e_list=e_list
  out$Sigma=Sigma
  out$Sigma_matrix=Sigma_matrix

  out$ystar1=ystar1
  out$ystar2=ystar2

  out$update_rate = rapply(update_rate, function(z) z/m , how = "replace")
  out$flip_sign_rate = flip_sign_rate / m 

  out$tau=tau

  out$index=index
  out$iteration_index= iteration_index
  out$ID=ID

  out$initial_seed = initial_seed
  out$seed_for_chain = seed_for_chain



  cat("Summary of update rate:\n")
  for (i in 1:length(out$update_rate)) {
    if (is.list( out$update_rate[[i]])){
      cat("$",names(out$update_rate[i]), "\n")
      tmp = do.call(rbind, out$update_rate[[i]])
      print(apply(tmp,2,quantile))
      cat("\n")
    } else{
      print(out$update_rate[i])
    }
  }

  cat( "sigma=", tail(Sigma_matrix,1), " | delta=", tail(delta,1), " | rho=", rho, "\n")

  cat("Summary of flip sign:")
  print(out$flip_sign_rate)
  cat("\n\n")
  class(out) = "simmen"
  out

}


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################



#' merge.simmen
#' @name merge.simmen
#' @aliases merge.simmen
#' @title merge.simmen
#' @detail merge.simmen
#' @param x x
#' @param y y 
#' @param ... ... 
#' @method merge simmen
#' @S3method merge simmen
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
merge.simmen = function(x,y,...){
  if (length(list(...))>0){
    list_args = list(...)
    out = list_args[[1]]
    for (i in 2:length(list_args)){
      out = merge(out, list_args[[i]])
    }
    return(out)
  }
  if (y$index < x$index){
    temp = x
    x = y
    y = temp
    rm(temp)
  }
  out = y
  out$phi = rbind(x$phi,y$phi)
  out$phi_rho = rbind(x$phi_rho,y$phi_rho)
  out$lambda = rbind(x$lambda,y$lambda)
  i=NULL
  j=NULL
  out$delta = foreach(i=x$delta, j=y$delta ) %do% {
    rbind(i,j)
  }

  out$Sigma_matrix = rbind(x$Sigma_matrix,y$Sigma_matrix)
  out$rho = rbind(x$rho,y$rho)
  out
}



#' plotmcmc2.simmen
#' @name plotmcmc2.simmen
#' @aliases plotmcmc2.simmen
#' @title plotmcmc2.simmen
#' @detail plotmcmc2.simmen
#' @param x x
#' @param name name
#' @param file filename of the plot, Default is "mcmc_plot.pdf"
#' @param tail How many observation of tail to use. tail<0 means remove first -tail iteration. -1<tail<1 means percentage 
#' @param ... ...
#' @return value
#' @method plotmcmc2 simmen
#' @S3method plotmcmc2 simmen
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
plotmcmc2.simmen = function(x, name, file="mcmc_plot.pdf", tail=10000,...){
  xx = lapply(x, getParameterMatrix)
  plotmcmc2.default(xx, name=name, tail=tail )
}

#' Get a matrix of all parameters
#' @name getParameterMatrix.simmen
#' @aliases getParameterMatrix.simmen
#' @title getParameterMatrix.simmen
#' @detail getParameterMatrix.simmen
#' @param x x
#' @param tail tail
#' @param ... ...
#' @method getParameterMatrix simmen
#' @S3method getParameterMatrix simmen
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
getParameterMatrix.simmen = function(x,tail,...){
  delta_matrix = do.call(cbind, x$delta)

  out = cbind(x$phi,x$rho,x$lambda,delta_matrix,x$Sigma_matrix)
  # colnames(out) = x$name
  if (!missing(tail) & is.numeric(tail) ) {
    out = tail(out,tail)
  }
  out
}



#' plotmcmc.simmen
#' @name plotmcmc.simmen
#' @aliases plotmcmc.simmen
#' @title plotmcmc.simmen
#' @detail plotmcmc.simmen
#' @param x x
#' @param tail tail
#' @param ... ...
#' @method plotmcmc simmen
#' @S3method plotmcmc simmen
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
plotmcmc.simmen = function(x, tail=-0.2,...){
  #x11(20,20)
  plotmcmc(cbind(x$phi,x$rho,x$lambda), tail=tail)
  #x11(20,20)
  plotmcmc.default(cbind(do.call(cbind, x$delta),x$Sigma_matrix), tail=tail )
}


# plotmcmc3 = function(x,tail){
#   tmp = cbind(x$rho, x$delta[[1]][,27+4], x$delta[[2]][,27+4] )
#   colnames(tmp) = c( colnames(x$rho), colnames(x$delta[[1]])[27+4], colnames(x$delta[[2]])[27+4] )
#   if (!missing(tail)){
#     tmp = tail(tmp,tail)
#   }
#   plotmcmc.default(tmp)
# }


# plotmcmc4 = function(x,tail=NA){
#   parameter_obj = lapply(x, function(z) {
#     out =cbind(z$rho, z$delta[[1]][,27+4], z$delta[[2]][,27+4] )
#     colnames(out) = c( colnames(z$rho), colnames(z$delta[[1]])[27+4], colnames(z$delta[[2]])[27+4] )
#     if (!is.na(tail)){
#       out = tail(out,tail)
#     }
#     out
#   }
#   )

#   par(mfrow=c(2,2))
#   for (i in 1:4){
#     tmp = foreach( j=parameter_obj, .combine=cbind) %do% {
#       j[,i]
#     }
#     matplot(tmp,type='l')
#     title(colnames(j)[i])
#   }
# }

# plotmcmc5 = function(x,tail=NA){
#   parameter_obj = lapply(x, function(z) {
#     out =cbind(z$rho, z$delta[[1]][,27], z$delta[[2]][,27] )
#     colnames(out) = c( colnames(z$rho), colnames(z$delta[[1]])[27], colnames(z$delta[[2]])[27] )
#     if (!is.na(tail)){
#       out = tail(out,tail)
#     }
#     out=apply(out, 2, cumsum) / (1:nrow(out))
#     out
#   }
#   )

#   par(mfrow=c(2,2))
#   for (i in 1:4){
#     tmp = foreach( j=parameter_obj, .combine=cbind) %do% {
#       j[,i]
#     }
#     matplot(tmp,type='l',lty=1)
#     title(colnames(j)[i])
#   }
# }


#' Compute the likelihood basic on the parameter of last iteration.
#' @name findLikelihood
#' @aliases findLikelihood
#' @title findLikelihood
#' @detail findLikelihood
#' @param data data
#' @param last_out last_out
#' @return value
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
findLikelihood = function(data, last_out) {
  number_of_network = length(data[[1]]$D_list)
  network_data = genNetworkData(data)
  data_matrix = genDataMatrix(data)

  n = sapply(data, "[[", "n")
  n2 = sapply(data, function(z) NROW(z$response_self ))
  nn = sum(n2)
  G = length(data)


  self_data_matrix = do.call(rbind, lapply(data, "[[", i="self_data_matrix" ))
  friends_data_matrix = do.call(rbind, lapply(data , "[[", i="friends_data_matrix"))

  y = do.call(rbind, lapply(data, "[[", "response_self"))
  y_not = !y

  delta = matrix(0, nrow=NCOL(self_data_matrix)+1, ncol=number_of_network)
  Sigma = diag(number_of_network)
  phi = rep(0, ncol(data_matrix$X))
  lambda = rep(0, number_of_network)
  rho = rep(0, number_of_network)

  number_col_Sigma_matrix = number_of_network*(number_of_network-1)/2
  sigma_name = genPairwiseIndex(number_of_network)
  sigma_name = sigma_name[,1] %+% sigma_name[,2]


  network_name = data[[1]]$network_name

  delta = sapply(last_out$delta, tail, n=1)
  Sigma = last_out$Sigma
  phi = as.vector(tail(last_out$phi,1))
  lambda = as.vector(tail(last_out$lambda,1))
  rho = as.vector(tail(last_out$rho,1))

  e_list = last_out$e_list

  demean_e = scale( do.call(rbind,e_list) )
  demean_e_list = splitBy( demean_e, n)

  ystar1=last_out$ystar1
  ystar2=last_out$ystar2

  # tau=updateTau(last_out$tau, last_out$update_rate)
  # index = last_out$index+1
  # name = last_out$name
  # ID = last_out$ID
  # m=10
  # phi_matrix = matrix(NA, nrow=m, ncol=length(phi) )
  # lambda_matrix = matrix(NA, nrow=m, ncol=length(lambda))
  # rho_matrix = matrix(NA, nrow=m, ncol=length(rho))
  # Sigma_matrix = matrix(NA, nrow=m, ncol = number_col_Sigma_matrix )

  # colnames(Sigma_matrix) = "Sigma_" %+%  sigma_name

  # colnames(phi_matrix) = c("phi_" %+% colnames(data[[1]]$x_wx))
  # colnames(rho_matrix) = "rho_" %+% network_name
  # colnames(lambda_matrix) = "lambda_" %+% network_name
  # colnames(Sigma_matrix) = "Sigma_" %+% sigma_name

  # delta_matrix = vector("list",number_of_network)
  # for (i in 1:length(delta_matrix)){
  #   temp = matrix(NA, nrow=m, ncol=nrow(delta))
  #   colnames(temp) = "delta_" %+% network_name[i] %+% c(colnames(self_data_matrix), "diff_e")
  #   delta_matrix[[i]] =temp
  # }


  # update_rate =list(lambda=0,rho=0)

  # for (i in 1: (number_of_network)){
  #   update_rate["e" %+% i] = 0
  # }

  diff_e = do.call(rbind, transform_e2(data,e_list) )

  R1 = self_data_matrix %*% head(delta,-1)
  R2 = friends_data_matrix %*% head(delta,-1)
  xb1 = R1 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1)))  
  xb2 = R2 + diff_e %*% if (number_of_network==1) as.vector(tail(delta, 1)) else diag(as.vector(tail(delta, 1)))  
  
  ystar1_demean = ystar1 - xb1
  ystar2_demean = ystar2 - xb2

  F1 = pnorm(xb1)
  F2 = pnorm(xb2)

  p_y = F1*F2
  not_py = 1- p_y

  sum(log(p_y[y])) + sum(log(not_py[y_not]))

  X = rbind(self_data_matrix,friends_data_matrix)


  lik_outcome = lik_multi_exogenous(phi = c(phi,rho),
      W_list = data_matrix$W_list,
      X = cbind(data_matrix$X,demean_e),
      Y = data_matrix$Y ,
      WY = data_matrix$WY,
      lambda=lambda,
      n_vector = n
  )

  # lik_network = sum( dmvnorm(ystar1_demean, sigma=Sigma, log=T) ) + sum( dmvnorm(ystar2_demean, sigma=Sigma, log=T) )
  lik_network =   sum(log(p_y[y])) + sum(log(not_py[y_not]))
  c(lik_outcome,lik_network,lik_outcome+lik_network)
}

# lapply(loaded_files, findLikelihood ,data=data)
