# drawingPhi_single = function(lambda, rho_e, data_matrix){
#   new_y = data_matrix$Y - rho_e
#   for (i in 1:length(lambda))
#     new_y = new_y - lambda[i] * data_matrix$WY[,i]

#   lm_fit = myFastLm(X=data_matrix$X, y=new_y)
#   mvrnorm(n=1, mu=lm_fit$coefficients, Sigma=lm_fit$cov)
# })

# loglikelihood_outcome_vector_single = function(e, new_e, lambda, phi, Y, X, WY, rho){

#   n = NROW(Y)

#   new_y = Y 
#   for (i in 1:length(lambda))
#     new_y = new_y - lambda[i] * WY[,i]

#   xb = X %*% phi

#   varepsilon_old = new_y - xb - rho*e
#   varepsilon_new = new_y - xb - rho*new_e

#   sigma2 = sum(varepsilon_old^2)/(n-1)

#   out = (- varepsilon_new^2 + varepsilon_old^2) / 2 / sigma2

#   out[!is.finite(out)] = -1e+20

#   out
# })


lik_single_exogenous_parser_rho = function(lambda,phi, rho, data_matrix, e){
  lik_multi_exogenous(lambda=lambda, phi=phi, X=data_matrix$X, Y=data_matrix$Y-e*rho, WY=data_matrix$WY, W_list=data_matrix$W_list, n_vector=data_matrix$n_vector)
}


update_sign_e = function(lambda,phi,rho, data_matrix,e){
  n = data_matrix$n_vector

  new_y = data_matrix$Y - data_matrix$X %*% phi

  for (i in 1:length(lambda))
    new_y = new_y - lambda[i] * data_matrix$WY[,i]

  unlist_e = unlist(e)
  xb1 = (new_y - rho *unlist_e)^2
  xb2 = (new_y + rho *unlist_e)^2

  xb1_list = splitBy(xb1,n)
  xb2_list = splitBy(xb2,n)

  p1 = sapply(xb1_list , sum ) 
  p2 = sapply(xb2_list , sum)
  index = which(  p1-p2 >0    )

  if (length(index)>0){
    cat("flip sign :")
    for (i in index ){
      cat(i,":", p1[i]-p2[i],"  ")
      e[[i]] = -e[[i]]
    }
    cat("\n")
  }
  list(e=e, update=length(index)>0)
}



# all_data = data
# all_network_data = network_data
# all_tau = tau
# all_e = e
# all_ystar1 = ystar1
# all_ystar2 = ystar2

# data =data[[2]]
# network_data = network_data[[2]]
# tau = 0.5
# e = e[[2]]
# ystar1 = ystar1_list[[2]]
# ystar2 = ystar2_list[[2]]

# delta = delta_matrix[i,]
# lambda = lambda_matrix[i,]
# rho = rho_matrix[i,]
# phi = phi_matrix[i,]



# data = all_data
# network_data = all_network_data
# tau = all_tau
# e = all_e
# ystar1=all_ystar1
# ystar2=all_ystar2

update_e_single_new = function(data, e, delta, network_data, tau, lambda, rho, phi, ystar1, ystar2){

  Y = data$y
  X = data$x_wx
  WY = data$wy

  n = data$n
  update_count = rep(FALSE, n)

  e_i = e[data$group_index[,1]]
  e_j = e[data$group_index[,2]]

  diff_e = (e_i-e_j)^2

  xb1 = data$self_data_matrix  %*% head(delta,-1)
  xb2 = data$friends_data_matrix  %*% head(delta,-1)


  prob_old1 = dnorm(ystar1 , mean= xb1 + diff_e * tail(delta,1), log=TRUE) 

  prob_old2 = dnorm(ystar2 , mean= xb2 + diff_e * tail(delta,1), log=TRUE) 

  new_e = e + rnorm(data$n, sd=tau )

# loglikelihood_outcome_vector_single(e, ee, lambda, phi, Y, X, WY, rho)

  prob_diff = dnorm(new_e,log=TRUE) - dnorm(e,log=TRUE) + loglikelihood_outcome_vector_multi(e, new_e, lambda, phi, Y, X, WY, rho, network_id=1)
  
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
    ystar1=ystar1,
    ystar2=ystar2,
    e_i=e_i,
    e_j=e_j,
    uniform_rv=uniform_rv,
    prob_diff=prob_diff,
    sigma2=1
  )
}

   # microbenchmark({

  # for (i in 1:n){
  #   location_index_all = network_data[[i]]$location_index_all
  #   location_index1 = network_data[[i]]$location_index1
  #   location_index2 = network_data[[i]]$location_index2

  #   ## create new error vector
  #   e_i_new = e_i
  #   e_j_new = e_j

  #   ## update the error vector of the i element

  #   e_i_new[location_index1] = new_e[i] 
  #   e_j_new[location_index2] = new_e[i] 

  #   ## find out the updated response vector
  #   # response_new = response[location_index_all]

  #   ## create the poly of error
  #   # poly_e1_new  = genPoly(e_i_new[location_index1],e_j_new[location_index1], tail(delta,1))

  #   # poly_e2_new  = genPoly(e_i_new[location_index2],e_j_new[location_index2], tail(delta,1))

  #   diff_e_new = (e_i_new[location_index_all]- e_j_new[location_index_all])^2


  #   ## create the prob of y=1
  #   # p_R_new = pnorm(R_ij[location_index_all,]+poly_e_new[,1]) * pnorm(R_ji[location_index_all,]+poly_e_new[,2])
  #   # Ri_new = xb1[location_index_all,]+diff_e_new*tail(delta,1)
  #   # Rj_new = xb2[location_index_all,]+diff_e_new*tail(delta,1)

  #   ## take log
  #   # prob_new = log( p_R_new^response_new * (1-p_R_new)^(1-response_new) )

  #   # ystar1_new = ystar1[location_index_all]
  #   # ystar2_new = ystar2[location_index_all]

  #   # prob_new1 = dnorm(ystar1_new-Ri_new, log=TRUE)
  #   # prob_new2 = dnorm(ystar2_new-Rj_new, log=TRUE)


  #   prob_new1 = dnorm(ystar1[location_index_all],       mean=xb1[location_index_all,]+diff_e_new*tail(delta,1) , log=TRUE)

  #   prob_new2 = dnorm(ystar2[location_index_all],       mean=xb2[location_index_all,]+diff_e_new*tail(delta,1) , log=TRUE)


  #   prob_new = sum(prob_new1,prob_new2)

  #   prob_old = sum(prob_old1[location_index_all]) + sum(prob_old2[location_index_all])

  #   alpha = exp( prob_new - prob_old + prob_diff[i] )
  #       alpha_vec[i] = alpha

  #   if (!is.finite(alpha))
  #     alpha=0

  #   if ( alpha > uniform_rv[i]){
  #     e_i = e_i_new
  #     e_j = e_j_new 
  #     prob_old1[location_index_all] = prob_new1
  #     prob_old2[location_index_all] = prob_new2
  #     update_count[i] = TRUE
  #   }
  # }
  # update_count




  # e[update_count] = new_e[update_count]
  # list(e=e, update_rate= mean(update_count))
# })

update_e_single_parallel_new = function(cl, e, delta, network_data, tau, lambda, phi, G,ystar1_list,ystar2_list, rho){
  parLapply(cl, 1:G, 
    function(i, e, delta, lambda, phi, tau,ystar1_list,ystar2_list,rho){
      update_e_single_new(data=data[[i]], e=e[[i]], delta=delta, tau=tau[[i]], lambda=lambda, phi=phi, rho=rho,ystar1=ystar1_list[[i]],ystar2=ystar2_list[[i]] )
    },
    delta=delta, 
    tau, 
    e=e,
    lambda=lambda, 
    phi=phi, 
    rho=rho,
    ystar1_list= ystar1_list ,
    ystar2_list= ystar2_list
  )

}

update_e_single_new_parser = function(e, delta, network_data, tau, lambda, phi, G,ystar1_list,ystar2_list, rho){
  lapply( 1:G, 
    function(i, e, delta, lambda, phi,  tau,ystar1_list,ystar2_list,rho){
      update_e_single_new(data=data[[i]], e=e[[i]], delta=delta, tau=tau[[i]], lambda=lambda, phi=phi, rho=rho, ystar1=ystar1_list[[i]],ystar2=ystar2_list[[i]], network_data= network_data[[i]])
    },
    e=e,
    delta=delta, 
    lambda=lambda, 
    phi=phi, 
    tau=tau, 
    ystar1_list= ystar1_list ,
    ystar2_list= ystar2_list,
    rho=rho
  )

}

# initialization_single = function(data,network_id=1){
#   n2 = sapply(data, function(z) length(z$response[[network_id]]))

#   data_matrix = genDataMatrix(data,network_id)

#   # some initialization
#   k_phi = ncol( data_matrix$X )
#   k_lambda = 1

#   #  coefficient for (ei-ej)^2, no coef for ei
#   k_delta = ( ncol(data[[1]]$self_data_matrix) + 1) 

#   k_e = data_matrix$n_vector

#   # define the initial value
#   parameter=list()
#   parameter$tau = list(lambda=0.1, rho=0.1, e=rep(0.5,length(data)))
#   parameter$phi = rep(0,k_phi)
#   parameter$lambda = rep(0,k_lambda)
#   parameter$delta = rep(0,k_delta)
#   parameter$sigma2e = 1
#   parameter$rho = 0

#   parameter$e =lapply(k_e, function(z) rnorm(z,sd =1))
#   parameter$ystar1 = rep(0, sum(n2))
#   parameter$ystar2 = rep(0, sum(n2))
#   parameter
# })

# getLastParameter_single = function(x){
#   parameter=list()
#   parameter$phi = as.vector(tail(x$phi,1))
#   parameter$lambda = as.vector(tail(x$lambda,1))
#   parameter$delta = as.vector(tail(x$delta,1))
#   parameter$sigma2e = as.vector(tail(x$sigma2e,1))
#   parameter$rho = as.vector(tail(x$rho,1))

#   parameter$e = x$e
#   parameter$ystar1 = x$ystar1
#   parameter$ystar2 = x$ystar2

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


merge.single_endogenous = function(x,y,...){
  if (length(list(...))>0){
    list_args = list(...)
    out = list_args[[1]]
    for (i in 2:length(list_args)){
      out = merge(out, list_args[[i]])
    }
    return(out)
  }
  if (x$ID != y$ID){
    stop("Two object have different ID!")
  }
  if (y$index < x$index){
    if (y$index+1 != x$index){
      stop("Indices are not consecutive")
    }
    temp = x
    x = y
    y = temp
    rm(temp)
  }
  out = y
  out$phi = rbind(x$phi,y$phi)
  out$lambda = rbind(x$lambda,y$lambda)
  out$delta = rbind(x$delta,y$delta)
  # out$sigma2e = rbind(x$sigma2e,y$sigma2e)
  out$rho = rbind(x$rho,y$rho)
  out
}

plotmcmc.single_endogenous = function(x, tail=-0.2){
  out = cbind(x$phi,x$rho,x$lambda,x$delta,x$sigma2e)
  plotmcmc.default(out, name=x$name, tail=tail )
}


# createVariableName_single = function(data,data_matrix){
#   c("phi_" %+% colnames(data_matrix$X)
#     ,c("phi_rho")
#     , paste0("lambda")
#     ,"delta_" %+%  c(colnames(data[[1]]$self_data_matrix),  "e_diff")
#     ,"sigma2e"
#   )
# }


getParameterMatrix.single_endogenous = function(x, tail=0){
  # out = cbind(x$phi,x$rho,x$lambda,x$delta,x$sigma2e)
  out = cbind(x$phi,x$rho,x$lambda,x$delta)
  # colnames(out) = x$name

  if (!missing(tail) & tail!=0) {
    out = tail(out,tail)
  }
  out
}

plotmcmc2.single_endogenous = function(x, name, file="mcmc_plot.pdf", tail=10000){
  xx = lapply(x, getParameterMatrix.single_endogenous)
  plotmcmc2.default(xx, name=name, tail=tail )
}

# code <- '
#   Rcpp::List location_index_all_list(location_index_all);
#   Rcpp::List location_index1_list(location_index1);
#   Rcpp::List location_index2_list(location_index2);
#   arma::vec ee = Rcpp::as<arma::vec>(e);
#   arma::vec new_ee = Rcpp::as<arma::vec>(new_e);
#   arma::vec xb11 = Rcpp::as<arma::vec>(xb1);
#   arma::vec xb22 = Rcpp::as<arma::vec>(xb2);
#   arma::vec deltaa = Rcpp::as<arma::vec>(delta);
#   arma::vec prob_old11 = Rcpp::as<arma::vec>(prob_old1);
#   arma::vec prob_old22 = Rcpp::as<arma::vec>(prob_old2);
#   arma::vec ystar11 = Rcpp::as<arma::vec>(ystar1);
#   arma::vec ystar22 = Rcpp::as<arma::vec>(ystar2);
#   arma::vec e_ii = Rcpp::as<arma::vec>(e_i);
#   arma::vec e_jj = Rcpp::as<arma::vec>(e_j);
#   arma::vec uniform_rvv = Rcpp::as<arma::vec>(uniform_rv);
#   arma::vec prob_difff = Rcpp::as<arma::vec>(prob_diff);
#   double sigma22 = Rcpp::as<double>(sigma2);

#   int n = location_index_all_list.size();


#   int update_count=0;
#   for (int i=0; i <n ; i++){
#     SEXP l1 = location_index_all_list[i]; 
#     SEXP l2 = location_index1_list[i]; 
#     SEXP l3 = location_index2_list[i]; 

#     arma::uvec location_index_all_vec = Rcpp::as<arma::uvec>(l1)-1;
#     arma::uvec location_index1_vec = Rcpp::as<arma::uvec>(l2)-1;
#     arma::uvec location_index2_vec= Rcpp::as<arma::uvec>(l3)-1;

#     arma::vec e_i_new = e_ii;
#     arma::vec e_j_new = e_jj;

#     e_i_new(location_index1_vec).fill(new_ee(i)) ;
#     e_j_new(location_index2_vec).fill(new_ee(i)) ; 

#     arma::vec diff_e_new = square((e_i_new(location_index_all_vec)- e_j_new(location_index_all_vec)));



#     arma::vec m1 = ystar11(location_index_all_vec) - xb11(location_index_all_vec) - diff_e_new * deltaa;
#     arma::vec m2 = ystar22(location_index_all_vec) - xb22(location_index_all_vec) - diff_e_new * deltaa;

#     // compute normal density 
#     arma::vec prob_new1 = log( 1/sqrt(2*M_PI*sigma22) * exp(-0.5* square(m1)/sigma22 ) ) ; 
#     arma::vec prob_new2 = log( 1/sqrt(2*M_PI*sigma22) * exp(-0.5* square(m2)/sigma22 ) ) ; 


#     arma::vec prob_old11_subset = prob_old11(location_index_all_vec);
#     arma::vec prob_old22_subset = prob_old22(location_index_all_vec);

#     double a1 = std::accumulate(prob_new1.begin(),prob_new1.end(), 0.0);
#     double a2 = std::accumulate(prob_new2.begin(),prob_new2.end(), 0.0);
#     double a3 = std::accumulate(prob_old11_subset.begin(),prob_old11_subset.end(), 0.0);
#     double a4 = std::accumulate(prob_old22_subset.begin(),prob_old22_subset.end(), 0.0);


#     double alphaa = exp( a1+a2-a3-a4 + prob_difff(i)); 


#     if (alphaa > arma::as_scalar(uniform_rvv(i)) ){
#       e_ii = e_i_new;
#       e_jj = e_j_new;
#       prob_old11(location_index_all_vec) = prob_new1 ;
#       prob_old22(location_index_all_vec) = prob_new2;
#       ee(i) = new_ee(i);
#       update_count++;
#     }
#   }

#   return Rcpp::List::create(Rcpp::Named("e")=ee, Rcpp::Named("update_rate")=update_count/(n+0.0) );

# '

# update_e_internal_single <- cxxfunction(signature(
#   location_index_all="list",
#   location_index1="list",
#   location_index2="list",
#   e="numeric",
#   new_e="numeric",
#   xb1="numeric",
#   xb2="numeric",
#   delta="numeric",
#   prob_old1="numeric",
#   prob_old2="numeric",
#   ystar1="numeric",
#   ystar2="numeric",
#   e_i="numeric",
#   e_j="numeric",
#   uniform_rv="numeric",
#   prob_diff="numeric",
#   sigma2="numeric"
#   ),
#   body=code, 
#   plugin="RcppArmadillo")



# ff <- cxxfunction(signature(
#   x="numeric",
#   sigma2=1
#   ),
#   body=
# '
#   arma::vec xx = Rcpp::as<arma::vec>(x);
#   double sigma22 = Rcpp::as<double>(sigma2);


#   arma::vec prob_new2 = log( 1/sqrt(2*M_PI*sigma22) * exp(-0.5* square(xx)/arma::as_scalar(sigma22) ) ) ; 

#   return wrap(prob_new2);

# '
#   , 
#   plugin="RcppArmadillo")
# ff(1,sigma2=2)




# location_index_all_list = lapply(network_data,"[[","location_index_all")
# location_index1_list = lapply(network_data,"[[","location_index1")
# location_index2_list = lapply(network_data,"[[","location_index2")




#   e_i = e[data$group_index[,1]]
#   e_j = e[data$group_index[,2]]
#   diff_e = (e_i-e_j)^2
#   prob_old1 = dnorm(ystar1 , mean= xb1 + diff_e * tail(delta,1), log=TRUE) 

#   prob_old2 = dnorm(ystar2 , mean= xb2 + diff_e * tail(delta,1), log=TRUE) 
#   prob_diff = dnorm(new_e,log=log) - dnorm(e,log=log) + loglikelihood_outcome_vector_single(e, new_e, lambda, phi, Y, X, WY, rho)

# var(new_e)

# microbenchmark(
# qq=
# update_e_internal_single(
#   location_index_all=location_index_all_list,
#   location_index1=location_index1_list,
#   location_index2=location_index2_list,
#   e=e,
#   new_e=new_e,
#   xb1=xb1,
#   xb2=xb2,
#   delta=tail(delta,1),
#   prob_old1=prob_old1,
#   prob_old2=prob_old2,
#   ystar1=ystar1,
#   ystar2=ystar2,
#   e_i=e_i,
#   e_j=e_j,
#   uniform_rv=uniform_rv,
#   prob_diff=prob_diff
# )
# )


# # qq

# qe=e
# qe[update_count] = new_e[update_count]
# sum(update_count)

# which(qq$e- qe !=0)
# qq[[2]]
# sum(update_count)


single_endogenous_mcmc = function(m, data, last_out){

  network_data = genNetworkData(data)
  data_matrix = genDataMatrix(data)


  n = sapply(data, "[[", "n")
  n2 = sapply(data, function(z) length(z$response_self ))
  G = length(data)

  k_phi = ncol( data_matrix$X )
  k_delta = ( ncol(data[[1]]$self_data_matrix) + 1) 
  k_e = data_matrix$n_vector

  network_name = data[[1]]$network_name

  phi_matrix = matrix(0, nrow=m+1, ncol=k_phi)
  colnames(phi_matrix) = "phi_" %+% colnames(data_matrix$X)
  lambda_matrix = matrix(0, nrow=m+1, ncol=1)
  colnames(lambda_matrix) = "lambda_" %+% network_name
  delta_matrix = matrix(-0.1, nrow=m+1, ncol=k_delta)
  colnames(delta_matrix) = "delta_" %+%network_name %+% "_" %+% c(colnames(data[[1]]$self_data_matrix), "diff_e") 
  # sigma2e_matrix = matrix(1, nrow=m+1, ncol=1)
  rho_matrix = matrix(0, nrow=m+1, ncol=1)
  colnames(rho_matrix) = "rho_" %+% network_name

  switch_sign = 0

  if (!missing(last_out)){
    phi_matrix[1,] = as.vector(tail(last_out$phi,1))
    lambda_matrix[1,] = as.vector(tail(last_out$lambda,1))
    delta_matrix[1,] = as.vector(tail(last_out$delta,1))
    # sigma2e_matrix[1,] = as.vector(tail(last_out$sigma2e,1))
    rho_matrix[1,] = as.vector(tail(last_out$rho,1))
    e= last_out$e
    ystar1 = last_out$ystar1
    ystar2 = last_out$ystar2
    tau=updateTau(last_out$tau, last_out$update_rate)
    index = last_out$index+1
    name = last_out$name
    ID = last_out$ID
  } else {
    e = lapply(data_matrix$n_vector, rep , x=0) 
    ystar1 = rep(0,sum(n2))
    ystar2 = rep(0,sum(n2))

    tau = list(lambda=0.1, rho=0.1)
    tau["e_" %+% network_name] = list(rep(0.5,G))

    index = 1
    ID = genUniqueID()
    cat("Start new instance with ID ", ID, "\n")
    # name =   
    #   c("phi_" %+% colnames(data_matrix$X)
    #   ,"phi_rho", "lambda"
    #   ,"delta_" %+%  c(colnames(data[[1]]$self_data_matrix),  "e_diff")
    #   ,"sigma2e")
  }

  update_rate = list(lambda=0, rho=0)
  update_rate["e_" %+% network_name] = list(rep(0,5))

  self_data_matrix = do.call(rbind, lapply(data, "[[", i="self_data_matrix" ))
  friends_data_matrix = do.call(rbind, lapply(data , "[[", i="friends_data_matrix"))

  y = do.call(rbind,lapply(data, function(z) z$response_self))
  y_not = !y


  tic()
  for (i in 1:m ){
    # if (i %% 1000 == 0 ){
    #   cat(i,">\n")
    # }

    # update e by network
    ystar1_list = splitBy(ystar1,n2)
    ystar2_list = splitBy(ystar2,n2)


    metro_e =  
      update_e_single_new_parser( e, delta=delta_matrix[i,], network_data, tau=tau[["e_"%+%network_name]], lambda=lambda_matrix[i,], phi=phi_matrix[i,], G=G,ystar1_list=ystar1_list, ystar2_list=ystar2_list, rho=rho_matrix[i,])

    e = lapply(metro_e, "[[", i = "e")
    update_rate[["e_"%+%network_name]] = update_rate[["e_"%+%network_name]] + sapply(metro_e, "[[", i="update_rate") 

    update_sign_obj = update_sign_e(lambda_matrix[i,],phi_matrix[i,],rho_matrix[i,], data_matrix,e)
    e = update_sign_obj$e
    switch_sign = switch_sign + update_sign_obj$update


    # sigma2e_matrix[i+1,] = var(unlist(e))



    demean_e = scale(unlist(e),scale=FALSE)
    demean_e_list = splitBy(demean_e,n)
    diff_e = transform_e(data, demean_e_list)

    ## update phi
    phi_matrix[i+1,] = drawingPhi(lambda=lambda_matrix[i,], rho_e=demean_e*rho_matrix[i,],  data_matrix=data_matrix)

    ## update rho
    metro_obj = metropolis(
      beta_previous=rho_matrix[i,], 
      tau=tau$rho, 
      likelihoodFunction=lik_single_exogenous_parser_rho, 
      prior_type="trunc_normal", 
      sampling_type="trunc_normal", 
      phi = phi_matrix[i+1,],
      lambda=lambda_matrix[i,],
      data_matrix = data_matrix,
      e=demean_e
      )
    rho_matrix[i+1,] = metro_obj$beta
    update_rate$rho = update_rate$rho + metro_obj$update

    ## update lambda
    metro_obj = metropolis(
      beta_previous=lambda_matrix[i,], 
      tau=tau$lambda, 
      likelihoodFunction=lik_multi_exogenous, 
      prior_type="normal", 
      sampling_type="normal", 
      phi = phi_matrix[i+1,],
      W_list = data_matrix$W_list,
      X = data_matrix$X,
      Y = data_matrix$Y-demean_e*rho_matrix[i+1,],
      WY = data_matrix$WY,
      n_vector = data_matrix$n_vector
      )
    lambda_matrix[i+1,] = metro_obj$beta
    update_rate$lambda = update_rate$lambda + metro_obj$update



    # q = vector("list",5)
    # for (i in 1:5){
    #   q[[i]] = (e[[i]][data[[i]]$group_index[,1]] -  e[[i]][data[[i]]$group_index[,2]])^2
    # }
    # q = unlist(q)
    # all.equal(q,diff_e,check.attributes = FALSE)

    x1 = cbind(self_data_matrix, diff_e)
    x2 = cbind(friends_data_matrix, diff_e)

    ## update ystar 
    xb1 = x1  %*% delta_matrix[i, ]
    xb2 = x2  %*% delta_matrix[i, ]

    ystar1 = drawYstar(y=y , ystar_other=ystar2, mean=xb1, y_not= y_not)
    ystar2 = drawYstar(y=y, ystar_other=ystar1, mean=xb2, y_not= y_not)

    # update delta
    lm_fit = myFastLm(rbind(x1,x2),c(ystar1, ystar2))

    delta_matrix[i+1, ] =   mvrnorm(n=1, mu=lm_fit$coef, lm_fit$cov/lm_fit$s^2)

    ## convert ystar1 ystar2 to ystar1_list and ystar2_list
    # stopifnot(all.equal(unlist(ystar1_list),ystar1))
    # stopifnot(all.equal(unlist(ystar2_list),ystar2))


  }
  toc()

  out=list()
  out$phi=head(phi_matrix,-1)  
  out$lambda=head(lambda_matrix,-1)
  out$delta=head(delta_matrix,-1)
  out$e=e
  out$ystar1=ystar1
  out$ystar2=ystar2
  out$update_rate=lapply(update_rate, function(z) z/m)
  out$tau=tau
  # out$sigma2e=head(sigma2e_matrix,-1)
  out$rho=head(rho_matrix,-1)
  out$index=index
  out$ID=ID
  # out$name=name
  out$switch_sign=switch_sign/m

   print(out$update_rate)

  class(out) = "single_endogenous"
  out
}



computeLikelihood=function(data, out, network_id,type=0){
  data_matrix = genDataMatrix(data,network_id)


  phi = as.vector( tail(out$phi,1) )
  lambda = as.vector( tail(out$lambda,1) )
  rho = as.vector( tail(out$rho,1) )
  delta = as.vector( tail(out$delta,1) )
  e = out$e

  n = sapply(data,"[[","n")

  demean_e = scale(unlist(e))
  demean_e_list = splitBy(demean_e,n)
  diff_e = transform_e(data, demean_e_list)


  lik1 = lik_single_exogenous_parser_rho(lambda, phi, rho, data_matrix,  demean_e)


  self_data_matrix = do.call(rbind, lapply(data, "[[", i="self_data_matrix" ))
  friends_data_matrix = do.call(rbind, lapply(data , "[[", i="friends_data_matrix"))

  y = unlist(lapply(data, function(z) z$response[[network_id]]))

  p1 = pnorm( cbind(self_data_matrix, diff_e) %*% delta )
  p2 = pnorm( cbind(friends_data_matrix, diff_e) %*% delta )

  p = p1 * p2

  lik2 = (y)*log(p) + (1-y) * log(1-p)

  if (type==1)
    return(sum(lik1))
  if (type==2)
    return(sum(lik2))
  if (type==0)
    return(sum(lik1,lik2)
  )
}

# computeLikelihood(data,chain_list[[1]], 2)
