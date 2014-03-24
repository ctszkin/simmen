# lik_single_exogenous2 = function(lambda, phi, data, network_id=1){
#   if (abs(lambda)>1)
#     return(-1e+20)
#   log_det_Gamma = 0 
#   new_y = NULL

#   if (network_id==1){
#     W = data$W
#     wy = data$wy
#     wx = data$wx
#   } else if (network_id==2){
#     W = data$W2
#     wy = data$wy2
#     wx = data$wx2
#   }

#   log_det_Gamma = determinant(diag(data$n) - lambda * W,log=TRUE)$modulus

#   new_y = data$y - lambda*wy 

#   if (missing(phi)){
#     lm_fit = fastLmPure(x=X, y=new_y, method=2)
#     phi = lm_fit$coefficients
#     varepsilon = lm_fit$residuals
#   } else{
#     varepsilon = new_y - X %*% phi
#   }

#   out = - (data$n-1)/2* log(sum(varepsilon^2)/(data$n-1)) + log_det_Gamma - log(1 - lambda)

#   if (!is.finite(out)){
#     out = -1e+20
#   }
#   out
# })



# lik_single_exogenous = function(lambda, phi, W_list, n_vector, X, Y, WY){
#   if (abs(lambda)>1)
#     return(-1e+20)

#   log_det_Gamma = 0
#   for (i in 1:length(W_list) ){
#     log_det_Gamma = log_det_Gamma + genLogDetGamma(W_list[[i]], lambda=lambda)
#   }

#   # log_det_Gamma = sum(sapply(W_list, genLogDetGamma, lambda=lambda)) 

#   G = length(n_vector)

#   new_y = Y - lambda*WY 
#   if (missing(phi)){
#     lm_fit = myFastLm(XX=X, yy=new_y)
#     phi = lm_fit$coefficients
#     varepsilon = lm_fit$residuals
#   } else {
#     varepsilon = new_y - X %*% phi
#   }


#   splitted_varepsilon = splitBy(varepsilon,  n_vector )
#   loglik_varepsilon = 0
#   for (i in 1:G) {
#     nn = length(splitted_varepsilon[[i]])
#     loglik_varepsilon = loglik_varepsilon -(nn-1)/2*log(sum(splitted_varepsilon[[i]]^2)/(nn-1))
#   }
#   # loglik_varepsilon = sum(
#   #   sapply(splitted_varepsilon, function(z) {
#   #     -(length(z)-1)/2*log(sum(z^2)/(length(z)-1))
#   #   })
#   # )
#   # loglik_varepsilon = -(n-1)/2 * log(sum(varepsilon^2)/(n-1))

#   out = loglik_varepsilon + log_det_Gamma - G *log(1 - lambda)

#   if (!is.finite(out)){
#     out = -1e+20
#   }
#   out
# })



# single_exogenous =function(data, network_id=1){
#   data_matrix = genDataMatrix(data,network_id)

#   lambda = 0 

#   for (i in 1:10){
#     new_y = data_matrix$Y - lambda * data_matrix$WY
#     phi = myFastLm(X=data_matrix$X, y=new_y)$coef

#     lik_obj_initial = maxLik(function(z, ...) lik_single_exogenous(lambda=z, ...), 
#       method="BFGS",
#       start=0, 
#       W_list=data_matrix$W_list, 
#       Y=data_matrix$Y, 
#       WY=data_matrix$WY, 
#       X=data_matrix$X, 
#       n_vector = data_matrix$n_vector,
#       phi=phi
#     )

#     lambda = lik_obj_initial$estimate
#   }

#   start0 = c(phi,lambda=lambda)

#   # yy= unlist(lapply(data, function(z) z$y ))
#   # wyy = unlist(lapply(data, function(z) z$wy ))
#   # summary(lm(yy~X+wyy))
#   # q=
#   out = maxLik(function(z, ...) lik_single_exogenous(phi=head(z,-1), lambda=tail(z,1), ...), 
#     method="BFGS",
#     start=start0, 
#     W_list=data_matrix$W_list, 
#     Y=data_matrix$Y, 
#     WY=data_matrix$WY, 
#     X=data_matrix$X, 
#     n_vector = data_matrix$n_vector

#   )
#   row_names = c(colnames(data_matrix$X) , data[[1]]$network_name[network_id] %+% "_y")
#   summary_table = generateSignificance(summary(out)$estimate[,1:2], row_names=row_names)

#   c(list(out),list(summary_table))
# })

# ols_single = function(data, network_id=1){
#   if (network_id==1){
#     X = do.call(rbind, lapply(data,function(z) cbind(z$x,z$wx)) )
#     Y = unlist(lapply(data, function(z) z$y))
#     WY = unlist(lapply(data, function(z) z$wy ))
#   } else if (network_id==2){
#     X = do.call(rbind, lapply(data,function(z) cbind(z$x,z$wx2)) )
#     Y = unlist(lapply(data, function(z) z$y))
#     WY = unlist(lapply(data, function(z) z$wy2 ))
#   }
#   out = generateSignificance(summary(lm(Y~X+WY-1))$coef[,1:2])
#   row.names(out) = c(colnames(X),data[[1]]$network_name[network_id] %+% "_y")
#   out
# })

# ols_multi = function(data){
#   X = do.call(rbind, lapply(data,function(z) cbind(z$x,z$wx,z$wx2,z$wy,z$wy2)) )
#   Y = unlist(lapply(data, function(z) z$y))

#   out = generateSignificance(summary(lm(Y~X-1))$coef[,1:2])
#   row.names(out) = c(head(colnames(X),-2),data[[1]]$network_name %+% "_y")
#   out
# })

peer_effect_ols = function(data){
  data_matrix = genDataMatrix(data)
  lm_fit = myFastLm(cbind(data_matrix$X,data_matrix$WY), data_matrix$Y)
  out = generateSignificance(cbind(lm_fit$coefficients, sqrt(diag(lm_fit$cov)) ) )
  row.names(out) = c( "phi_" %+% colnames(data_matrix$X) , "lambda_" %+% data[[1]]$network_name )
  out
}

ols = function(data){
  data_matrix = genDataMatrix(data, any_wx=FALSE)
  lm_fit = myFastLm(cbind(data_matrix$X), data_matrix$Y)
  out = generateSignificance(cbind(lm_fit$coefficients, sqrt(diag(lm_fit$cov)) ) )
  row.names(out) = c( "phi_" %+% colnames(data_matrix$X) )
  out
}

peer_effect_ols_exogenous = function(data){
  data_matrix = genDataMatrix(data)
  lm_fit = myFastLm(cbind(data_matrix$X), data_matrix$Y)
  out = generateSignificance(cbind(lm_fit$coefficients, sqrt(diag(lm_fit$cov)) ) )
  row.names(out) = c( "phi_" %+% colnames(data_matrix$X) )
  out
}




genLogDetGamma = function(W_list, lambda){
  out = 0
  for (i in 1:length(W_list)){
    out = out - lambda[i] * W_list[[i]]
  }
  diag(out) = 1
  determinant(out, logarithm=TRUE)$modulus
}


lik_multi_exogenous = function(lambda, phi, X, Y, WY, W_list, n_vector){

  if (any(sum(abs(lambda))>=1)  ) 
    return(-1e+20)

  log_det_Gamma= 0 
  new_y = NULL

  G = length(n_vector)

  for (i in 1:G){
    log_det_Gamma = log_det_Gamma+ genLogDetGamma(W_list[[i]], lambda)
  }

  new_y = Y - WY %*% lambda

  if (missing(phi)){
    lm_fit = myFastLm(X=X, y=new_y)
    phi = lm_fit$coef
    varepsilon = lm_fit$residuals
  } else {
    varepsilon = new_y - X %*% phi
  }

  splitted_varepsilon = split(varepsilon, rep(1:G, n_vector))
  loglik_varepsilon = sum(
    sapply(splitted_varepsilon, function(z) {
      -(length(z)-1)/2*log(sum(z^2)/(length(z)-1))
    })
  )

  out = loglik_varepsilon + log_det_Gamma - G * log(1 - sum(lambda))

  # out = -(n-1)/2 * log(sum(varepsilon^2)/(n-1)) + log_det_Gamma - length(data) * log(1 - sum(lambda))
  if (!is.finite(out)){
    out = -1e+20
  }
  out
}


peer_effect_exogenous =function(data){
  data_matrix = genDataMatrix(data)

  number_of_network = length(data_matrix$W_list[[1]])

  lambda= rep(0, number_of_network )
  for (i in 1:10){

    new_y = data_matrix$Y - data_matrix$WY %*% lambda

    phi = lm.fit(x=data_matrix$X, y=new_y)$coef
    q1=maxLik(function(z, ...) lik_multi_exogenous(lambda=z,...), 
      method="BFGS",
      phi=phi,
      start=lambda, 
      X = data_matrix$X, 
      Y = data_matrix$Y,
      WY = data_matrix$WY,
      W_list = data_matrix$W_list, 
      n_vector = data_matrix$n_vector
    )
    lambda = q1$estimate
  }

  start0 = c(phi,lambda=lambda)

  tic()
  out = maxLik(function(z, ...) lik_multi_exogenous(phi=head(z,-number_of_network), lambda=tail(z,number_of_network), ...) , 
    method="BFGS",
    start=start0, 
    X = data_matrix$X, 
    Y = data_matrix$Y,
    WY = data_matrix$WY,
    W_list = data_matrix$W_list, 
    n_vector = data_matrix$n_vector
  )
  toc()
  
  row_names = c("phi_" %+% colnames(data_matrix$X), "lambda_" %+% data[[1]]$network_name)

  summary_table = generateSignificance(summary(out)$estimate[,1:2], row_names=row_names)

  c(list(out),list(summary_table))
}

