## Given a seq, beta, D, compute the likelihood
## Given a seq m  n*(n-1)/2 x 2 matrix
## D: n by n network matrix
## U_xb : an n by n utility matrix, i,j element is the utility of i to make friends with j. (Xbeta)
## delta1 delta2 


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

# computeNetworkSummary = function(seq_m,D){
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
# 		q1= f(seq_m=seq_m, D=D)
#   }
#   ,b={
# 		q2= g(seq_m=seq_m, D=D)
#   }
# )

# all(q1[[1]]==q2[[1]])
# all(q1[[2]]==q2[[2]])


computeNetworkSummary2=function(seq_m=seq_m_new, D=D){
	# out = list()
	# for (i in 1:length(D)){
	# 	temp = computeNetworkSummary(g(seq_m=seq_m[[i]], D=D[[i]]))
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

# q1 = computeNetworkSummary2(seq_m,D)



#############################################




SNF_single_mcmc = function(m, data,  last_out, update_tau=TRUE,tau=0.0005){

  # if (network_id==1){
		# D = lapply(data, function(z) z$W!=0 )
		# y = do.call(c, lapply(data, "[[", "response1"))
		# y_not = !y
  # } else{
		# D = lapply(data, function(z) z$W2!=0 )
		# y = do.call(c, lapply(data, "[[", "response2"))
		# y_not = !y
  # }

  D = lapply(data, function(z) z$D_list[[1]])
  y = do.call(c, lapply(data,"[[","response_self") )
  y_not = !y

	n= sapply(data,"[[","n")
	n2 = sapply(data, function(z) length(z$response_self))

	x1 = do.call(rbind, lapply(data, "[[" , "self_data_matrix") )
	x2 = do.call(rbind, lapply(data, "[[" , "friends_data_matrix") )
	
	number_of_network_variable = 3

	postition = mapply(seq, c(0,head(cumsum(n2),-1)) +1,cumsum(n2))
	ystar1 = rep(0,length(y))
	ystar2 = rep(0,length(y))
	seq_m = lapply(data,"[[","group_index")

	delta_matrix = matrix(0, nrow=m+1, ncol= ncol(x1) + number_of_network_variable )
	update_rate = 0

	## initialization 

	if (!missing(last_out) && !is.null(last_out) ){
		cat("Using last_out \n")
		ystar1 = last_out$ystar1
		ystar2 = last_out$ystar2
		seq_m = last_out$seq_m
		delta_matrix[1,] = as.vector(tail(last_out$delta,1))

    index = last_out$index+1
    # name = last_out$name
    ID = last_out$ID

		if (update_tau){
    	tau=updateTau(last_out$tau, last_out$update_rate, lower_bound=0.2, upper_bound=0.4,optim_rate=0.3,min_rate=0.00001)
		} else{
			tau=last_out$tau
		}
	} else {
    index = 1
    ID = genUniqueID()
    cat("Start new instance with ID ", ID, "\n")
  }

	network_summary = computeNetworkSummary2(seq_m=seq_m, D=D)

	xx1 = cbind(x1,network_summary$self)
	xx2 = cbind(x2,network_summary$friends)

	xb1 = xx1 %*% delta_matrix[1,]
	xb2 = xx2 %*% delta_matrix[1,]

	colnames(delta_matrix) = colnames(xx1)
  name = colnames(xx1) 

	delta_x_index = 1:ncol(x1)
	delta_network_index = 1:number_of_network_variable + ncol(x1)

	tic()
	## start the gibbs
	for (i in 1:m){
		## base on the seq, compute the network summary
		## draw ystar

		ystar1 = drawYstar_single(y=y , ystar_other=ystar2, mean=xb1, y_not= y_not)
		ystar2 = drawYstar_single(y=y, ystar_other=ystar1, mean=xb2, y_not= y_not)

		## draw delta
		lm_fit = myFastLm(X= rbind(xx1,xx2), y = c(ystar1,ystar2))
		delta_matrix[i+1, ] =   mvrnorm(n=1, mu=lm_fit$coef, lm_fit$cov/lm_fit$s^2)

    R1 = x1 %*% delta_matrix[i+1, delta_x_index]
    R2 = x2 %*% delta_matrix[i+1, delta_x_index]

    xb1 = R1 + network_summary$self %*% delta_matrix[i+1,delta_network_index]
    xb2 = R2 + network_summary$friends %*% delta_matrix[i+1,delta_network_index]


		## update sequence 
		seq_m_new = DrawSeqSample(seq_m,p=tau)

		# sapply(1:5, function(i) sum(seq_m_new[[i]]!=seq_m[[i]]) )

		network_summary_new = computeNetworkSummary2(seq_m=seq_m_new, D=D)

    xb1_new = R1 + network_summary_new$self %*% delta_matrix[i+1,delta_network_index]
    xb2_new = R2 + network_summary_new$friends %*% delta_matrix[i+1,delta_network_index]


		p1 = splitBy(dnorm(ystar1 - xb1, log=TRUE),by=n2)
		p2 = splitBy(dnorm(ystar2 - xb2, log=TRUE),by=n2)
		p1_new = splitBy(dnorm(ystar1 - xb1_new, log=TRUE),by=n2)
		p2_new = splitBy(dnorm(ystar2 - xb2_new, log=TRUE),by=n2)

		p1 = sapply(p1, sum)
		p2 = sapply(p2, sum)
		p1_new = sapply(p1_new, sum)
		p2_new = sapply(p2_new, sum)

		alpha = exp( p1_new+ p2_new - p1- p2  )
		update_index = alpha  > runif(5)
		seq_m[update_index] = seq_m_new[update_index]
		update_rate = update_rate + update_index


		update_position = unlist(postition[update_index])
    network_summary$self[update_position,] = network_summary_new$self[update_position,]
    network_summary$friends[update_position,] = network_summary_new$friends[update_position,]
	  xb1[update_position] = xb1_new[update_position]
	  xb2[update_position] = xb2_new[update_position]

		xx1[update_position,delta_network_index] = network_summary$self[update_position,]
		xx2[update_position,delta_network_index] = network_summary$friends[update_position,]

		# test
		# xx1_q = cbind(x1,network_summary$self)
		# xx2_q = cbind(x2,network_summary$friends)

		# network_summary_q = computeNetworkSummary2(seq_m=seq_m, D=D)
		# xx1_q = cbind(x1,network_summary_q$self)
		# xx2_q = cbind(x2,network_summary_q$friends)
		# xb1_q = xx1_q %*% delta_matrix[i+1,]
		# xb2_q = xx2_q %*% delta_matrix[i+1,]

		# identical(xx1,xx1_q)
		# identical(xx2,xx2_q)

		# identical(xb1,xb1_q)
		# identical(xb2,xb2_q)


	}
	toc()


	update_rate = update_rate/m
	cat("Update rate : \n")
	print(update_rate)

	out = list(delta=tail(delta_matrix,-1) , seq_m=seq_m,ystar1=ystar1,ystar2=ystar2, tau=tau, update_rate=update_rate, index=index,ID=ID, name=name)

	class(out) = "SNF_single"
	out
}

merge.SNF_single = function(x,y,...){
  out = y 
  out$delta_matrix = rbind(x$delta_matrix, y$delta_matrix)
  out 
}

getParameterMatrix.SNF_single = function(x ){
  x$delta
}






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


# network_summary = computeNetworkSummary2(seq_m=seq_m, D=D)

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
# 	network_summary = computeNetworkSummary2(seq_m=seq_m, D=D)

# 	xx1 = cbind(x1,network_summary$self)
# 	xx2 = cbind(x2,network_summary$friends)
# 	xb1 = xx1 %*% delta_matrix[1,]
# 	xb2 = xx2 %*% delta_matrix[1,]

# 	ystar1 = drawYstar_single(y=y , ystar_other=ystar2, mean=xb1, y_not= y_not)
# 	ystar2 = drawYstar_single(y=y, ystar_other=ystar1, mean=xb2, y_not= y_not)

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
