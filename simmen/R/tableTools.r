# exogenous_model
# exogenous_model2 = exogenous_model
# exogenous_model2[[2]] = exogenous_model2[[2]] [[2]]
# exogenous_model2[[4]] = exogenous_model2[[4]] [[2]]
# exogenous_model2[[6]] = exogenous_model2[[6]] [[2]]

# all_rownames_list = lapply(exogenous_model2, rownames)

# all_rownames= unlist(all_rownames_list)

# names(all_rownames) = NULL

# unique_all_rownames= all_rownames[!duplicated(all_rownames)]


## x is a list of summary table of model
combineTable = function(x, grouping = c("phi","lambda","friends","studymates") ,sort = FALSE , remove=c("std","bfi","I((play","matrix","test")){
	model_name = names(x)
	variable_names_list = lapply(x, rownames)

	variable_names = unlist(variable_names_list)

	names(variable_names) = NULL

	unique_all_rownames= variable_names[!duplicated(variable_names)]

	if (missing(grouping)){
		sort( table( unlist( strsplit(unique_all_rownames, "[._-]") ) ) ,decreasing=T)
		# do something
	}

	new_name_order = seq(length(unique_all_rownames))

	for (i in seq(grouping)){
		pattern = grouping[i]
		pattern_matched = str_detect(unique_all_rownames, pattern)
		new_name_order = c(which( pattern_matched ) , which( !pattern_matched ))
		  
	}

}	

library(stringr)
reordering = function(x ,by){
	pattern = by[1]
	pattern_matched = str_detect(x, pattern)
	if (length(by)==1){
		x[ c(which( pattern_matched ) , which( !pattern_matched )) ]
	} else{
		c(reordering(x[which(pattern_matched)], by[-1]), reordering(x[which(!pattern_matched)], by[-1]) )
	}
}

# reordering(x,c("phi","friends","studymates"))


## need something to check for nested group

