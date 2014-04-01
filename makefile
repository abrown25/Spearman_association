spearman : spearman.d arg_parse.d calculation.d
	dmd -L-lgsl -L-lgslcblas spearman.d arg_parse.d calculation.d
