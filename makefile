spearman : spearman.d arg_parse.d calculation.d
	dmd -L-lgsl -L-lgslcblas spearman.d arg_parse.d calculation.d
.PHONY : perm.p.calc perm
perm.p.calc:
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 1000000 -pval
perm :
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 4,12


