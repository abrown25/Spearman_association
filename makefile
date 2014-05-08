spearman : spearman.d arg_parse.d calculation.d run_analysis.d
	dmd -O -release -noboundscheck -inline -L-lgsl -L-lgslcblas calculation.d run_analysis.d arg_parse.d spearman.d

.PHONY : perm.p.calc perm tabix.perm fwer clean release gdc

clean:
	rm -f *.o spearman spearman_gdc

release : spearman.d arg_parse.d calculation.d run_analysis.d
	dmd -O -release -noboundscheck -L-lgsl -L-lgslcblas spearman.d arg_parse.d calculation.d run_analysis.d

gdc : spearman.d arg_parse.d calculation.d run_analysis.d
	gdc spearman.d arg_parse.d calculation.d run_analysis.d -L/usr/include/gsl/ -l gsl -L/usr/include/gsl/gsl_cblas -l gslcblas -O3 -o spearman_gdc

perm.p.calc:
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 1000000,12 -pval

perm :
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 4,12

tabix.perm:
	tabix genotype.txt.gz chr1:1-2000 | ./spearman -p phenotype.txt -pi -pc 3 -gs 2 -perm 4,12

fwer:
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 100000,12 -fwer

