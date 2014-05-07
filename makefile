spearman : spearman.d arg_parse.d calculation.d run_analysis.d
	dmd -L-lgsl -L-lgslcblas spearman.d arg_parse.d calculation.d run_analysis.d

.PHONY : perm.p.calc perm tabix.perm fwer clean release

clean:
	rm -f spearman.o spearman

release : spearman.d arg_parse.d calculation.d run_analysis.d
	dmd -O -release -noboundscheck -L-lgsl -L-lgslcblas spearman.d arg_parse.d calculation.d run_analysis.d

perm.p.calc:
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 1000000,12 -pval

perm :
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 4,12

tabix.perm:
	tabix genotype.txt.gz chr1:1-2000 | ./spearman -p phenotype.txt -pi -pc 3 -gs 2 -perm 4,12

fwer:
	cat genotype.txt | ./spearman -p phenotype.txt -pi -gi -pc 3 -gs 2 -perm 100000,12 -fwer

