spearman : spearman.d arg_parse.d calculation.d run_analysis.d regress.o
	gcc -c regress.c -o regress.o -lgsl -lgslcblas
	dmd -O -release -noboundscheck -inline -L-lgsl -L-lgslcblas calculation.d run_analysis.d arg_parse.d spearman.d regress.o -ofspearman

.PHONY : perm.p.calc perm tabix.perm fwer clean release gdc time cov

clean:
	rm -f *.o spearman spearman_gdc

gdc : spearman.d arg_parse.d calculation.d run_analysis.d regress.o
	gcc -c regress.c -o regress.o -lgsl -lgslcblas
	gdc spearman.d arg_parse.d calculation.d run_analysis.d regress.o -L/usr/include/gsl/ -l gsl -L/usr/include/gsl/gsl_cblas -l gslcblas -O3 -o spearman_gdc

simple:
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2

perm.p.calc:
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 1000000,12 -pval

perm :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 4,12

tabix.perm:
	tabix genotype.txt.gz chr1:1-2000 | ./spearman -p phenotype.txt -pid -pc 3 -gs 2 -perm 4,12

fwer:
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -fwer

time:
	time -f "Real : %E, User : %U, System : %S\n" ./spearman_gdc -g large_genotype.txt -perm 100,4 -o out -fwer large_phenotype.txt

cov:
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -pval -cov cov.txt
