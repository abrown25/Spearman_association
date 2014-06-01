spearman : main.d arg_parse.d calculation.d run_analysis.d regress.o setup_all.d
	gcc -c regress.c -o regress.o -lgsl -lgslcblas
	dmd -O -release -noboundscheck -inline -L-lgsl -L-lgslcblas calculation.d run_analysis.d arg_parse.d main.d setup_all.d regress.o -ofspearman

gdc : main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.o
	gcc -c regress.c -o regress.o -lgsl -lgslcblas
	gdc main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.o -L/usr/include/gsl/ -l gsl -L/usr/include/gsl/gsl_cblas -l gslcblas -O3 -o spearman_gdc

unittest : main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.o
	gcc -c regress.c -o regress.o -lgsl -lgslcblas
	dmd -unittest -L-lgsl -L-lgslcblas calculation.d setup_all.d run_analysis.d arg_parse.d main.d regress.o -ofspearman

.PHONY : perm.p.calc perm tabix.perm fwer clean release gdc time cov all.tests

clean :
	rm -f *.o spearman spearman_gdc

simple :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2

perm :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 4,12

perm.p.calc :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 1000000,12 -pval

fwer :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -fwer

cov :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -pval -cov cov.txt

tabix.perm :
	tabix genotype.txt.gz chr1:1-2000 | ./spearman -p phenotype.txt -pid -pc 3 -gs 2 -perm 4,12

time :
	time -f "Real : %E, User : %U, System : %S\n" ./spearman_gdc -g large_genotype.txt -perm 1000,4 -o out -fwer large_phenotype.txt

all.tests :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 > output/simple
	diff output/simple.txt output/simple
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 1000000,12 -pval > output/perm.p.calc
	diff output/perm.p.calc.txt output/perm.p.calc
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 4,12 > output/perm
	diff output/perm.txt output/perm
	tabix genotype.txt.gz chr1:1-2000 | ./spearman -p phenotype.txt -pid -pc 3 -gs 2 -perm 4,12 > output/tabix
	diff output/tabix.txt output/tabix
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -fwer > output/fwer
	diff output/fwer.txt output/fwer
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -pval -cov cov.txt > output/cov
	diff output/cov.txt output/cov
	rm output/cov output/fwer output/tabix output/perm output/perm.p.calc output/simple
