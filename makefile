spearman : main.d arg_parse.d calculation.d run_analysis.d regress.c setup_all.d
	gcc -c regress.c -o regress.o
	dmd -O -release -noboundscheck -inline -L-lgsl -L-lgslcblas calculation.d run_analysis.d arg_parse.d main.d setup_all.d regress.o -ofspearman

gdc : main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.c
	gcc -c regress.c -o regress.o
	gdc main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.o -lgsl -lgslcblas -O3 -o spearman_gdc

ldc : main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.c
	gcc -c regress.c -o regress.o
	ldc2 main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.o -L-lgsl -L-lgslcblas -O3 -of="spearman_ldc"
	rm *.o

unittest : main.d arg_parse.d calculation.d setup_all.d run_analysis.d regress.c
	gcc -c regress.c -o regress.o
	dmd -unittest -L-lgsl -L-lgslcblas calculation.d setup_all.d run_analysis.d arg_parse.d main.d regress.o -ofunittest
	./unittest
	rm -f unittest

.PHONY : perm.p.calc perm tabix.perm fwer clean release gdc time cov all.tests fdr

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

fdr :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -fdr

cov :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -pval -cov cov.txt

tabix.perm :
	tabix genotype.txt.gz chr1:1-2000 | ./spearman -p phenotype.txt -pid -pc 3 -gs 2 -perm 4,12

time :
	time -f "Real : %E, User : %U, System : %S\n" ./spearman_gdc -g large_genotype.txt -perm 1000,4 -o out -fwer large_phenotype.txt

all.tests :
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -o output/simple
	diff output/simple.txt output/simple
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 1000000,12 -pval -o output/perm.p.calc
	diff output/perm.p.calc.txt output/perm.p.calc
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 4,12 -o output/perm
	diff output/perm.txt output/perm
	tabix genotype.txt.gz chr1:1-2000 | ./spearman -p phenotype.txt -pid -pc 3 -gs 2 -perm 4,12 -o output/tabix
	diff output/tabix.txt output/tabix
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -fwer -o output/fwer
	diff output/fwer.txt output/fwer
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -fdr -o output/fdr
	diff output/fdr.txt output/fdr
	cat genotype.txt | ./spearman -p phenotype.txt -pid -gid -pc 3 -gs 2 -perm 100000,12 -pval -cov cov.txt -o output/cov
	diff output/cov.txt output/cov
	rm output/cov output/fwer output/tabix output/perm output/perm.p.calc output/simple output/fdr
