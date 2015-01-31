GSL = /usr/lib/libgsl.so /usr/lib/libgslcblas.so

spearman : src/main.d src/arg_parse.d src/calculation.d src/run_analysis.d src/regress.c src/setup_all.d
	gcc -c src/regress.c -o regress.o
	dmd -O -release -noboundscheck -inline ${GSL} src/calculation.d src/run_analysis.d src/arg_parse.d src/main.d src/setup_all.d regress.o -ofspearman

gdc : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o regress.o
	gdc src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d regress.o ${GSL} -O3 -o spearman_gdc

ldc : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o regress.o
	ldc src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d regress.o -L-lgsl -L-lgslcblas -O3 -of="spearman_ldc"
	rm *.o

unittest : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o regress.o
	dmd -unittest -L-lgsl -L-lgslcblas src/calculation.d src/setup_all.d src/run_analysis.d src/arg_parse.d src/main.d regress.o -ofunittest
	./unittest
	rm -f unittest

.PHONY : perm.p.calc perm tabix.perm fwer clean release gdc time cov all.tests fdr

clean :
	rm -f *.o spearman spearman_gdc spearman_ldc

simple :
	cat genotype.txt | ./spearman -pphenotype.txt --pid --gid --pc 3 --gs 2

perm :
	cat genotype.txt | ./spearman --pheno phenotype.txt --pid --gid --pc=3 --gs 2 --perm 4,12

perm.p.calc :
	cat genotype.txt | ./spearman -pphenotype.txt --pid --gid --pc 3 --gs 2 --perm 1000000,12 --pval

fwer :
	cat genotype.txt | ./spearman --p phenotype.txt --pid --gid --pc 3 --gs 2 --perm 100000,12 --fwer

fdr :
	cat genotype.txt | ./spearman --p phenotype.txt --pid --gid --pc 5 --gs 2 --perm 100000,12 --fdr

cov :
	cat genotype.txt | ./spearman -p phenotype.txt --pid --gid nonsense --pc 3 --gs 2 --perm 100000,12 --pval --cov cov.txt

tabix.perm :
	tabix genotype.txt.gz chr1:1-2000 | ./spearman --p phenotype.txt --pid --pc 3 --gs 2 --perm 4,12

time :
	time -f "Real : %E, User : %U, System : %S\n" ./spearman_gdc -glarge_genotype.txt --perm 1000,4 -oout --fwer large_phenotype.txt

all.tests :
	cat genotype.txt | ./spearman --p phenotype.txt --pid --gid --pc 3 --gs 2 -ooutput/simple
	diff output/simple.txt output/simple
	cat genotype.txt | ./spearman -pphenotype.txt --pid --gid --pc 3 --gs 2 --perm 1000000,12 --pval --o output/perm.p.calc
	diff output/perm.p.calc.txt output/perm.p.calc
	cat genotype.txt | ./spearman --p phenotype.txt --pid --gid --pc 3 --gs 2 --perm 4,12 --o output/perm
	diff output/perm.txt output/perm
	tabix genotype.txt.gz chr1:1-2000 | ./spearman --p phenotype.txt --pid --pc 3 --gs 2 --perm 4,12 --o output/tabix
	diff output/tabix.txt output/tabix
	cat genotype.txt | ./spearman --p phenotype.txt --pid --gid --pc 3 --gs 2 --perm 100000,12 --fwer --o output/fwer
	diff output/fwer.txt output/fwer
	cat genotype.txt | ./spearman --p phenotype.txt --pid --gid --pc 5 --gs 2 --perm 100000,12 --fdr --o output/fdr
	diff output/fdr.txt output/fdr
	cat genotype.txt | ./spearman --p phenotype.txt --pid --gid --pc 3 --gs 2 --perm 100000,12 --pval --cov cov.txt --o output/cov
	diff output/cov.txt output/cov
	rm output/cov output/fwer output/tabix output/perm output/perm.p.calc output/simple output/fdr
