GSL = /usr/lib/libgsl.so /usr/lib/libgslcblas.so

np_gwas : src/main.d src/arg_parse.d src/calculation.d src/run_analysis.d src/regress.c src/setup_all.d
	gcc -c src/regress.c -o regress.o
	ldc -release -enable-inlining -O -w -oq src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d regress.o -L-lgsl -L-lgslcblas -of="np_gwas"
	rm *.o

gdc : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o regress.o
	gdc src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d regress.o ${GSL} -O3 -o np_gwas_gdc

dmd : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o regress.o
	dmd -O -release -noboundscheck -inline ${GSL} src/calculation.d src/run_analysis.d src/arg_parse.d src/main.d src/setup_all.d regress.o -ofnp_gwas_dmd

unittest : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o regress.o
	ldc -d-debug -g -unittest -w -L-lgsl -L-lgslcblas src/calculation.d src/setup_all.d src/run_analysis.d src/arg_parse.d src/main.d regress.o -of="unittest"
	./unittest
	rm -f unittest *.o

genotype_utilities : src/genotype_utilities.d
	ldc -release -enable-inlining -O -w -oq src/genotype_utilities.d -of="genotype_utilities"
	rm -f *.o

.PHONY : perm.p.calc perm tabix.perm fwer clean release gdc time cov all.tests fdr

clean :
	rm -f *.o np_gwas np_gwas_gdc np_gwas_dmd

simple :
	cat data/genotype.txt | np_gwas -pdata/phenotype.txt --pid --gid --pc 3 --gs 2

perm :
	cat data/genotype.txt | np_gwas --pheno data/phenotype.txt --pid --gid --pc=3 --gs 2 --perm 4,12

perm.p.calc :
	cat data/genotype.txt | np_gwas -pdata/phenotype.txt --pid --gid --pc 3 --gs 2 --perm 1000000,12 --pval

fwer :
	cat data/genotype.txt | np_gwas --p data/phenotype.txt --pid --gid --pc 3 --gs 2 --perm 100000,12 --fwer

fdr :
	cat data/genotype.txt | np_gwas --p data/phenotype.txt --pid --gid --pc 5 --gs 2 --perm 100000,12 --fdr

cov :
	cat data/genotype.txt | np_gwas -p data/phenotype.txt --pid --gid nonsense --pc 3 --gs 2 --perm 100000,12 --pval --cov cov.txt

tabix.perm :
	tabix genotype.txt.gz chr1:1-2000 | np_gwas --p data/phenotype.txt --pid --pc 3 --gs 2 --perm 4,12

time :
	time -f "Real : %E, User : %U, System : %S\n" np_gwas -gdata/large_genotype.txt --perm 1000,4 -oout --fwer data/large_phenotype.txt

all.tests :
	cat data/genotype.txt | np_gwas --p data/phenotype.txt --pid --gid --pc 3 --gs 2 -ooutput/simple
	diff output/simple.txt output/simple
	cat data/genotype.txt | np_gwas -pdata/phenotype.txt --pid --gid --pc 3 --gs 2 --perm 1000000,12 --pval --o output/perm.p.calc
	diff output/perm.p.calc.txt output/perm.p.calc
	cat data/genotype.txt | np_gwas --p data/phenotype.txt --pid --gid --pc 3 --gs 2 --perm 4,12 --o output/perm
	diff output/perm.txt output/perm
	tabix data/genotype.txt.gz chr1:1-2000 | np_gwas --p data/phenotype.txt --pid --pc 3 --gs 2 --perm 4,12 --o output/tabix
	diff output/tabix.txt output/tabix
	cat data/genotype.txt | np_gwas --p data/phenotype.txt --pid --gid --pc 3 --gs 2 --perm 100000,12 --fwer --o output/fwer
	diff output/fwer.txt output/fwer
	cat data/genotype.txt | np_gwas --p data/phenotype.txt --pid --gid --pc 5 --gs 2 --perm 100000,12 --fdr --o output/fdr
	diff output/fdr.txt output/fdr
	cat data/genotype.txt | np_gwas --p data/phenotype.txt --pid --gid --pc 3 --gs 2 --perm 100000,12 --pval --cov cov.txt --o output/cov
	diff output/cov.txt output/cov
	rm output/cov output/fwer output/tabix output/perm output/perm.p.calc output/simple output/fdr
