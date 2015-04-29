GSL = /usr/lib/libgsl.so /usr/lib/libgslcblas.so

np_gwas : src/main.d src/arg_parse.d src/calculation.d src/run_analysis.d src/regress.c src/setup_all.d
	gcc -c src/regress.c -o bin/regress.o
	ldc -release -enable-inlining -O -w -oq src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d bin/regress.o -L-lgsl -L-lgslcblas -of="bin/np_gwas"
	rm -f *.o
	rm -f bin/*.o

gdc : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o bin/regress.o
	gdc -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_np_gwas src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d bin/regress.o ${GSL} -O3 -o bin/np_gwas
	strip bin/np_gwas
	rm -f bin/*.o

dmd : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o bin/regress.o
	dmd -O -release -noboundscheck -inline ${GSL} src/calculation.d src/run_analysis.d src/arg_parse.d src/main.d src/setup_all.d bin/regress.o -ofbin/np_gwas_dmd
	rm bin/*.o

unittest : src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d src/regress.c
	gcc -c src/regress.c -o regress.o
	ldc -d-debug -g -unittest -w -L-lgsl -L-lgslcblas src/calculation.d src/setup_all.d src/run_analysis.d src/arg_parse.d src/main.d regress.o -of="unittest"
	./unittest
	rm -f unittest *.o

genotype_utilities : src/genotype_utilities.d
	ldc -release -enable-inlining -O -w -oq src/genotype_utilities.d -of="bin/genotype_utilities"
	rm -f bin/*.o

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
