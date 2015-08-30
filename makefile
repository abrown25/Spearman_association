#Replace with location of static GSL libbrary if you want to make statically linked version
GSL = /usr/lib/libgsl.a /usr/lib/libgslcblas.a
DSOURCES = src/main.d src/arg_parse.d src/calculation.d src/setup_all.d src/run_analysis.d

gdc : ${DSOURCES} src/regress.o
	gdc -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_np_gwas ${DSOURCES} src/regress.o -lgsl -lgslcblas -O3 -o bin/np_gwas
	strip bin/np_gwas
	rm -f src/*.o

static : ${DSOURCES} src/regress.o
	gdc -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_np_gwas ${DSOURCES} src/regress.o ${GSL} -O3 -o bin/np_gwas
	strip bin/np_gwas
	rm -f src/*.o

test : ${DSOURCES} src/regress.o
	gdc -fdebug -g -funittest -Werror -Wall ${DSOURCES} src/regress.o -lgsl -lgslcblas -o unittest
	./unittest
	rm -f unittest src/*.o *.o

ldc : src/main.d src/arg_parse.d src/calculation.d src/run_analysis.d src/regress.o src/setup_all.d
	ldc -release -enable-inlining -O -w -oq ${DSOURCES} src/regress.o -L-lgsl -L-lgslcblas -of="bin/np_gwas"
	rm -f *.o
	rm -f src/*.o

ldc_test : ${DSOURCES} src/regress.o
	ldc -d-debug -g -unittest -w -L-lgsl -L-lgslcblas ${DSOURCES} src/regress.o -of="unittest"
	./unittest
	rm -f unittest src/*.o *.o

dmd : ${DSOURCES} src/regress.o
	dmd -O -release -noboundscheck -inline -L-lgsl -L-lgslcblas ${DSOURCES} src/regress.o -ofbin/np_gwas
	rm src/*.o

dmd_test : ${DSOURCES} src/regress.o
	dmd -debug -g -unittest -w -L-lgsl -L-lgslcblas ${DSOURCES} src/regress.o -ofunittest
	./unittest
	rm -f unittest src/*.o *.o

genotype_utilities : src/genotype_utilities.d
	gdc -frelease -finline-functions -O3 -Werror -Wall src/genotype_utilities.d -o bin/genotype_utilities

.PHONY : test static ldc dmd ldc_test dmd_test clean install

clean :
	rm -f src/*.o bin/np_gwas

install : bin/largeQvalue largeQvalue.1
	ln -s $(shell pwd)/bin/np_gwas /usr/local/bin/
	ln -s $(shell pwd)/np_gwas.1 /usr/local/man/man1/
