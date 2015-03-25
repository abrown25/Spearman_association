##np_gwas : Fast, non parametric tests for genetic associations with permutations

###Introduction:

The purpose of this program is to allow genome wide analyses of non normally distributed data to be performed quickly and easily. Fast routines to calculate permutations are included, with options to calculate p values adjusted for multiple testing, either the family wise error rate (FWER) or the typically less conservative false discovery rate (FDR). The program is designed to be flexible with the data formats it will accept: genotype data should come with individuals are columns and SNPs in rows, with an optional header line and number of preceding columns giving SNP information. The program has been designed to work with linux pipes, allowing it to be seamlessly chained with programs to manipulate SNP data such as tabix and vcf-tools. Phenotype information should be given in a specified column of a text file: if there are known confounders, these can be passed as a matrix and controlled out using linear regression.
###Usage:

    Usage: no_gwas [options]:

    Options:
    --help    Display help information.
    --version Display version information.
    --pheno, --p CHAR
        Phenotype file [last argument].
    --geno, --g CHAR
        Genotype file [stdin].
    --out, --o CHAR
        Output file [stdout].
    --cov, --c CHAR
        Optional covariates file, if specified analysis will be performed on the residuals,
        after controlling for covariates with least squares regression.
    --ttest   Runs a test of standard parametric correlation between genotype and phenotype.
    --pheno-id, --pid
        Phenotype IDs are in the first column, if genotype IDs are also present then we
        check for mismatches.
    --geno-id, --gid
        Genotype IDs are in the first row, if phenotype IDs are also present then we check
        for mismatches.
    --pheno-col, --pc INT
        Column for phenotype values, default is 1 if phenotype IDs are not present, 2 other-
        wise.
    --geno-skip, --gs INT
        Column at which genotype values start, preceding columns are printed.
    --perm INT(,INT)
        Calculated permuted p values, one following number indicates the number of permuta-
        tions, two comma separated numbers gives the number of permutations and the seed.
    --pval    Report permutation p values for each test (needs perm options to be specified).
    --fwer    Calculates the Family Wise Error Rate (FWER) based on permutations, corrected P val-
        ues in the last column.
    --fdr     Calculates the False Discovery Rate (FDR) based on permutations, corrected P values
        in the last column.
    --pi DOUBLE
        Specify π₀, proportion of null hypotheses, to more accurately estimate False Discov-
        ery Rate.
    --nocheck Skip check of IDs when both genotype and phenotype IDs are present.
    --match   The program will rearrange the phenotype data so that the genotype and phenotype IDs
        match. If individuals are present in the genotype file only, the analysis will be
        halted.

    FILE FORMATS

    INPUT FILE FORMATS
      Phenotype
        Tab or whitespace separated file with phenotype values in column specified by --pc, and
	    optional subject IDs in column 1.

      Genotype
        Tab or whitespace separated file where each row corresponds to single SNP, optional header
	    line can contain subject IDs, number of columns specified by --gs are copied to results file.

    OUTPUT FILE FORMAT
      Output contains the SNP information columns from the genotype file, followed by spearman correlation,
	  t statistic, p value columns. If permutations are analysed, then without the --pval, --fwer, --fdr
	  flags p values calculated on permuted datasets are reported next. Otherwise, the p value calculated
	  by permutations is printed followed by the p value adjusted for multiple testing is shown if --fwer
	  or --fdr flag is used.

###Binaries:

Binaries are compiled for 64bit linux, using dmd, gdc and ldc compilers, available here:

https://www.dropbox.com/s/vj96sci7ozvtjdm/spearman (dmd version)

https://www.dropbox.com/s/mzpcoumgdeg8mnq/spearman_gdc (gdc version)

https://www.dropbox.com/s/44j60a70kqvgj88/spearman_ldc (ldc version)

Binaries compiled with ldc or gdc are faster than those compiled with the reference implementation (dmd).
