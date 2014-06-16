##NP-GWAS : Fast, non parametric tests for genetic associations with permutations

###Introduction:

The purpose of this program is to allow genome wide analyses of non normally distributed data to be performed quickly and easily. Fast routines to calculate permutations are included, with options to calculate p values adjusted for multiple testing, either the family wise error rate (FWER) or the typically less conservative false discovery rate (FDR). The program is designed to be flexible with the data formats it will accept: genotype data should come with individuals are columns and SNPs in rows, with an optional header line and number of preceding columns giving SNP information. The program has been designed to work with linux pipes, allowing it to be seamlessly chained with programs to manipulate SNP data such as tabix and vcf-tools. Phenotype information should be given in a specified column of a text file: if there are known confounders, these can be passed as a matrix and controlled out using linear regression.

###Usage:

    Usage: spearman [options]:

    Options:
        --help             : display help file
        --pheno, --p       : phenotype file [default: last argument]
        --geno, --g        : genotype file [default stdin]
        --out, --o         : output file [default stdout]
        --cov, --c         : optional covariates file, if specified analysis will be performed on the residuals, after controlling for covariates with least squares regression
        --pheno-id, --pid  : phenotype IDs are in the first column, if genotype IDs are also present then we check for mismatches
        --geno-id, --gid   : genotype IDs are in the first row, if phenotype IDs are also present then we check for mismatches
        --pheno-col, --pc  : column for phenotype values, default is 1 if phenotype IDs are not present, 2 otherwise
        --geno-skip, --gs  : column at which genotype values start, preceding columns are printed
        --perm             : calculated permuted p values, one following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed
        --pval             : report permutation p values for each test (needs perm options to be specified)
        --fwer             : calculates the Family Wise Error Rate (FWER) based on permutations, corrected P values in the last column
        --fdr              : calculates the False Discovery Rate (FDR) based on permutations, corrected P values in the last column
        --nocheck          : skip check of IDs when both genotype and phenotype IDs are present

    Input file formats:
        phenotype          : Tab or whitespace separated file with phenotype values in column specified by --pc, and optional subject IDs in column 1
        genotype           : Tab or whitespace separated file where each row corresponds to single SNP, optional header line can contain subject IDs, number of columns specified by --gs are copied to results file
    Output:
        Output contains the first info columns from the genotype file, followed by spearman correlation, t statistic, p value columns. When permutations are analysed, the p value calculated by permutations is printed if the --pval flag is used. The p value calculated by permutations and then the p value adjusted for multiple testing is shown if --fwer or --fdr flag is used. If none of these flags are present, then p values for calculated on permuted datasets are reported next.

###Binaries:

Binaries compiled with ldc or gdc are faster than those copmiled with the reference implementation, but may have slightly biased permutations.

###Issues:

The D language changed its random number algorithms to ensure better coverage of the chosen region. The gdc and ldc compilers have yet to catch up to this change, though code compiled with these is significantly faster. These issues mean that the code will fail unittests when using these compilers.