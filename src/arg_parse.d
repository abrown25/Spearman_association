module arg_parse;

import std.array : split;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.stdio : writeln, stderr;
import std.string : chompPrefix, startsWith;

class Opts{
  import std.getopt;
  //write appropriate string and quit
  bool help = false;
  bool version_ = false;
  bool ttest = false;
  //run true means run permutations
  bool run = false;
  bool give_seed = false;
  //min means calculate fwer
  bool min = false;
  bool pval = false;
  bool fdr = false;
  //phenotype and genotype ids are given
  bool pid = false;
  bool gid = false;
  //don't check ids
  bool nocheck = false;
  //rearrange to match phenotype and genotype ids;
  bool match = false;
  //number of genotype columns to skip, and phenotype column
  int skip = 0;
  int phenC = 0;
  //permutation numbers and seeds
  int number = 0;
  int seed = 0;
  //file names
  string genotype = "";
  string phenotype = "";
  string output = "";
  string cov = "";

  this(string[] args){
    try{
      getopt(
	     args,
	     "pheno|p", &phenotype,
	     "geno|g", &genotype,
	     "out|o", &output,
	     "cov|c", &cov,
	     "pheno-col|pc", &getPhenColumn,
	     "geno-skip|gs", &skip,
	     "perm", &getPerms,
	     "pheno-id|pid", &pid,
	     "geno-id|gid", &gid,
	     "fwer", &min,
	     "pval", &pval,
	     "fdr", &fdr,
	     "ttest", &ttest,
	     "nocheck", &nocheck,
	     "match", &match,
	     "help", &help,
	     "version", &version_
	     );
    } catch(Exception e){
      writeln(e.msg);
      exit(0);
    }

      checkOptions();
      if (pid && phenC == 0)
	phenC = 1;
      if (phenotype=="" && args.length > 0)
	phenotype = args[$ - 1];
 }
    private void checkOptions(){
      //write help or version string and quit
      if (help)
	giveHelp(helpString);
      if (version_)
	giveHelp(versionString);
      //check permutation strings are consistent
      if (min && pval)
	{
	  stderr.writeln("Failed to run analysis: Both -fwer and -pval flag specified");
	  exit(0);
	}
      if (min && fdr)
	{
	  stderr.writeln("Failed to run analysis: Both -fwer and -fdr flag specified");
	  exit(0);
	}
      if (fdr && pval)
	{
	  stderr.writeln("Failed to run analysis: Both -fdr and -pval flag specified");
	  exit(0);
	}
      if ((min || pval || fdr) && !run)
	{
	  stderr.writeln("Failed to run analysis: Permutations must be specified with the -perm flag");
	  exit(0);
	}
      //check ID options are consistent
      if (match && !(pid && gid))
	{
	  stderr.writeln("Failed to run analysis: --match specified but missing genotype or phenotype IDs");
	  exit(0);
	}
      if (nocheck && match)
	{
	  stderr.writeln("Failed to run analysis: Both --nocheck and --match specified");
	  exit(0);
	}
    }
    //phenotype column is column - 1
    private void getPhenColumn(string opt, string val){
      try{
	phenC = to!int(val) - 1;
	  } catch(ConvException e){
	stderr.writeln("Failed to run analysis: Non-numeric argument to -pheno-col");
	exit(0);
      }
    }
    //with perms flag, either number of perms is given (single number) or seed as well (a,b)
    private void getPerms(string opt, string val){
      run = true;
      string[] value = split(val, ",");
      try{
	number = to!int(value[0]);
      } catch(ConvException e){
	stderr.writeln("Failed to run analysis: Non-integer argument to -perm");
	exit(0);
      }
      //get seed
      if (value.length == 2)
	{
	  give_seed = true;
	  try{
	    seed = to!int(value[1]);
	  } catch(ConvException e){
	    stderr.writeln("Failed to run analysis: Non-integer argument to seed");
	    exit(0);
	  }
	}
    }
}

static immutable string helpString = "Usage: NP-GWAS [options]:
Options:
    --help             : display help information.
    --version          : display version information.
    --pheno, --p       : phenotype file [default: last argument].
    --geno, --g        : genotype file [default stdin].
    --out, --o         : output file [default stdout].
    --cov, --c         : optional covariates file, if specified analysis will be performed on the residuals, after controlling for covariates with least squares regression.
    --ttest            : runs a test of standard parametric correlation between genotype and phenotype.
    --pheno-id, --pid  : phenotype IDs are in the first column, if genotype IDs are also present then we check for mismatches.
    --geno-id, --gid   : genotype IDs are in the first row, if phenotype IDs are also present then we check for mismatches.
    --pheno-col, --pc  : column for phenotype values, default is 1 if phenotype IDs are not present, 2 otherwise.
    --geno-skip, --gs  : column at which genotype values start, preceding columns are printed.
    --perm             : calculated permuted p values, one following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed.
    --pval             : report permutation p values for each test (needs perm options to be specified).
    --fwer             : calculates the Family Wise Error Rate (FWER) based on permutations, corrected P values in the last column.
    --fdr              : calculates the False Discovery Rate (FDR) based on permutations, corrected P values in the last column.
    --nocheck          : skip check of IDs when both genotype and phenotype IDs are present.
    --match            : the program will rearrange the phenotype data so that the genotype and phenotype IDs match. If individuals are present in the genotype file only, the analysis will be halted.

Input file formats:
    phenotype          : Tab or whitespace separated file with phenotype values in column specified by --pc, and optional subject IDs in column 1.
    genotype           : Tab or whitespace separated file where each row corresponds to single SNP, optional header line can contain subject IDs, number of columns specified by --gs are copied to results file.

Output:
    Output contains the first info columns from the genotype file, followed by spearman correlation, t statistic, p value columns. When permutations are analysed, the p value calculated by permutations is printed if the --pval flag is used. The p value calculated by permutations and then the p value adjusted for multiple testing is shown if --fwer or --fdr flag is used. If none of these flags are present, then p values for calculated on permuted datasets are reported next.
";

static immutable string versionString = "NP-GWAS, version 0.9";

void giveHelp(immutable string quitString){
  writeln(quitString);
  exit(0);
}
