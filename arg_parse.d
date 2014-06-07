module arg_parse;

import std.array : split;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.stdio : writeln, stderr;
import std.string : chompPrefix, startsWith;

class Opts{
  bool run = false;
  bool give_seed = false;
  bool min = false;
  bool pval = false;
  bool pid = false;
  bool gid = false;
  bool nocheck = false;
  bool fdr = false;
  int number = 0;
  int seed = 0;
  int skip = 0;
  int phenC = 0;
  string genotype = "";
  string phenotype = "";
  string output = "";
  string cov = "";

  this(string[string] option){
    getGenotypeSkip(option);
    getPhenColumn(option);
    getFlags(option);
    getPerms(option);
    getFiles(option);
  }

  private void getGenotypeSkip(string[string] option){
    auto p = "gs" in option;
    try{
      skip = p ? to!int(*p) : 0;
    } catch (ConvException e){
      stderr.writeln("Failed to run analysis: Non-numeric argument to -geno-skip");
      exit(0);
    }
  }

  private void getPhenColumn(string[string] option){
    auto p = "pc" in option;
  try{
    phenC = p ? to!int(*p) - 1
      : "pid" in option ? 1 : 0;
  } catch (ConvException e){
    stderr.writeln("Failed to run analysis: Non-numeric argument to -pheno-col");
    exit(0);
  }
}

  private void getFlags(string[string] option){
    if ("nocheck" in option)
      nocheck = true;
    if ("pid" in option)
      pid = true;
    if ("gid" in option)
      gid = true;
    if ("fwer" in option)
      min = true;
    if ("pval" in option)
      pval = true;
    if ("fdr" in option)
      fdr = true;
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
    if ((min || pval || fdr) && !("perm" in option))
      {
	stderr.writeln("Failed to run analysis: Permutations must be specified with the -perm flag");
	exit(0);
      }
  }

  private void getPerms(string[string] option){
    if (auto p = "perm" in option)
      {
	run = true;
	string[] value = split(*p, ",");
	try{
	  number = to!int(value[0]);
	} catch (ConvException e){
	  stderr.writeln("Failed to run analysis: Non-integer argument to -perm");
	  exit(0);
	}
	if (value.length == 2)
	  {
	    give_seed = true;
	    try{
	      seed = to!int(value[1]);
	    } catch (ConvException e){
	      stderr.writeln("Failed to run analysis: Non-integer argument to seed");
	      exit(0);
	    }
	  }
      }
  }

  private void getFiles(string[string] option){
    phenotype = option["p"];
    if (auto p1 = "g" in option)
      genotype = *p1;
    if (auto p2 = "o" in option)
      output = *p2;
    if (auto p3 = "cov" in option)
      cov = *p3;
  }
}

static immutable string helpString = "Usage: spearman [options]:
Options:
    --help           : display help file
    -pheno, -p       : phenotype file [default: last argument]
    -geno, -g        : genotype file [default stdin]
    -output, -o      : output file [default stdout]
    -cov             : optional covariates file, if specified analysis will be performed on the residuals, after controlling for covariates with least squares regression
    -pheno-id, -pid  : phenotype IDs are in the first column, if genotype IDs are also present then we check for mismatches
    -geno-id, -gid   : genotype IDs are in the first row, if phenotype IDs are also present then we check for mismatches
    -pheno-col, -pc  : column for phenotype values, default is 1 if phenotype IDs are not present, 2 otherwise
    -geno-skip, -gs  : column at which genotype values start, preceding columns are printed
    -perm            : calculated permuted p values, one following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed
    -pval            : report permutation p values for each test (needs perm options to be specified)
    -fwer            : calculates Family Wise Error Rate (FWER) based on permutations
    -fdr             : calculates False Discovery Rate (FDR) based on permutations
    -nocheck         : skip check of IDs when both genotype and phenotype IDs are present

Input file formats:
    phenotype        : Tab or whitespace separated file with phenotype values in column specified by -pc, and optional subject IDs in column 1
    genotype         : Tab or whitespace separated file where each row corresponds to single SNP, optional header line can contain subject IDs, number of columns specified by -gs are copied to results file

Output:
    Output contains the first info columns from the genotype file, followed by spearman correlation, t statistic, p value columns. When permutations are analysed, the p value calculated by permutations is printed if the -pval flag is used. The p value calculated by permutations and then the p value adjusted for multiple testing is shown if -fwer is used. If neither flag is present, then p values for calculated on permuted datasets are reported next.
";

static immutable string versionString = "NP-GWAS, version 0.1.";

void giveHelp(immutable string quitString){
  writeln(quitString);
  exit(0);
}

string[string] getOpts(in string[] args){
  string[string] opts;
  string prefix;
  immutable string[string] optsDictParam = ["p" : "p", "phenotype" : "p", "g" : "g", "genotype" : "g",
					    "cov" : "cov", "pc" : "pc", "pheno-col" : "pc", "gs" : "gs",
					    "geno-skip" :"gs", "perm" :"perm", "output" : "o", "o" : "o"];

  immutable string[string] optsDictFlag = ["gid" : "gid", "geno-id" : "gid", "pid" : "pid", "nocheck" : "nocheck",
					   "pheno-id" : "pid", "pval" : "pval", "fwer" : "fwer", "fdr" : "fdr"];

  foreach(i, arg; args)
    {
      if (arg.startsWith("-"))
	{
	  prefix = chompPrefix(arg.idup, "-");
	  if (prefix == "-help")
	    giveHelp(helpString);
	  if (prefix == "-version")
	    giveHelp(versionString);
	  if (auto pParam = prefix in optsDictParam)
	    {
	      if (i == (args.length - 1))
		{
		  stderr.writeln("Failed to run analysis: Missing parameter for -", prefix, " option");
		  exit(0);
		}
	      opts[*pParam] = args[i+1].idup;
	    }
	  else
	    {
	      if (auto pFlag = prefix in optsDictFlag)
		opts[*pFlag] = "T";
	      else
		{
		  stderr.writeln("Failed to run analysis: Unknown command -", prefix);
		  exit(0);
		}
	    }
	}
    }
  if (!("p" in opts))
    opts["p"] = args[args.length - 1];
  return opts;
}
