import std.stdio, std.file, std.string, std.c.stdlib, std.conv;

struct Opts{
  bool run = false;
  bool give_seed = false;
  bool min = false;
  bool pval = false;
  bool pid = false;
  bool gid = false;

  int number = 0;
  int seed = 0;
  int skip = 0;
  int phenC = 0;
}

immutable auto helpString = "Usage: spearman [options]:
Options:
    --help           : Display help file
    -pheno, -p       : phenotype file [default: last argument]
    -geno, -g        : genotype file [default stdin]
    -output, -o      : output file [default stdout]
    -pheno-id, -pid  : phenotype IDs are in the first column, if genotype IDs are also present then we check for mismatches
    -geno-id, -gid   : genotype IDs are in the first row, if phenotype IDs are also present then we check for mismatches
    -pheno-col, -pc  : column for phenotype values, default is 1 if phenotype IDs are not present, 2 otherwise
    -geno-skip, -gs  : column at which genotype values start, preceding columns are printed
    -perm            : calculated permuted p values, one following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed
    -pval            : report permutation p values for each test (needs perm options to be specified)
    -fwer            : calculates Family Wise Error Rate based on permutations

Input file formats:
    phenotype        : Tab or whitespace separated file with phenotype values in column specified by -pc, and optional subject IDs in column 1
    genotype         : Tab or whitespace separated file where each row corresponds to single SNP, optional header line can contain subject IDs, number of columns specified by -gs are copied to results file

Output:
    Output contains the first info columns from the genotype file, followed by spearman correlation, t statistic, p value columns. When permutations are analysed, the p value calculated by permutations is printed if the -pval flag is used. The p value calculated by permutations and then the p value adjusted for multiple testing is shown if -fwer is used. If neither flag is present, then p values for calculated on permuted datasets are reported next.
";

void giveHelp(){
  writeln(helpString);
  exit(0);
}

 
string[string] getOpts(string[] args){
  string[string] opts;
  string prefix;
  immutable string[string] optsDictParam = ["p" : "p", "phenotype" : "p", "g" : "g", "genotype" : "g", 
				  "pc" : "pc", "pheno-col" : "pc", "gs" : "gs", "geno-skip" :"gs", 
				  "perm" :"perm", "output" : "o", "o" : "o"];

  immutable string[string] optsDictFlag = ["gid" : "gid", "geno-id" : "gid", "pid" : "pid", 
				 "pheno-id" : "pid", "pval" : "pval", "fwer" : "fwer"];
 
  foreach(i, arg; args)
    {
      if (arg.startsWith("-"))
	{
	  prefix = chompPrefix(arg.idup, "-");
	  if (prefix=="-help")
	    giveHelp();
	  if (prefix in optsDictParam)
	    {
	      if (i==(args.length - 1))
		{
		  writeln("Failed to run analysis: Missing parameter for -", prefix, " option");
		  exit(0);
		}
	      opts[optsDictParam[prefix]] = args[i+1].idup;
	    }
	  else 
	    {
	      if (prefix in optsDictFlag)
		opts[optsDictFlag[prefix]] = "T";
	      else
		{
		  writeln("Failed to run analysis: Unknown command -", prefix);
		  exit(0);
		}
	    }
	}
    }
  if (!("p" in opts))
    opts["p"] = args[args.length - 1];
  if (!opts["p"].exists)
    {
      writeln("Failed to run analysis: Phenotype file missing");
      exit(0);
    }
  if ("g" in opts && !opts["g"].exists)
    {
      writeln("Failed to run analysis: Genotype file missing");
      exit(0);
    }
  return opts;
}

int getGenotypeSkip(string[string] opt){
  int skip = 0;
  if ("gs" in opt)
    try{
      skip = to!int(opt["gs"]);
    } catch (ConvException e){
      writeln("Failed to run analysis: Non-numeric argument to -geno-skip");
      exit(0);
    }
  return skip;
}

int getPhenColumn(string[string] opt){
  int phenCol;
  if ("pc" in opt)
    {
      try{
	phenCol = (to!int(opt["pc"]) - 1);
      } catch (ConvException e){
	writeln("Failed to run analysis: Non-numeric argument to -pheno-col");
	exit(0);
      }
    }
  else 
    {
      if ("pid" in opt)
	phenCol = 1;
      else 
	phenCol = 0;
    }
  return phenCol;
}

Opts getOptions(string[string] option){
  Opts opts;
  
  opts.skip = getGenotypeSkip(option);
  opts.phenC = getPhenColumn(option);

  if ("pid" in option)
    opts.pid = true;
  if ("gid" in option)
    opts.gid = true;
  
  if ("perm" in option)
    {
      opts.run = true;
      string[] value = split(option["perm"], ",");
      try{
	opts.number = to!int(value[0]);
      } catch (ConvException e){
	writeln("Failed to run analysis: Non-numeric argument to -perm");
	exit(0);
      }
      if (value.length==2)
	{
	  opts.give_seed = true;
	  try{
	    opts.seed = to!int(value[1]);
	  } catch (ConvException e){
	    writeln("Failed to run analysis: Non-numeric argument to -perm");
	    exit(0);
	  }
	}
      if ("fwer" in option)
	opts.min = true;
      if ("pval" in option)
	opts.pval = true;
    } 
  else if (("fwer" in option) || ("pval" in option))
    {
      writeln("Failed to run analysis: Permutations must be specified with the -perm flag");
      exit(0);
    }

  return opts;
}
