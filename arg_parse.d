import std.stdio, std.algorithm, std.file, std.string, std.c.stdlib, std.conv;

struct PermOpts{
  bool run = false;
  bool give_seed = false;
  int number = 0;
  int seed = 0;
}

auto helpString = "Usage: spearman [options]:
Options:
--help          : Display help file
-pheno, -p      : phenotype file [mandatory]
-geno, -g       : genotype file [default stdin]
-pheno-id, -pi  : phenotype IDs are in the first column, if genotype IDs are also present then we check for mismatches
-geno-id, -gi   : genotype IDs are in the first row, if phenotype IDs are also present then we check for mismatches
-pheno-col, -pc : column for phenotype values, default is 1 if phenotype IDs are not present, 2 otherwise
-geno-skip, -gs : column at which genotype values start, preceding columns are printed
-perm           : calculated permuted p values, one following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed

Input file formats:
phenotype       : Tab or whitespace separated file with phenotype values in column specified by -pc, and optional subject IDs in column 1
genotype        : Tab or whitespace separated file where each row corresponds to single SNP, optional header line can contain subject IDs, number of columns specified by -gs are copied to results file
Output:
Output is sent to the stdout, contains the first info columns from the genotype file, followed by spearman correlation, t statistic, p value columns and then a p value column for every calculated permutation of the data
";

void giveHelp(){
  writeln(helpString);
  exit(0);
}

 
string[string] getOpts(string[] args){
  string[string] opts;
  string prefix;
  string[string] optsDictParam = ["p" : "p", "phenotype" : "p", "g" : "g", "genotype" : "g", 
		      "pc" : "pc", "pheno-col" : "pc", "gs" : "gs", "geno-skip" :"gs", 
		      "perm" :"perm"];

  string[string] optsDictFlag = ["gi" : "gi", "geno-id" : "gi", "pi" : "pi", "pheno-id" : "pi"];
 
  foreach(i, arg; args)
    {
      if (arg.startsWith("-"))
	{
	  prefix = chompPrefix(arg.idup,"-");
	  if (prefix=="-help")
	    giveHelp();
	  if (prefix in optsDictParam)
	    {
	      if (i==(args.length - 1))
		{
		  writeln("Missing parameter for -", prefix, " option");
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
		  writeln("Unknown parameter");
		  exit(0);
		}
	    }
	}
    }
  if (!("p" in opts) || !opts["p"].exists)
    {
      writeln("Phenotype file missing");
      exit(0);
    }
  if ("g" in opts && !opts["g"].exists)
    {
      writeln("Genotype file missing");
      exit(0);
    }
  return opts;
}

int getGenotypeSkip(string[string] opt){
  int skip = 0;
  if ("gs" in opt)
    skip = to!int(opt["gs"]);
  return skip;
}

int getPhenColumn(string[string] opt){
  int phenCol;
  if ("pc" in opt)
    phenCol = (to!int(opt["pc"]) - 1);
  else 
    {
      if ("pi" in opt)
	phenCol = 1;
      else 
	phenCol = 0;
    }
  return phenCol;
}

PermOpts getPermOptions(string[string] option){
  PermOpts opts;
  string[] value;
  if ("perm" in option)
    {
      opts.run = true;
      value = split(option["perm"], ",");
      opts.number = to!int(value[0]);
      if (value.length==2)
	{
	  opts.give_seed = true;
	  opts.seed = to!int(value[1]);
	}
    }
  return opts;
}
