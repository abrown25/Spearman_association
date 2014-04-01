import std.stdio, std.algorithm, std.file, std.string, std.c.stdlib;

auto helpString = "Usage: spearman [options]:
Options:
-p : phenotype file [mandatory]
-g : genotype file [default stdin]
-pi : phenotype IDs are in the first column, if genotype IDs are also present then we check for mismatches
-gi : genotype IDs are in the first row, if phenotype IDs are also present then we check for mismatches
-pc : column for phenotype values, default is 1 if phenotype IDs are not present, 2 otherwise
-gs : column at which genotype values start, preceding columns are printed
-perm : calculated permuted p values, one following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed

Input file formats:
phenotype : Tab or whitespace separated file with phenotype values in column specified by -pc, and optional subject IDs in column 1
genotype : Tab or whitespace separated file where each row corresponds to single SNP, optional header line can contain subject IDs, number of columns specified by -gs are copied to results file
Output:
Output is sent to the stdout, contains the first info columns from the genotype file, followed by spearman correlation, t statistic, p value columns and then a p value column for every calculated permutation of the data
";

// to do: handle permutations, handle errors
void giveHelp(){
  writeln(helpString);
}

string[string] getOpts(string[] args){
  // options: -p phenotype, -g genotype, -pi phen ids, -gi gen ids, -perm generate permuations -pc phenotype column, -gs genotype skip
  string[string] opts;
  string prefix;

  foreach(i, arg; args){
    if (arg.startsWith("-")){
      prefix = chompPrefix(arg.idup,"-");
      switch (prefix){
      case "p":
	opts["p"] = args[i+1].idup;
	break;
      case "g":
	opts["g"] = args[i+1].idup;
	break;
      case "gi":
	opts["gi"] = "T";
	break;
      case "pi":
	opts["pi"] = "T";
	break;
      case "pc":
	opts["pc"] = args[i+1].idup;
	break;
      case "gs":
	opts["gs"] = args[i+1].idup;
	break;
      case "perm":
	opts["perm"] = args[i+1].idup;
	break;
      default:
	writeln("Unknown option: ", prefix);
	break;
      }
    }
  }
  if ("p" in opts && !opts["p"].exists){
    writeln("Phenotype file missing");
    exit(0);
  }
  if ("g" in opts && !opts["g"].exists){
    writeln("Genotype file missing");
    exit(0);
  }
  return opts;
}
