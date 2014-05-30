module setup_all;

import std.algorithm : reduce;
import std.array : split;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.exception : enforce;
import std.file : File, exists;
import std.range : iota;
import std.stdio : stdin, writeln, stdout;
import std.string : join;

import arg_parse : Opts;
import calculation : rank, transform, VarianceException, covariates;
import run_analysis : InputException;

class FileExistsException : Exception {
  pure nothrow this(string s) {super(s);}
}

enum{
  phenF = 0, genF = 1, outF = 2
}

void fileSetup(ref File[3] fileArray, Opts opts, string[string] options){

  try{
    fileArray[phenF] = File(options["p"]);
  } catch(Exception e){
    writeln(e.msg);
    exit(0);
  }

  auto pGen = "g" in options;
  if (pGen)
    try{
      fileArray[genF] = File(*pGen);
    } catch(Exception e){
      writeln(e.msg);
      exit(0);
    }
  else
    fileArray[genF] = stdin;

  auto pOut = "o" in options;
  if (!pOut && !opts.min)
    fileArray[outF] = stdout;
  else
    {
      try{
	if (pOut && !opts.min)
	  fileArray[outF] = File(*pOut, "w");
	else if (pOut)
	  {
	    string outName = (*pOut ~ "temp");
	    enforce(!outName.exists, new FileExistsException(("Failed to run analysis: file called " ~ outName ~ " already exists.
Please choose a different name for output file or delete temp file.")));
	    fileArray[outF] = File(*pOut ~ "temp", "w");
	  }
	else
	  {
	    enforce(!"temp".exists, new FileExistsException("Failed to run analysis: file called temp already exists.
Please choose a different name for output file or delete temp file."));
	    fileArray[outF] = File("temp", "w");
	  }
      } catch(Exception e){
	writeln(e.msg);
	exit(0);
      }
    }
}

double[] setup(ref File[3] fileArray, Opts opts, string[string] options){
  double[] phenotype;
  string[] phenId;

  foreach(line; fileArray[phenF].byLine())
    {
      auto phenLine = split(line);
      if (phenLine.length < opts.phenC)
	{
	  writeln("Failed to run analysis: column ", opts.phenC + 1, " in phenotype file doesn't exist");
	  exit(0);
	}
      try{
	phenotype ~= to!double(phenLine[opts.phenC]);
      } catch(ConvException e){
	writeln("Failed to run analysis: Non-numeric data in phenotype");
	exit(0);
      }
      if (opts.pid)
	phenId ~= phenLine[0].idup;
    }

  string headerLine;
  string[] genId;
  string[] splitLine;

  if (opts.gid)
    {
      splitLine = split(fileArray[genF].readln());
      genId = splitLine[opts.skip..$];
      headerLine ~= join(splitLine[0..opts.skip], "\t");
      headerLine ~= "\t";
    } 
  else if(opts.skip > 0)
    headerLine ~= "".reduce!((a, b) => a ~ "F" ~ to!string(b + 1) ~ "\t")(iota(0, opts.skip));

  headerLine ~= "Cor\tT_stat\tP";

  if (opts.run)
    {
      headerLine ~= opts.pval ? "\tPermP"
	: opts.min ? "\tPermP\tFWER"
	: "".reduce!((a, b) => a ~ "\tP" ~ to!string(b + 1))(iota(0, opts.number));
    }
    
  if (opts.pid && opts.gid && !opts.nocheck && genId != phenId)
    {
      writeln("Failed to run analysis: Mismatched IDs");
      exit(0);
    }

  auto p = "cov" in options;
  if (p)
    try{
      covariates(*p, phenotype);
    } catch(InputException e){
      writeln(e.msg);
      exit(0);
    } catch(ConvException){
      writeln("Failed to run analysis, non-numeric data in covariates file");
      exit(0);
    } catch(Exception e){
      writeln(e.msg);
      exit(0);
    }

  try {
    transform(rank(phenotype));
  } catch(VarianceException e){
    writeln("Failed to run analysis: Phenotype is constant");
    exit(0);
  }
  
  fileArray[outF].writeln(headerLine);

  return phenotype;
}

