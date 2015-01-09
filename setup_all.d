module setup_all;

import std.algorithm : reduce;
import std.array : split;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.exception : enforce;
import std.file : File, exists;
import std.range : iota;
import std.stdio : stdin, stderr, writeln, stdout;
import std.string : join;

import arg_parse : Opts;
import calculation : rank, transform, VarianceException, covariates;
import run_analysis : InputException;

class FileExistsException : Exception
{
  pure nothrow this(string s) {super(s);}
}
//files are stored in array, enum defines phenotype, genotype and output slots
enum F
  {
    phen, gen, out_
  }

void fileSetup(ref File[3] fileArray, Opts opts)
{
  try{
    fileArray[F.phen] = File(opts.phenotype);
  } catch(Exception e){
    stderr.writeln(e.msg);
    exit(0);
  }
  //   version(WINDOWS)
  //     {
  //       try{
  // 	fileArray[F.gen] = File(opts.genotype);
  //       } catch(Exception e){
  // 	stderr.writeln(e.msg);
  // 	exit(0);
  //       }
  //       try{
  // 	if (!opts.min || !opts.fdr)
  // 	  fileArray[F.out_] = File(opts.output, "w");
  // 	else
  // 	  {
  // 	    string outName = (opts.output ~ "temp");
  // 	    enforce(!outName.exists, new FileExistsException(("Failed to run analysis: file called " ~ outName ~ " already exists.
  // Please choose a different name for output file or delete temp file.")));
  // 	    fileArray[F.out_] = File(opts.output ~ "temp", "w");
  // 	  }
  //       } catch(Exception e){
  // 	stderr.writeln(e.msg);
  // 	exit(0);
  //       }
  //     }
  //   else
  {
    //if genotype file not specified then use stdin
    if (opts.genotype != "")
      {
	try{
	  fileArray[F.gen] = File(opts.genotype);
	} catch(Exception e){
	  stderr.writeln(e.msg);
	  exit(0);
	}
      }
    else
      fileArray[F.gen] = stdin;
    //if fwer is needed, make temp file for output, output goes to stdout if not specified
    if (opts.output == "" && !(opts.min || opts.fdr))
      fileArray[F.out_] = stdout;
    else
      {
	try{
	  if (opts.output != "" && !(opts.min || opts.fdr))
	    fileArray[F.out_] = File(opts.output, "w");
	  else if (opts.output != "")
	    {
	      string outName = (opts.output ~ "temp");
	      enforce(!outName.exists, new FileExistsException(("Failed to run analysis: file called " ~ outName ~ " already exists.
Please choose a different name for output file or delete temp file.")));
	      fileArray[F.out_] = File(opts.output ~ "temp", "w");
	    }
	  else
	    {
	      enforce(!"temp".exists, new FileExistsException("Failed to run analysis: file called temp already exists.
Please choose a different name for output file or delete temp file."));
	      fileArray[F.out_] = File("temp", "w");
	    }
	} catch(Exception e){
	  stderr.writeln(e.msg);
	  exit(0);
	}
      }
  }
}

T[] setup(T)(ref File[3] fileArray, Opts opts)
{
  import std.file : remove;
  T[] phenotype;
  string[] phenId;
  //reading in phenotype from column opts.phenC, if --pid specified, get IDs
  foreach(line; fileArray[F.phen].byLine)
    {
      auto phenLine = split(line);
      try{
	enforce(phenLine.length >= opts.phenC, new InputException(""));
	phenotype ~= to!T(phenLine[opts.phenC]);
      } catch(ConvException e){
	stderr.writeln("Failed to run analysis: Non-numeric data in phenotype");
	if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
	  (opts.output ~ "temp").remove;
	exit(0);
      } catch(InputException e){
	stderr.writeln("Failed to run analysis: column ", opts.phenC + 1, " in phenotype file doesn't exist");
	if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
	  (opts.output ~ "temp").remove;
	exit(0);
      }
      if (opts.pid)
	phenId ~= phenLine[0].idup;
    }

  string headerLine;
  string[] genId;
  string[] splitLine;
  //--gid means get the genotype IDs, we're building the header line
  if (opts.gid)
    {
      splitLine = split(fileArray[F.gen].readln);
      genId = splitLine[opts.skip..$];
      headerLine ~= join(splitLine[0..opts.skip], "\t");
      headerLine ~= "\t";
    }
  // no header line then write your own
  else if(opts.skip > 0)
    headerLine ~= "".reduce!((a, b) => a ~ "F" ~ to!string(b + 1) ~ "\t")(iota(0, opts.skip));

  headerLine ~= "Cor\tT_stat\tP";

  if (opts.run)
    {
      headerLine ~= opts.pval ? "\tPermP"
	                      : opts.min ? "\tPermP\tFWER"
	                      : opts.fdr ? "\tPermP\tFDR"
	                      : "".reduce!((a, b) => a ~ "\tP" ~ to!string(b + 1))(iota(0, opts.number));
    }
  //check IDs match
  if (!opts.nocheck && !opts.match && opts.pid && opts.gid)
    {
      if (genId != phenId)
	{
	  stderr.writeln("Failed to run analysis: Mismatched IDs");
	  if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
	    (opts.output ~ "temp").remove;
	  exit(0);
	}
      import std.algorithm : count, map, max, reduce;
      if (phenId.map!(x => phenId.count(x)).reduce!(max) > 1)
	stderr.writeln("Warning, duplicate phenotype IDs.");
      if (genId.map!(x => genId.count(x)).reduce!(max) > 1)
	stderr.writeln("Warning, duplicate genotype IDs.");
    }
  //if a covariates file is specified, regress covariates out of phenotype
  if (opts.cov != "")
    {
      try{
	static if (is(T == double))
	  covariates(opts.cov, phenotype);
	else
	  {
	    import std.algorithm : map;
	    import std.array : array;
	    double[] tempPhen = phenotype.map!(a => to!double(a)).array;
	    covariates(opts.cov, tempPhen);
	    phenotype = tempPhen.map!(a => to!T(a)).array;
	  }
      } catch(InputException e){
	stderr.writeln(e.msg);
	if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
	  (opts.output ~ "temp").remove;
	exit(0);
      } catch(ConvException){
	stderr.writeln("Failed to run analysis, non-numeric data in covariates file");
	if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
	  (opts.output ~ "temp").remove;
	exit(0);
      } catch(Exception e){
	stderr.writeln(e.msg);
	if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
	  (opts.output ~ "temp").remove;
	exit(0);
      }
    }
  if (opts.match)
    {
      //Check no individuals only in genotype file, if so, stop, otherwise rearrange phenotype so sample IDs match
      import std.range : indexed, zip;
      import std.algorithm : array, count, countUntil, filter, map, max;
      auto orderPhen = genId.map!(x => phenId.countUntil(x));
      if (orderPhen.countUntil(-1)!=-1)
	{
	  auto genIndNotPhen = zip(orderPhen, genId).filter!(a => a[0]==-1)
	                                            .map!(a => a[1])
	                                            .array;
	  if (genIndNotPhen.length==1)
	    stderr.writeln("Failed to run analysis: individual ", genIndNotPhen[0], " present in genotype file but not phenotype. Please remove this individual from genotype file.");
	  else
	    stderr.writeln("Failed to run analysis: individuals ", genIndNotPhen[0 .. ($-1)].join(", "), " and ", genIndNotPhen[$ - 1], " present in genotype file but not phenotype. Please remove these individuals from genotype file.");
	  exit(0);
	}

      enforce(phenId.indexed(orderPhen).array == genId);
      phenotype = phenotype.indexed(orderPhen).array;

      if (phenId.map!(x => phenId.count(x)).reduce!(max) > 1)
	stderr.writeln("Warning, duplicate phenotype IDs.");
      if (genId.map!(x => genId.count(x)).reduce!(max) > 1)
	stderr.writeln("Warning, duplicate genotype IDs.");
    }
  //return ranked phenotype, normalised to mean 0, sum of squares 1
  try{
    if (opts.ttest)
      transform(phenotype);
    else
      transform(rank(phenotype));
  } catch(VarianceException e){
    stderr.writeln("Failed to run analysis: Phenotype is constant");
    if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
      (opts.output ~ "temp").remove;
    exit(0);
  }

  fileArray[F.out_].writeln(headerLine);

  return phenotype;
}
