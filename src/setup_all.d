module setup_all;

import std.algorithm : reduce;
import std.array : split;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.exception : enforce;
import std.file : exists;
import std.range : iota;
import std.stdio : File, stdin, stderr, writeln, stdout;
import std.string : join;

import arg_parse : Opts;
import calculation : rank, transform, VarianceException, covariates;
import run_analysis : InputException;

class FileExistsException : Exception
{
  pure nothrow this(string s)
  {
    super(s);
  }
}

//files are stored in array, enum defines phenotype, genotype and output slots
enum F
{
  phen,
  gen,
  out_
}

void fileSetup(ref File[3] fileArray, Opts opts)
{
  try
  {
    fileArray[F.phen] = File(opts.phenotype);
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }
  //if genotype file not specified then use stdin
  if (opts.genotype != "")
  {
    try
    {
      fileArray[F.gen] = File(opts.genotype);
    }
    catch (Exception e)
    {
      stderr.writeln(e.msg);
      exit(0);
    }
  }
  else
    fileArray[F.gen] = stdin;
  //if fwer is needed, make temp file for output, output goes to stdout if not specified
  if (opts.output == "")
    fileArray[F.out_] = stdout;
  else
    try
  {
    fileArray[F.out_] = File(opts.output, "w");
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }
}

T[] setup(T)(ref File[3] fileArray, Opts opts)
{
  import std.array : array;
  import std.algorithm : map, splitter;
  import std.range : drop, take;

  T[] phenotype;
  string[] phenId;
  //reading in phenotype from column opts.phenC, if --pid specified, get IDs
  try
  {
    if (opts.row)
    {
      int phenColumn = opts.phenC;
      if (opts.pid)
      {
        phenId = fileArray[F.phen].readln.idup.split[opts.rowskip .. $];
        phenColumn--;
      }

      phenotype = fileArray[F.phen].byLine.drop(phenColumn).front.idup.splitter.drop(opts.rowskip).map!(
        a => a.to!T).array;

      enforce(phenotype.length != 0,
        new InputException(
        "Failed to run analysis: row " ~ to!string(opts.phenC + 1) ~ " in phenotype file doesn't exist."));

      if (opts.pid)
        enforce(phenotype.length == phenId.length,
          new InputException(
          "Failed to run analysis: Phenotype IDs and data have different lengths."));
    }
    else
    {
      foreach (line; fileArray[F.phen].byLine)
      {
        auto phenLine = split(line);

        enforce(phenLine.length > opts.phenC,
          new InputException(
          "Failed to run analysis: column " ~ to!string(opts.phenC + 1) ~ " in phenotype file doesn't exist."));
        phenotype ~= to!T(phenLine[opts.phenC]);
        if (opts.pid)
          phenId ~= phenLine[0].idup;
      }
    }
  }
  catch (ConvException e)
  {
    stderr.writeln("Failed to run analysis: Non-numeric data in phenotype");
    exit(0);
  }
  catch (InputException e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }

  string headerLine;
  string[] genId;
  string[] splitLine;
  //--gid means get the genotype IDs, we're building the header line
  if (opts.gid)
  {
    splitLine = split(fileArray[F.gen].readln);
    genId = splitLine[opts.skip .. $];
    headerLine = join(splitLine[0 .. opts.skip], "\t") ~ "\t";
  }
  // no header line then write your own
  else if (opts.skip > 0)
    headerLine ~= "".reduce!((a, b) => a ~ "F" ~ to!string(b + 1) ~ "\t")(iota(0, opts.skip));

  headerLine ~= "Cor\tT_stat\tP";

  if (opts.run)
    headerLine ~= opts.pval ? "\tPermP" : opts.min ? "\tPermP\tFWER" : opts.fdr
      ? "\tPermP\tFDR" : "".reduce!((a, b) => a ~ "\tP" ~ to!string(b + 1))(iota(0,
      opts.number));

  //check IDs match
  import std.algorithm : count, countUntil, filter, map, max;
  import std.array : array;

  if (!opts.nocheck && !opts.match && opts.pid && opts.gid)
  {
    if (genId != phenId)
    {
      stderr.writeln("Failed to run analysis: Mismatched IDs");
      exit(0);
    }
    if (phenId.map!(x => phenId.count(x)).reduce!(max) > 1)
      stderr.writeln("Warning, duplicate phenotype IDs.");
    if (genId.map!(x => genId.count(x)).reduce!(max) > 1)
      stderr.writeln("Warning, duplicate genotype IDs.");
  }
  //if a covariates file is specified, regress covariates out of phenotype
  if (opts.cov != "")
  {
    try
    {
      static if (is(T == double))
        covariates(opts.cov, phenotype);
      else
      {
        double[] tempPhen = phenotype.map!(a => to!double(a)).array;
        covariates(opts.cov, tempPhen);
        phenotype = tempPhen.map!(a => to!T(a)).array;
      }
    }
    catch (InputException e)
    {
      stderr.writeln(e.msg);
      exit(0);
    }
    catch (ConvException)
    {
      stderr.writeln("Failed to run analysis, non-numeric data in covariates file");
      exit(0);
    }
    catch (Exception e)
    {
      stderr.writeln(e.msg);
      exit(0);
    }
  }
  if (opts.match)
  {
    //Check no individuals only in genotype file, if so, stop, otherwise rearrange phenotype so sample IDs match
    import std.range : indexed, zip;

    auto orderPhen = genId.map!(x => phenId.countUntil(x));
    if (orderPhen.countUntil(-1) != -1)
    {
      auto genIndNotPhen = zip(orderPhen, genId).filter!(a => a[0] == -1).map!(a => a[1]).array;
      if (genIndNotPhen.length == 1)
        stderr.writeln("Failed to run analysis: individual ",
          genIndNotPhen[0],
          " present in genotype file but not phenotype. Please remove this individual from genotype file.");
      else
        stderr.writeln("Failed to run analysis: individuals ",
          genIndNotPhen[0 .. ($ - 1)].join(", "), " and ",
          genIndNotPhen[$ - 1],
          " present in genotype file but not phenotype. Please remove these individuals from genotype file.");
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
  try
  {
    if (opts.ttest)
      transform(phenotype);
    else
    {
      auto orderBuffer = new size_t[phenotype.length];
      transform(rank(phenotype, orderBuffer));
    }
  }
  catch (VarianceException e)
  {
    stderr.writeln("Failed to run analysis: Phenotype is constant");
    exit(0);
  }

  if (!opts.noheader)
    fileArray[F.out_].writeln(headerLine);

  return phenotype;
}

version (unittest)
{
  import std.math : approxEqual;
}

unittest
{
  string[] options = [
    "dummy", "--p", "data/phenotype.txt", "--g", "data/genotype.txt", "--o",
    "/dev/null", "--ttest", "--pid", "--pc", "3"
  ];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);

  const double[] rankPhenotype = setup!(double)(fileArray, opts);

  double[] valueFromR = [
    -0.288243461282365, -0.113328088218807, 0.45012512455932,
    -0.113328088218807, 0.418848325283882, -0.250669619340339,
    -0.113872251655051, -0.483212251570077, 0.0443231235511582, 0.4493571868910
  ];
  assert(approxEqual(rankPhenotype, valueFromR));
}

unittest
{
  string[] options = [
    "dummy", "--p", "data/tphenotype.txt", "--g", "data/genotype.txt", "--o",
    "/dev/null", "--ttest", "--pid", "--pc", "3", "--row"
  ];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);

  const double[] rankPhenotype = setup!(double)(fileArray, opts);

  double[] valueFromR = [
    -0.288243461282365, -0.113328088218807, 0.45012512455932,
    -0.113328088218807, 0.418848325283882, -0.250669619340339,
    -0.113872251655051, -0.483212251570077, 0.0443231235511582, 0.4493571868910
  ];
  assert(approxEqual(rankPhenotype, valueFromR));
}

unittest
{
  string[] options = [
    "dummy", "--p", "data/tphenotype_row", "--g", "data/genotype.txt", "--o",
    "/dev/null", "--ttest", "--pid", "--pc", "3", "--row", "--row-skip", "2"
  ];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);

  const double[] rankPhenotype = setup!(double)(fileArray, opts);

  double[] valueFromR = [
    -0.288243461282365, -0.113328088218807, 0.45012512455932,
    -0.113328088218807, 0.418848325283882, -0.250669619340339,
    -0.113872251655051, -0.483212251570077, 0.0443231235511582, 0.4493571868910
  ];
  assert(approxEqual(rankPhenotype, valueFromR));
}

unittest
{
  string[] options = [
    "dummy", "--p", "data/tphenotype.txt", "--g", "data/genotype.txt", "--o",
    "/dev/null", "--ttest", "--pid", "--pc", "3", "--row", "--match", "--gid", "--gs",
    "2"
  ];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);

  const double[] rankPhenotype = setup!(double)(fileArray, opts);

  double[] valueFromR = [
    -0.288243461282365, -0.113328088218807, 0.45012512455932,
    -0.113328088218807, 0.418848325283882, -0.250669619340339,
    -0.113872251655051, -0.483212251570077, 0.4493571868910, 0.0443231235511582
  ];
  assert(approxEqual(rankPhenotype, valueFromR));
}
