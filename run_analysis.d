module run_analysis;

import std.array : split;
import std.conv : to, ConvException;
import std.numeric: dotProduct;
import std.range : repeat;
import std.stdio : File, stderr, stdout, writeln;
import std.string : join;
import std.array : appender;
import std.range : chunks;

version(unittest){
  import setup_all;
  import std.digest.sha;
  import std.file : remove;
  import std.range : put;
}

import calculation;

enum double EPSILON = 0.00000000001;

enum{
  phenF, genF, outF
}

class InputException : Exception {
  pure this(string s) {super(s);}
}

template readGenotype()
{
  const char[] readGenotype = "auto splitLine = split(line);

  if (skip > 0)
    fileArray[outF].write(join(splitLine[0..skip], \"\t\"), \"\t\");

  enforce(splitLine.length == nInd + skip, new InputException(\"\"));

  auto rankGenotype = to!(T[])(splitLine[skip..$]);
  transform(rank(rankGenotype));
";
}

string genErrorMsg(int x)
{
  string y = to!string(x);
  string results = "string varErr = join(\"NaN\".repeat(" ~ y ~ "), \"\t\");
  string inputErr = join(\"NA\".repeat(" ~ y ~ "), \"\t\");
  string convErr = join(\"Idiot\".repeat(" ~ y ~ "), \"\t\");
";
 return results;
}

void noPerm(T)(ref File[3] fileArray, in size_t skip, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;

  mixin(genErrorMsg(3));

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	fileArray[outF].writeln(join(to!(string[])(cor), "\t"));
      } catch(VarianceException e){
	fileArray[outF].writeln(varErr);
      } catch(InputException e){
	fileArray[outF].writeln(inputErr);
      } catch(ConvException e){
	fileArray[outF].writeln(convErr);
      }
    }
}

unittest{
  string[string] options = ["p" : "phenotype.txt", "g" : "genotype.txt",
			    "o" : "testtemp", "pid" : "T",
			    "gid" : "T", "pc" : "3", "gs" : "2"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit){
    if ("testtemp".exists)
      "testtemp".remove;
  }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  noPerm!(double)(fileArray, opts.skip, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "C7BA06FE182202627D7B882F890133171A2F0E78");
}

void simplePerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  string varErr = join("NaN".repeat(3 + nPerm), "\t");
  string inputErr = join("NA".repeat(3 + nPerm), "\t");
  string convErr = join("Idiot".repeat(3 + nPerm), "\t");

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	fileArray[outF].write(join(to!(string[])(cor), "\t"), "\t");
	foreach(ref perm; chunks(perms, nInd))
	  {
	    auto singlePerm = dotProduct(rankGenotype, perm);
	    corPvalue!(T)(singlePerm, nInd);
	    fileArray[outF].write(singlePerm, "\t");
	  }
	fileArray[outF].write("\n");
      } catch(VarianceException e){
	fileArray[outF].writeln(varErr);
      } catch(InputException e){
	fileArray[outF].writeln(inputErr);
      } catch(ConvException e){
	fileArray[outF].writeln(convErr);
      }
    }
}

unittest{
  string[string] options = ["p" : "phenotype.txt", "g" : "genotype.txt",
			    "o" : "testtemp", "pid" : "T", "perm" : "4,12",
			    "gid" : "T", "pc" : "3", "gs" : "2"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit){
    if ("testtemp".exists)
      "testtemp".remove;
  }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  simplePerm!(double)(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "9927C02CEB488C2315C099642661D99782249537");
}

void pvalPerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  mixin(genErrorMsg(4));

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	fileArray[outF].write(join(to!(string[])(cor), "\t"));
	T countBetter = 0.0;
	foreach(ref perm; chunks(perms, nInd))
	  {
	    if (fabs(dotProduct(rankGenotype, perm)) > corReal)
	      ++countBetter;
	  }
	fileArray[outF].writeln("\t", countBetter / nPerm);
      } catch(VarianceException e){
	fileArray[outF].writeln(varErr);
      } catch(InputException e){
	fileArray[outF].writeln(inputErr);
      } catch(ConvException e){
	fileArray[outF].writeln(convErr);
      }
    }
}

unittest{
  string[string] options = ["p" : "phenotype.txt", "g" : "genotype.txt",
			    "o" : "testtemp", "pid" : "T", "perm" : "1000000,12",
			    "gid" : "T", "pc" : "3", "gs" : "2", "pval" : "T"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit){
    if ("testtemp".exists)
      "testtemp".remove;
  }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  pvalPerm!(double)(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "8644F5EFDB30F9466911BF692F5E5BE9ACD38878");
}

T[] minPerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  T[] maxCor = new T[opts.number];
  mixin(genErrorMsg(4));
  maxCor[] = 0.0;
  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	fileArray[outF].writef("%a\t%g\t%g", cor[0], cor[1], cor[2]);
	T countBetter = 0.0;
	size_t i = 0;
	foreach(ref perm; chunks(perms, nInd))
	  {
	    auto singlePerm = fabs(dotProduct(rankGenotype, perm));
	    if (singlePerm > corReal)
	      ++countBetter;
	    if (singlePerm > maxCor[i])
	      maxCor[i] = singlePerm;
	    i++;
	  }
	fileArray[outF].writeln("\t", countBetter/ nPerm);
      } catch(VarianceException e){
	fileArray[outF].writeln(varErr);
      } catch(InputException e){
	fileArray[outF].writeln(inputErr);
      } catch(ConvException e){
	fileArray[outF].writeln(convErr);
      }
    }
  return maxCor;
}

void writeFWER(T)(in Opts opts, ref T[] maxCor){
  import std.algorithm : sort;
  import std.c.stdlib : exit;
  import std.range : SearchPolicy;

  File oldFile = File(opts.output ~ "temp", "r");
  File newFile;

  version(WINDOWS)
    {
      try{
	newFile = File(opts.output, "w");
      } catch(Exception e){
	stderr.writeln(e.msg);
	exit(0);
      }
    }
  else
    {
      try{
	if (opts.output != "")
	  newFile = File(opts.output, "w");
	else
	  newFile = stdout;
      } catch(Exception e){
	stderr.writeln(e.msg);
	exit(0);
      }
    }

  auto sortMax = sort!()(maxCor);
  T len = sortMax.length;

  auto headerLine = oldFile.readln();
  newFile.write(headerLine);

  auto pvalCol = split(headerLine).length - 5;
  T corStat;
  T adjusted;
  foreach(line; oldFile.byLine())
    {
      auto splitLine = split(line);
      auto tString = splitLine[pvalCol];
      if (splitLine.length > 4)
	newFile.write(join(splitLine[0..$-4], "\t"), "\t");
      if (tString == "NaN" || tString == "NA" || tString == "Idiot")
	newFile.writeln(join(tString.repeat(5), "\t"));
      else
	{
	  corStat = to!T(tString);
	  adjusted = sortMax.upperBound!(SearchPolicy.gallop)(fabs(corStat) - EPSILON).length / len;
	  newFile.writefln("%g\t%s\t%g", corStat, join(splitLine[$-3..$], "\t"), adjusted);
	}
    }
}

unittest{
  string[string] options = ["p" : "phenotype.txt", "g" : "genotype.txt",
			    "o" : "testtemp", "pid" : "T", "perm" : "100000,12",
			    "gid" : "T", "pc" : "3", "gs" : "2", "fwer" : "T"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit){
    if ("testtemp".exists)
      "testtemp".remove;
    if ("testtemptemp".exists)
      "testtemptemp".remove;
  }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  double[] minPvalues = minPerm!(double)(fileArray, opts, rankPhenotype);

  fileArray[outF].close();

  writeFWER!(double)(opts, minPvalues);

  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "B018C9BEC3EAD53106456397ED5699562490B978");
}

void fdrCalc(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  import std.algorithm : sort;
  import std.c.stdlib : exit;
  import std.math : fmin;

  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  auto permCor = appender!(T[])();
  auto realCor = appender!(T[])();
  mixin(genErrorMsg(4));

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	realCor.put(cor[0]);
	fileArray[outF].writef(join(to!(string[])(cor), "\t"));
	T countBetter = 0.0;
	foreach(ref perm; chunks(perms, nInd))
	  {
	    auto singlePerm = fabs(dotProduct(rankGenotype, perm));
	    if (singlePerm > corReal)
	      ++countBetter;
	    permCor.put(singlePerm);
	  }
	fileArray[outF].writeln("\t", countBetter/ nPerm);
      } catch(VarianceException e){
	fileArray[outF].writeln(varErr);
      } catch(InputException e){
	fileArray[outF].writeln(inputErr);
      } catch(ConvException e){
	fileArray[outF].writeln(convErr);
      }
    }

  fileArray[outF].close();

  File oldFile = File(opts.output ~ "temp", "r");
  File newFile;

  version(WINDOWS)
    {
      try{
	newFile = File(opts.output, "w");
      } catch(Exception e){
	stderr.writeln(e.msg);
	exit(0);
      }
    }
  else
    {
      try{
	if (opts.output != "")
	  newFile = File(opts.output, "w");
	else
	  newFile = stdout;
      } catch(Exception e){
	stderr.writeln(e.msg);
	exit(0);
      }
    }

  auto sortPerm = sort!()(permCor.data);
  size_t[] orderReal = new size_t[realCor.data.length];
  bestRank(orderReal, realCor.data);

  auto headerLine = oldFile.readln();
  newFile.write(headerLine);

  immutable T doubPerm = cast(immutable T) nPerm;

  size_t i = 0;
  T adjusted;
  foreach(ref line; oldFile.byLine())
    {
      auto lastString = split(line)[$-1];
      if (lastString == "NaN" || lastString == "NA" || lastString == "Idiot")
	newFile.writeln(line, "\t", lastString);
      else
	{
	  adjusted = sortPerm.upperBound!()(fabs(realCor.data[i]) - EPSILON).length / doubPerm;
	  adjusted = fmin(1, adjusted / orderReal[i]);
	  newFile.writeln(line, "\t", adjusted);
	  i++;
	}
    }
}

unittest{
  string[string] options = ["p" : "phenotype.txt", "g" : "genotype.txt",
			    "o" : "testtemp", "pid" : "T", "perm" : "100000,12",
			    "gid" : "T", "pc" : "3", "gs" : "2", "fdr" : "T"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit){
    if ("testtemp".exists)
      "testtemp".remove;
    if ("testtemptemp".exists)
      "testtemptemp".remove;
  }


  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  fdrCalc!(double)(fileArray, opts, rankPhenotype);

  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "C88E098A4CA35E036A967489C1DDA1BE5E7F51EF");
}
