module run_analysis;

import std.array : split;
import std.conv : to, ConvException;
import std.numeric: dotProduct;
import std.range : repeat;
import std.stdio : File, stderr, stdout, writeln;
import std.string : join;
import std.array : appender;

version(unittest){
  import setup_all;
  import std.digest.sha;
  import std.file : remove;
}

import calculation;

enum double EPSILON = 0.00000001;

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

  auto rankGenotype = to!(double[])(splitLine[skip..$]);
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

void noPerm(ref File[3] fileArray, in size_t skip, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;

  mixin(genErrorMsg(3));

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
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

  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  noPerm(fileArray, opts.skip, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "C7BA06FE182202627D7B882F890133171A2F0E78");
}

void simplePerm(ref File[3] fileArray, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  string varErr = join("NaN".repeat(3 + nPerm), "\t");
  string inputErr = join("NA".repeat(3 + nPerm), "\t");
  string convErr = join("Idiot".repeat(3 + nPerm), "\t");

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	fileArray[outF].write(join(to!(string[])(cor), "\t"), "\t");
	for(auto i = 0; i < nPerm; i++)
	  {
	    auto singlePerm = dotProduct(rankGenotype, perms[i * nInd..(i + 1) * nInd]);
	    corPvalue(singlePerm, nInd);
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

  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  simplePerm(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "9927C02CEB488C2315C099642661D99782249537");
}

void pvalPerm(ref File[3] fileArray, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  mixin(genErrorMsg(4));

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	double tReal = fabs(cor[1]) - EPSILON;
	fileArray[outF].write(join(to!(string[])(cor), "\t"));
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    auto singlePerm = dotProduct(rankGenotype, perms[i * nInd..(i + 1) * nInd]);
	    if (fabs(singlePerm * sqrt((nInd - 2) / (1 - singlePerm * singlePerm))) > tReal)
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

  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  pvalPerm(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "8644F5EFDB30F9466911BF692F5E5BE9ACD38878");
}

double[] minPerm(ref File[3] fileArray, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  double[] maxT = new double[opts.number];
  mixin(genErrorMsg(4));

  maxT[] = 0.0;
  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	double tReal = fabs(cor[1]) - EPSILON;
	fileArray[outF].writef("%g\t%a\t%g", cor[0], cor[1], cor[2]);
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    auto singlePerm = dotProduct(rankGenotype, perms[i * nInd..(i + 1) * nInd]);
	    singlePerm = fabs(singlePerm * sqrt((nInd - 2) / (1 - singlePerm * singlePerm)));
	    if (singlePerm > tReal)
	      ++countBetter;
	    if (singlePerm > maxT[i])
	      maxT[i] = singlePerm;
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
  return maxT;
}

void writeFWER(in Opts opts, ref double[] maxT){
  import std.algorithm : sort;
  import std.c.stdlib : exit;
  import std.file : remove;
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

  auto sortMax = sort!()(maxT);
  double len = sortMax.length;

  auto headerLine = oldFile.readln();
  newFile.write(headerLine);

  auto pvalCol = split(headerLine).length - 4;
  double tStat;
  double adjusted;
  foreach(line; oldFile.byLine())
    {
      auto splitLine = split(line);
      auto tString = splitLine[pvalCol];
      if (tString == "NaN" || tString == "NA" || tString == "Idiot")
	newFile.writeln(line, "\t", tString);
      else
	{
	  tStat = to!double(tString);
	  adjusted = sortMax.upperBound!(SearchPolicy.gallop)(fabs(tStat) - EPSILON).length / len;
	  newFile.write(join(splitLine[0..$-3], "\t"));
	  newFile.writefln("\t%g\t%s\t%s\t%g", tStat, splitLine[$-2], splitLine[$-1], adjusted);
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

  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  double[] minPvalues = minPerm(fileArray, opts, rankPhenotype);

  fileArray[outF].close();

  writeFWER(opts, minPvalues);

  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "B018C9BEC3EAD53106456397ED5699562490B978");
}

void fdrCalc(ref File[3] fileArray, in Opts opts, immutable(double[]) rankPhenotype){
  import std.algorithm : sort;
  import std.c.stdlib : exit;
  import std.file : remove;
  import std.math : fmin;

  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  auto permT = appender!(double[])();
  auto realT = appender!(double[])();
  mixin(genErrorMsg(4));

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	double tReal = fabs(cor[1]) - EPSILON;
	realT.put(cor[1]);
	fileArray[outF].writef(join(to!(string[])(cor), "\t"));
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    auto singlePerm = dotProduct(rankGenotype, perms[i * nInd..(i + 1) * nInd]);
	    singlePerm = fabs(singlePerm * sqrt((nInd - 2) / (1 - singlePerm * singlePerm)));
	    if (singlePerm > tReal)
	      ++countBetter;
	    permT.put(singlePerm);
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

  auto sortPerm = sort!()(permT.data);
  size_t[] orderReal = new size_t[realT.data.length];
  bestRank(orderReal, realT.data);

  auto headerLine = oldFile.readln();
  newFile.write(headerLine);

  immutable double doubPerm = cast(immutable double) nPerm;

  int i = -1;
  double adjusted;
  foreach(ref line; oldFile.byLine())
    {
      auto lastString = split(line)[$-1];
      if (lastString == "NaN" || lastString == "NA" || lastString == "Idiot")
	newFile.writeln(line, "\t", lastString);
      else
	{
	  i++;
	  adjusted = sortPerm.upperBound!()(fabs(realT.data[i]) - EPSILON).length / doubPerm;
	  adjusted = fmin(1, adjusted / orderReal[i]);
	  newFile.writeln(line, "\t", adjusted);
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


  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  fdrCalc(fileArray, opts, rankPhenotype);

  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "C88E098A4CA35E036A967489C1DDA1BE5E7F51EF");
}
