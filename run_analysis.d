module run_analysis;

import std.array : split;
import std.conv : to, ConvException;
import std.numeric: dotProduct;
import std.range : repeat;
import std.stdio : File, stderr, stdout, writeln;
import std.string : join;

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

void noPerm(ref File[3] fileArray, in size_t skip, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;

  foreach(line; fileArray[genF].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	fileArray[outF].writeln(join(to!(string[])(cor), "\t"));
      } catch(VarianceException e){
	fileArray[outF].writeln(join("NaN".repeat(3), "\t"));
      } catch(InputException e){
	fileArray[outF].writeln(join("NA".repeat(3), "\t"));
      } catch(ConvException e){
	fileArray[outF].writeln(join("Idiot".repeat(3), "\t"));
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
  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  noPerm(fileArray, opts.skip, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "C7BA06FE182202627D7B882F890133171A2F0E78");
  "testtemp".remove;
}

void simplePerm(ref File[3] fileArray, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

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
	fileArray[outF].writeln(join("NaN".repeat(3 + nPerm), "\t"));
      } catch(InputException e){
	fileArray[outF].writeln(join("NA".repeat(3 + nPerm), "\t"));
      } catch(ConvException e){
	fileArray[outF].writeln(join("Idiot".repeat(3 + nPerm), "\t"));
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
  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  simplePerm(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "9927C02CEB488C2315C099642661D99782249537");
  "testtemp".remove;
}

void pvalPerm(ref File[3] fileArray, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

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
	fileArray[outF].writeln(join("NaN".repeat(4), "\t"));
      } catch(InputException e){
	fileArray[outF].writeln(join("NA".repeat(4), "\t"));
      } catch(ConvException e){
	fileArray[outF].writeln(join("Idiot".repeat(4), "\t"));
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
  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);
  pvalPerm(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  auto buffer = cast(ubyte[]) std.file.read("testtemp");
  hash.put(buffer);
  assert(toHexString(hash.finish) == "8644F5EFDB30F9466911BF692F5E5BE9ACD38878");
  "testtemp".remove;
}

double[] minPerm(ref File[3] fileArray, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  double[] maxT = new double[opts.number];

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
	fileArray[outF].writeln(join("NaN".repeat(4), "\t"));
      } catch(InputException e){
	fileArray[outF].writeln(join("NA".repeat(4), "\t"));
      } catch(ConvException e){
	fileArray[outF].writeln(join("Idiot".repeat(4), "\t"));
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

  remove(opts.output ~ "temp");
}

unittest{
  string[string] options = ["p" : "phenotype.txt", "g" : "genotype.txt",
			    "o" : "testtemp", "pid" : "T", "perm" : "100000,12",
			    "gid" : "T", "pc" : "3", "gs" : "2", "fwer" : "T"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
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
  "testtemp".remove;
}
