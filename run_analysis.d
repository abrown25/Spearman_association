module run_analysis;

import std.array : split;
import std.conv : to, ConvException;
import std.numeric: dotProduct;
import std.range : repeat;
import std.stdio : File, stderr, stdout, writeln;
import std.string : join;
import std.array : appender;
import std.range : chunks;

import calculation;
import setup_all : F;

version(unittest){
  import setup_all;
  import std.digest.sha;
  import std.file : remove;
  import std.range : put;
}


class InputException : Exception {
  pure this(string s) {super(s);}
}

/*mixin code which reads one line of genotype file, writes out first few columns and stores genotypes
this throws errors if array too short
*/
template readGenotype()
{
  const char[] readGenotype = "auto splitLine = split(line);

  if (skip > 0)
    fileArray[F.out_].write(join(splitLine[0..skip], \"\t\"), \"\t\");

  enforce(splitLine.length == nInd + skip, new InputException(\"\"));

  auto rankGenotype = to!(T[])(splitLine[skip..$]);
  transform(rank(rankGenotype));
";
}

//generates error messages during compilation
string genErrorMsg(int x)
{
  string y = to!string(x);
  string results = "string varErr = join(\"NaN\".repeat(" ~ y ~ "), \"\t\");
  string inputErr = join(\"NA\".repeat(" ~ y ~ "), \"\t\");
  string convErr = join(\"NaN\".repeat(" ~ y ~ "), \"\t\");
";
 return results;
}

//simple analysis, gives corr, t stat and p value
void noPerm(T)(ref File[3] fileArray, in size_t skip, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;

  mixin(genErrorMsg(3));
  /*variance error if SNP is monomorphic (write NA), non numeric data give NaN
   */
  foreach(line; fileArray[F.gen].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	fileArray[F.out_].writeln(join(to!(string[])(cor), "\t"));
      } catch(VarianceException e){
	fileArray[F.out_].writeln(varErr);
      } catch(InputException e){
	fileArray[F.out_].writeln(inputErr);
      } catch(ConvException e){
	fileArray[F.out_].writeln(convErr);
      }
    }
}

unittest{
  string[] options = ["dummy", "--p", "phenotype.txt", "--g",  "genotype.txt",
			    "--o", "testtemp", "--pid",
			    "--gid", "--pc", "3", "--gs", "2"];
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
  assert(toHexString(hash.finish) == "C8DF91C2B6FBD3F5D9303AE9759E8737A5FFE248");
}

//writes out statistics as above, then p values from analysis of opts.number permuted datasets
void simplePerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  string varErr = join("NaN".repeat(3 + nPerm), "\t");
  string inputErr = join("NA".repeat(3 + nPerm), "\t");
  string convErr = join("NaN".repeat(3 + nPerm), "\t");

  foreach(line; fileArray[F.gen].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	fileArray[F.out_].write(join(to!(string[])(cor), "\t"), "\t");
	foreach(ref perm; chunks(perms, nInd))
	  {
	    auto singlePerm = dotProduct(rankGenotype, perm);
	    corPvalue!(T)(singlePerm, nInd);
	    fileArray[F.out_].write(singlePerm, "\t");
	  }
	fileArray[F.out_].write("\n");
      } catch(VarianceException e){
	fileArray[F.out_].writeln(varErr);
      } catch(InputException e){
	fileArray[F.out_].writeln(inputErr);
      } catch(ConvException e){
	fileArray[F.out_].writeln(convErr);
      }
    }
}

unittest{
  string[] options = ["dummy", "--p", "phenotype.txt", "--g", "genotype.txt",
			    "-otesttemp", "--pid", "--perm", "4,12",
			    "--gid", "--pc", "3", "--gs", "2"];
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
  assert(toHexString(hash.finish) == "DDF70C1EC8A7B9680EA039E701D7F5C07DC0EA82");
}

//calculates permutation p values
void pvalPerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  mixin(genErrorMsg(4));

  foreach(line; fileArray[F.gen].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	fileArray[F.out_].write(join(to!(string[])(cor), "\t"));
	/*countBetter stores number of times permuted datasets produced more significant statistic
	 perm p value = countBetter / number of perms*/
	T countBetter = 0.0;
	foreach(ref perm; chunks(perms, nInd))
	  {
	    if (fabs(dotProduct(rankGenotype, perm)) > corReal)
	      ++countBetter;
	  }
	fileArray[F.out_].writeln("\t", countBetter / nPerm);
      } catch(VarianceException e){
	fileArray[F.out_].writeln(varErr);
      } catch(InputException e){
	fileArray[F.out_].writeln(inputErr);
      } catch(ConvException e){
	fileArray[F.out_].writeln(convErr);
      }
    }
}

unittest{
  string[] options = ["dummy", "--p", "phenotype.txt", "--g", "genotype.txt",
			    "--o", "testtemp", "--pid", "--perm", "1000000,12",
			    "--geno-id", "--pc", "3", "--gs", "2", "--pval" ];
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
  assert(toHexString(hash.finish) == "AF535B844750CA413255FC2FDAD6C938521BEC66");
}

//calculates family wise error rate
T[] minPerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  T[] maxCor = new T[opts.number];
  mixin(genErrorMsg(4));
  //we need to store greatest statistic across all SNPs in maxCor
  maxCor[] = 0.0;
  foreach(line; fileArray[F.gen].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	fileArray[F.out_].writef("%a\t%g\t%g", cor[0], cor[1], cor[2]);
	//perm P value as before
	T countBetter = 0.0;
	size_t i = 0;
	foreach(ref perm; chunks(perms, nInd))
	  {
	    auto singlePerm = fabs(dotProduct(rankGenotype, perm));
	    if (singlePerm > corReal)
	      ++countBetter;
	    //if p value is greater than previous, store it
	    if (singlePerm > maxCor[i])
	      maxCor[i] = singlePerm;
	    i++;
	  }
	fileArray[F.out_].writeln("\t", countBetter/ nPerm);
      } catch(VarianceException e){
	fileArray[F.out_].writeln(varErr);
      } catch(InputException e){
	fileArray[F.out_].writeln(inputErr);
      } catch(ConvException e){
	fileArray[F.out_].writeln(convErr);
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

  // version(WINDOWS)
  //   {
  //     try{
  // 	newFile = File(opts.output, "w");
  //     } catch(Exception e){
  // 	stderr.writeln(e.msg);
  // 	exit(0);
  //     }
  //   }
  // else
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
    //sort stored maximum statistics
  auto sortMax = sort!()(maxCor);
  T len = sortMax.length;
  //read through old file and compare correlations to maxCor to calculate FWER
  auto headerLine = oldFile.readln();
  newFile.write(headerLine);

  auto corCol = split(headerLine).length - 5;
  T corStat;
  T adjusted;
  foreach(line; oldFile.byLine())
    {
      auto splitLine = split(line);
      auto corString = splitLine[corCol];
      if (splitLine.length > 4)
	newFile.write(join(splitLine[0..$-4], "\t"), "\t");
      if (corString == "NaN" || corString == "NA")
	newFile.writeln(join(corString.repeat(5), "\t"));
      else
	{
	  corStat = to!T(corString);
	  //this counts number of maxCor which are greater than current correlation
	  adjusted = sortMax.upperBound!(SearchPolicy.gallop)(fabs(corStat) - EPSILON).length / len;
	  newFile.writefln("%g\t%s\t%g", corStat, join(splitLine[$-3..$], "\t"), adjusted);
	}
    }
}

unittest{
  string[] options = ["dummy", "--p", "phenotype.txt", "--g", "genotype.txt",
			    "--o", "testtemp", "--pid", "--perm", "100000,12",
			    "--gid", "--pc", "3", "--gs", "2", "--fwer"];
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

  fileArray[F.out_].close();

  writeFWER!(double)(opts, minPvalues);

  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start();
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "0F4312B15A9C7903817CBB74CF8BA0BD29E07B68");
}

void fdrCalc(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype){
  import std.algorithm : makeIndex, sort;
  import std.c.stdlib : exit;

  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  const T[] perms = getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  T[] permCor;
  T[] realCor;
  mixin(genErrorMsg(4));

  foreach(line; fileArray[F.gen].byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	realCor ~= cor[0];
	fileArray[F.out_].writef(join(to!(string[])(cor), "\t"));
	T countBetter = 0.0;
	foreach(ref perm; chunks(perms, nInd))
	  {
	    auto singlePerm = fabs(dotProduct(rankGenotype, perm));
	    if (singlePerm > corReal)
	      ++countBetter;
	    permCor ~= singlePerm;
	  }
	fileArray[F.out_].writeln("\t", countBetter/ nPerm);
      } catch(VarianceException e){
	fileArray[F.out_].writeln(varErr);
      } catch(InputException e){
	fileArray[F.out_].writeln(inputErr);
      } catch(ConvException e){
	fileArray[F.out_].writeln(convErr);
      }
    }

  fileArray[F.out_].close();

  File oldFile = File(opts.output ~ "temp", "r");
  File newFile;
  // version(WINDOWS)
  //   {
  //     try{
  // 	newFile = File(opts.output, "w");
  //     } catch(Exception e){
  // 	stderr.writeln(e.msg);
  // 	exit(0);
  //     }
  //   }
  // else
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
     import std.array;
     auto sortPerm = sort!()(permCor);

     T[] adjusted = new T[realCor.length];
     auto orderIndex = new size_t[adjusted.length];
     makeIndex!("fabs(a) > fabs(b)")(realCor, orderIndex);
     adjusted[orderIndex[0]] = sortPerm.upperBound(fabs(realCor[orderIndex[0]]) - EPSILON)
                                       .length
                                       .to!T;

     foreach(e; 1..orderIndex.length)
       adjusted[orderIndex[e]] = sortPerm[0 .. (sortPerm.length - cast(size_t)adjusted[orderIndex[e - 1]])]
                                        .upperBound(fabs(realCor[orderIndex[e]]) - EPSILON)
                                        .length
                                        .to!T + adjusted[orderIndex[e - 1]];


     size_t dupcount = 0;

     foreach(i, ref e; orderIndex)
       {
	 dupcount++;
	 if (i == orderIndex.length - 1 || fabs(realCor[e]) - EPSILON > fabs(realCor[orderIndex[i + 1]]))
	   {
	     foreach(ref j; orderIndex[(i - dupcount + 1) .. (i + 1)])
	       adjusted[j] = adjusted[j] / nPerm / (i+1);
	     dupcount = 0;
	   }
       }

     adjusted[orderIndex[0]] = adjusted[orderIndex[0]] > 1 ? 1
                                                           : adjusted[orderIndex[0]];

     foreach(i, ref e; orderIndex[1 .. $])
       {
	 if (adjusted[e] < adjusted[orderIndex[i]])
	   adjusted[e] = adjusted[orderIndex[i]];
	 if (adjusted[e] > 1)
	   adjusted[e] = 1;
       }

     auto headerLine = oldFile.readln();
     newFile.write(headerLine);

     size_t i = 0;
     foreach(ref line; oldFile.byLine())
       {
	 auto lastString = split(line)[$-1];
	 if (lastString == "NaN" || lastString == "NA")
	   newFile.writeln(line, "\t", lastString);
	 else
	   {
	     newFile.writeln(line, "\t", adjusted[i]);
	     i++;
	   }
       }
}

unittest{
  string[] options = ["dummy", "--p", "phenotype.txt", "--g", "genotype.txt",
		      "--o", "testtemp", "--pid", "--perm", "100000,12",
			    "--gid", "--pc", "5", "--gs", "2", "--fdr"];
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
  assert(toHexString(hash.finish)=="E2CFD6561E8A6DCC890F017561D23956EAB5A50E");
}
