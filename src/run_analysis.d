module run_analysis;

import std.algorithm : count, map;
import std.array : array, split, splitter;
import std.conv : to, ConvException;
import std.file : exists, remove;
import std.numeric: dotProduct;
import std.range : repeat, retro, chunks, take, drop;
import std.stdio : File, stderr, stdout, writeln;
import std.string : join;

import calculation;
import setup_all : F;

version(unittest)
{
  import setup_all;
  import std.digest.sha;
  import std.range : put;
}


class InputException : Exception
{
  pure this(string s) {super(s);}
}

/*mixin code which reads one line of genotype file, writes out first few columns and stores genotypes
this throws errors if array too short
*/
string readGenotype(string x)
{
  const string readGenotype1 = "auto splitLine = splitter(line);

  if (skip > 0)
    ";

  const string readGenotype2 = ".write(join(splitLine.take(skip), \"\t\"), \"\t\");

  auto rankGenotype = splitLine.drop(skip).map!(to!T).array;

  enforce(rankGenotype.length == nInd, new InputException(\"\"));

  if (opts.ttest)
    transform(rankGenotype);
  else
    transform(rank(rankGenotype));
";
  return readGenotype1 ~ x ~ readGenotype2;
}

//generates error messages during compilation
string genErrorMsg(int x)
{
  string y = to!string(x);
  string results = "string varErr = join(\"NaN\".repeat(" ~ y ~ "), \"\t\");
  string inputErr = join(\"NA\".repeat(" ~ y ~ "), \"\t\");
  string convErr = join(\"NA\".repeat(" ~ y ~ "), \"\t\");
";
 return results;
}

//simple analysis, gives corr, t stat and p value
void noPerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype)
{
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  mixin(genErrorMsg(3));
  /*variance error if SNP is monomorphic (write NA), non numeric data give NaN
   */
  foreach(line; fileArray[F.gen].byLine)
    {
      try{
	mixin(readGenotype("fileArray[F.out_]"));
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

unittest
{
  string[] options = ["dummy", "--p", "data/phenotype.txt", "--g",  "data/genotype.txt",
			    "--o", "testtemp", "--pid",
			    "--gid", "--pc", "3", "--gs", "2"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit)
    {
      if ("testtemp".exists)
	"testtemp".remove;
    }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  noPerm!(double)(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start;
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "5BFF0DFB3357B92BD9A2C32FB200D5D2D0F3898F");
}

//writes out statistics as above, then p values from analysis of opts.number permuted datasets
void simplePerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype)
{
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  immutable T[] perms = cast(immutable)getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  string varErr = join("NaN".repeat(3 + nPerm), "\t");
  string inputErr = join("NA".repeat(3 + nPerm), "\t");
  string convErr = join("NA".repeat(3 + nPerm), "\t");

  foreach(line; fileArray[F.gen].byLine)
    {
      try{
	mixin(readGenotype("fileArray[F.out_]"));
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	fileArray[F.out_].write(join(to!(string[])(cor), "\t"), "\t");
	fileArray[F.out_].writeln(
				  chunks(perms, nInd).map!(a => to!string(corPvalue!(T)(dotProduct(rankGenotype, a), nInd)))
				                     .join("\t"));
      } catch(VarianceException e){
	fileArray[F.out_].writeln(varErr);
      } catch(InputException e){
	fileArray[F.out_].writeln(inputErr);
      } catch(ConvException e){
	fileArray[F.out_].writeln(convErr);
      }
    }
}

unittest
{
  string[] options = ["dummy", "--p", "data/phenotype.txt", "--g", "data/genotype.txt",
			    "-otesttemp", "--pid", "--perm", "4,12",
			    "--gid", "--pc", "3", "--gs", "2"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit)
    {
      if ("testtemp".exists)
	"testtemp".remove;
    }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  simplePerm!(double)(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start;
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "6EFF35B593B0FB4E62ECC2582E4817F601326675");
}

//calculates permutation p values
void pvalPerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype)
{
  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  immutable T[] perms = cast(immutable)getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  mixin(genErrorMsg(4));

  foreach(line; fileArray[F.gen].byLine)
    {
      try{
	mixin(readGenotype("fileArray[F.out_]"));
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	fileArray[F.out_].write(join(to!(string[])(cor), "\t"));
	
	//writes the perm p value = number of times permuted datasets produced more significant statistic / number of perms
	fileArray[F.out_].writeln("\t",
				  1.0 * chunks(perms, nInd).map!(a => fabs(dotProduct(rankGenotype, a)))
				                           .count!(a => a > corReal) / nPerm
				  );
      } catch(VarianceException e){
	fileArray[F.out_].writeln(varErr);
      } catch(InputException e){
	fileArray[F.out_].writeln(inputErr);
      } catch(ConvException e){
	fileArray[F.out_].writeln(convErr);
      }
    }
}

unittest
{
  string[] options = ["dummy", "--p", "data/phenotype.txt", "--g", "data/genotype.txt",
			    "--o", "testtemp", "--pid", "--perm", "1000000,12",
			    "--geno-id", "--pc", "3", "--gs", "2", "--pval" ];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit)
    {
      if ("testtemp".exists)
	"testtemp".remove;
    }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  pvalPerm!(double)(fileArray, opts, rankPhenotype);
  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start;
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "E7E61B089FA5E0C24700108D10DE5CEF69CE1CCC");
}

//calculates family wise error rate
void minPerm(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype)
{
  import std.algorithm : max, sort, zip;
  import std.range : iota, SearchPolicy;
  import std.c.stdlib : exit;

  auto tmpFile = File("AndrewWantsATempFile", "w");

  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  immutable T[] perms = cast(immutable)getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  T[] maxCor = new T[opts.number];
  mixin(genErrorMsg(4));
  //we need to store greatest statistic across all SNPs in maxCor
  maxCor[] = 0.0;
  foreach(line; fileArray[F.gen].byLine)
    {
      try{
	mixin(readGenotype("tmpFile"));
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	tmpFile.writef("%a\t%g\t%g", cor[0], cor[1], cor[2]);

	//	perm P value as before,and store maximum correlation for each permutation
	auto simplePerm = map!(a => to!T(fabs(dotProduct(rankGenotype, a))))(chunks(perms, nInd)).array;
	tmpFile.writeln("\t",
				  1.0 * simplePerm.count!(a => a > corReal) / nPerm);

	foreach(e; zip(iota(nPerm), simplePerm))
	  maxCor[e[0]] = max(maxCor[e[0]], e[1]);

      } catch(VarianceException e){
	tmpFile.writeln(varErr);
      } catch(InputException e){
	tmpFile.writeln(inputErr);
      } catch(ConvException e){
	tmpFile.writeln(convErr);
      }
    }

  tmpFile.close;
  auto newFile = File("AndrewWantsATempFile", "r");
  scope (exit)
    newFile.close;

  //sort stored maximum statistics
  auto sortMax = sort!()(maxCor);
  T len = sortMax.length;
  //read through old file and compare correlations to maxCor to calculate FWER
  auto corCol = opts.skip;

  T corStat;
  T adjusted;
  foreach(line; newFile.byLine)
    {
      auto splitLine = splitter(line);
      auto corString = splitLine.drop(corCol).front;
      if (corCol > 0)
	fileArray[F.out_].write(join(splitLine.take(corCol), "\t"), "\t");
      if (corString == "NaN" || corString == "NA")
	fileArray[F.out_].writeln(join(corString.repeat(5), "\t"));
      else
	{
	  corStat = to!T(corString);
	  adjusted = sortMax.upperBound!(SearchPolicy.gallop)(fabs(corStat) - EPSILON).length / len;
	  fileArray[F.out_].writefln("%g\t%s\t%g", corStat, join(splitLine.drop(corCol + 1), "\t"), adjusted);
	}
    }
}

unittest
{
  string[] options = ["dummy", "--p", "data/phenotype.txt", "--g", "data/genotype.txt",
			    "--o", "testtemp", "--pid", "--perm", "100000,12",
			    "--gid", "--pc", "3", "--gs", "2", "--fwer"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit)
    {
      if ("testtemp".exists)
	"testtemp".remove;
      if ("AndrewWantsATempFile".exists)
	"AndrewWantsATempFile".remove;
    }

  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);

  minPerm!(double)(fileArray, opts, rankPhenotype);

  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start;
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish) == "AF331E54550D37EFB955D9E8B19B2ABCA74EFB2E");
}

void fdrCalc(T)(ref File[3] fileArray, in Opts opts, immutable(T[]) rankPhenotype)
{
  import std.algorithm : makeIndex, min, reverse, sort;
  import std.c.stdlib : exit;
  import std.range : zip;

  auto tmpFile = File("AndrewWantsATempFile", "w");

  T[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;
  immutable T[] perms = cast(immutable)getPerm!(T)(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  T[] permCor;
  T[] realCor;
  mixin(genErrorMsg(4));

  foreach(line; fileArray[F.gen].byLine)
    {
      try{
	mixin(readGenotype("tmpFile"));
	cor = correlation!(T)(rankGenotype, rankPhenotype);
	T corReal = fabs(cor[0]) - EPSILON;
	realCor ~= cor[0];
	tmpFile.write(join(to!(string[])(cor), "\t"));

	auto simplePerm = map!(a => to!T(fabs(dotProduct(rankGenotype, a))))(chunks(perms, nInd)).array;
	permCor ~= simplePerm;
	tmpFile.writeln("\t",
				  1.0 * simplePerm.count!(a => a > corReal) / nPerm);
      } catch(VarianceException e){
	tmpFile.writeln(varErr);
      } catch(InputException e){
	tmpFile.writeln(inputErr);
      } catch(ConvException e){
	tmpFile.writeln(convErr);
      }
    }

  tmpFile.close;

  auto newFile = File("AndrewWantsATempFile");

  if (realCor.length == 0)
    {
      stderr.writeln("No P values to analyse.");
      exit(0);
    }

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

  reverse(orderIndex);

  if (adjusted[orderIndex[0]] > 1)
    adjusted[orderIndex[0]] = 1;

  foreach(ref e; zip(orderIndex[0 .. ($ - 1)], orderIndex[1 .. $]))
    adjusted[e[1]] = min(adjusted[e[1]], adjusted[e[0]], 1);

  size_t i = 0;
  foreach(ref line; newFile.byLine)
    {
      auto lastString = line.split.retro.front;
      if (lastString == "NaN" || lastString == "NA")
	fileArray[F.out_].writeln(line, "\t", lastString);
      else
	{
	  fileArray[F.out_].writeln(line, "\t", adjusted[i]);
	  i++;
	}
    }
}

unittest
{
  string[] options = ["dummy", "--p", "data/phenotype.txt", "--g", "data/genotype.txt",
		      "--o", "testtemp", "--pid", "--perm", "100000,12",
			    "--gid", "--pc", "5", "--gs", "2", "--fdr"];
  Opts opts = new Opts(options);
  File[3] fileArray;
  fileSetup(fileArray, opts);
  scope(exit)
    {
      if ("testtemp".exists)
	"testtemp".remove;
      if ("AndrewWantsATempFile".exists)
	"AndrewWantsATempFile".remove;
    }


  immutable(double[]) rankPhenotype = cast(immutable)setup!(double)(fileArray, opts);
  fdrCalc!(double)(fileArray, opts, rankPhenotype);

  foreach(ref e; fileArray)
    e.close;
  SHA1 hash;
  hash.start;
  put(hash, File("testtemp").byChunk(1024));
  assert(toHexString(hash.finish)=="02BE6651EEE53D7AA40A1E71B451B5E3F563D23D");
}
