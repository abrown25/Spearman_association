module run_analysis;

import std.algorithm : sort;
import std.c.stdlib : exit;
import std.exception : enforce;
import std.file : remove;
import std.range : repeat, SearchPolicy;
import std.stdio : File, stdout, writeln;
import std.string : join;

import calculation;

const double EPSILON = 0.00000001;

class InputException : Exception {
  pure this(string s) {super(s);}
}

template readGenotype()
{
  const char[] readGenotype = "auto splitLine = split(line);
    
  if (skip > 0)
    outFile.write(join(splitLine[0..skip], \"\t\"), \"\t\");

  enforce(splitLine.length == nInd + skip, new InputException(\"\"));

  auto rankGenotype = to!(double[])(splitLine[skip..$]);
  transform(rank(rankGenotype));
";
}

void noPerm(ref File phenFile, ref File genFile, ref File outFile, in size_t skip, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;

  foreach(line; genFile.byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.writeln(join(to!(string[])(cor), "\t"));
      } catch(VarianceException e){
	outFile.writeln(join("NaN".repeat(3), "\t"));
      } catch(InputException e){
	outFile.writeln(join("NA".repeat(3), "\t"));
      } catch(ConvException e){
	outFile.writeln(join("Idiot".repeat(3), "\t"));
      }
    }
}

void simplePerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  foreach(line; genFile.byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"), "\t");
	for(auto i = 0; i < nPerm; i++)
	  {
	    auto singlePerm = dotProduct(rankGenotype, perms[i * nInd..(i + 1) * nInd]);
	    corPvalue(singlePerm, nInd);
	    outFile.write(singlePerm, "\t");
	  }
	outFile.write("\n");
      } catch(VarianceException e){
	outFile.writeln(join("NaN".repeat(3 + nPerm), "\t"));
      } catch(InputException e){
	outFile.writeln(join("NA".repeat(3 + nPerm), "\t"));
      } catch(ConvException e){
	outFile.writeln(join("Idiot".repeat(3 + nPerm), "\t"));
      }
    }
}

void pvalPerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  
  foreach(line; genFile.byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	double tReal = fabs(cor[1]) - EPSILON;
	outFile.write(join(to!(string[])(cor), "\t"));
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    auto singlePerm = dotProduct(rankGenotype, perms[i * nInd..(i + 1) * nInd]);
	    if (fabs(singlePerm * sqrt((nInd - 2) / (1 - singlePerm * singlePerm))) > tReal)
	      ++countBetter;
	  }
	outFile.writeln("\t", countBetter / nPerm);
      } catch(VarianceException e){
	outFile.writeln(join("NaN".repeat(4), "\t"));
      } catch(InputException e){
	outFile.writeln(join("NA".repeat(4), "\t"));
      } catch(ConvException e){
	outFile.writeln(join("Idiot".repeat(4), "\t"));
      }
    }
}


double[] minPerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  const double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  double[] maxT = new double[opts.number];

  maxT[] = 0.0;
  foreach(line; genFile.byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	double tReal = fabs(cor[1]) - EPSILON;
	outFile.writef("%g\t%a\t%g", cor[0], cor[1], cor[2]);
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
	outFile.writeln("\t", countBetter/ nPerm);
      } catch(VarianceException e){
	outFile.writeln(join("NaN".repeat(4), "\t"));
      } catch(InputException e){
	outFile.writeln(join("NA".repeat(4), "\t"));
      } catch(ConvException e){
	outFile.writeln(join("Idiot".repeat(4), "\t"));
      }
    }
  return maxT;
}

void writeFWER(in string[string] options, ref double[] maxT){

  File oldFile = File(options.get("o", "") ~ "temp", "r");

  File newFile;
  auto p = "o" in options;
  try{
    if (p)
      newFile = File(*p, "w");
    else
      newFile = stdout;
  } catch(Exception e){
    writeln(e.msg);
    exit(0);
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
      if (tString=="NaN" || tString=="NA" || tString=="Idiot")
	newFile.writeln(line, "\t", tString);
      else
	{
	  tStat = to!double(tString);
	  adjusted = sortMax.upperBound!(SearchPolicy.gallop)(fabs(tStat) - EPSILON).length / len;
	  newFile.write(join(splitLine[0..$-3], "\t"));
	  newFile.writefln("\t%g\t%s\t%s\t%g", tStat, splitLine[$-2], splitLine[$-1], adjusted);
	}
    }

  if ("o" in options)
    remove(options["o"] ~ "temp");
  else
    remove("temp");
}
