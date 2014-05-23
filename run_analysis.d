import std.algorithm : sort;
import std.c.stdlib : exit;
import std.file : remove;
import std.range : repeat, SearchPolicy;
import std.stdio : File, stdout, writeln;
import std.string : join;

import calculation;

class InputException : Exception {
  this(string s) {super(s);}
}

template readGenotype()
{
const char[] readGenotype = "auto splitLine = split(line);
    
  if (skip > 0)
    outFile.write(join(splitLine[0..skip], \"\t\"), \"\t\");

  if (splitLine.length != nInd + skip)
    throw new InputException(\"\");

  auto rankGenotype = to!(double[])(splitLine[skip..$]);
  transform(rank(rankGenotype));
";
}

void noPerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

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
  double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  foreach(line; genFile.byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"));
	outFile.write("\t");
	for(auto i = 0; i < nPerm; i++)
	  {
	    double singlePerm = 0;
	    for(auto j = 0; j < nInd; j++)
	      singlePerm += rankGenotype[j] * perms[i * nInd + j];
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

  double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;
  
  foreach(line; genFile.byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"));
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    double singlePerm = 0;
	    for(auto j = 0; j < nInd; j++)
	      singlePerm += rankGenotype[j] * perms[i * nInd + j];
	    corPvalue(singlePerm, nInd);
	    if (singlePerm <= cor[2])
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
  double[] minPvalues = new double[opts.number];
  double[3] cor;
  immutable size_t nInd = rankPhenotype.length;
  immutable size_t skip = opts.skip;

  double[] perms = getPerm(opts, rankPhenotype);
  immutable size_t nPerm = perms.length / nInd;

  minPvalues[] = 1.0;
  foreach(line; genFile.byLine())
    {
      try {
	mixin(readGenotype!());
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.writef("%g\t%g\t%a", cor[0], cor[1], cor[2]);
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    double singlePerm = 0;
	    for(auto j = 0; j < nInd; j++)
	      singlePerm += rankGenotype[j] * perms[i * nInd + j];
	    corPvalue(singlePerm, nInd);
	    if (singlePerm <= cor[2])
	      ++countBetter;
	    if (singlePerm < minPvalues[i])
	      minPvalues[i] = singlePerm;
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
  return minPvalues;
}

void writeFWER(in string[string] options, ref double[] minPvalues){

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

  auto sortMin = sort!()(minPvalues);
  double len = sortMin.length;

  auto headerLine = oldFile.readln();
  newFile.write(headerLine);

  auto pvalCol = split(headerLine).length - 3;
  size_t countBetter;
  double pVal;
  double adjusted;
  foreach(line; oldFile.byLine())
    {
      auto splitLine = split(line);
      auto pValString = splitLine[pvalCol];
      if (pValString=="NaN" || pValString=="NA" || pValString=="Idiot")
	newFile.writeln(line, "\t", pValString);
      else
	{
	  pVal = to!double(pValString);
	  countBetter = sortMin.upperBound!(SearchPolicy.gallopBackwards)(pVal).length;
	  adjusted = (len - countBetter) / len;
	  newFile.write(join(splitLine[0..$-2], "\t"));
	  newFile.writefln("\t%g\t%s\t%g", pVal, splitLine[$-1], adjusted);
	}
    }

  if ("o" in options)
    remove(options["o"] ~ "temp");
  else
    remove("temp");
}
