import std.stdio;
import arg_parse: Opts;
import calculation;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.string : join, split;
import std.range : repeat, SearchPolicy;
import std.algorithm : sort;
import std.file : remove;

class InputException : Exception {
  this(string s) {super(s);}
}

void writeError(in string error, ref File outFile, in int count){
  outFile.writeln(join(error.repeat(count), "\t"));
}

double[] readGenotype(in char[] line, ref File outFile, in int skip, in size_t indCount){
  auto splitLine = split(line);
    
  if (skip > 0)
    outFile.write(join(splitLine[0..skip], "\t"), "\t");

  if (splitLine.length != indCount + skip)
    throw new InputException("");

  auto genotype = to!(double[])(splitLine[skip..$]);
  transform(rank(genotype));

  return(genotype);
}

void noPerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  const nInd = rankPhenotype.length;
  const skip = opts.skip;

  foreach(line; genFile.byLine())
    {
      try {
	auto rankGenotype = readGenotype(line, outFile, skip, nInd);
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.writeln(join(to!(string[])(cor), "\t"));
      } catch(VarianceException e){
	writeError("NaN", outFile, 3);
      } catch(InputException e){
	writeError("NA", outFile, 3);
      } catch(ConvException e){
	writeError("Idiot", outFile, 3);
      }
    }
}


void simplePerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  const nInd = rankPhenotype.length;
  const skip = opts.skip;
  immutable(double[]) perms = cast(immutable)getPerm(opts, rankPhenotype);
  const nPerm = perms.length / nInd;
  foreach(line; genFile.byLine())
    {
      try {
	auto rankGenotype = readGenotype(line, outFile, skip, nInd);
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
	writeError("NaN", outFile, 3 + opts.number);
      } catch(InputException e){
	writeError("NA", outFile, 3 + opts.number);
      } catch(ConvException e){
	writeError("Idiot", outFile, 3 + opts.number);
      }
    }
}

void pvalPerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[3] cor;
  const nInd = rankPhenotype.length;
  const skip = opts.skip;

  immutable(double[]) perms = cast(immutable)getPerm(opts, rankPhenotype);
  const nPerm = perms.length / nInd;
  
  foreach(line; genFile.byLine())
    {
      try {
	auto rankGenotype = readGenotype(line, outFile, skip, nInd);
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"));
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    double singlePerm = 0;
	    for(auto j = 0; j < nInd; j++)
	      singlePerm += rankGenotype[j] * perms[i * nInd + j];
	    corPvalue(singlePerm, nInd);
	    if (singlePerm < cor[2])
	      ++countBetter;
	  }
	outFile.writeln("\t", countBetter / nPerm);
      } catch(VarianceException e){
	writeError("NaN", outFile, 4);
      } catch(InputException e){
	writeError("NA", outFile, 4);
      } catch(ConvException e){
	writeError("Idiot", outFile, 4);
      }
    }
}


double[] minPerm(ref File phenFile, ref File genFile, ref File outFile, in Opts opts, immutable(double[]) rankPhenotype){
  double[] minPvalues = new double[opts.number];
  double[3] cor;
  const nInd = rankPhenotype.length;
  const skip = opts.skip;

  double[] perms = getPerm(opts, rankPhenotype);
  const nPerm = perms.length / nInd;

  minPvalues[] = 1.0;
  foreach(line; genFile.byLine())
    {
      try {
	auto rankGenotype = readGenotype(line, outFile, skip, nInd);
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"));
	double countBetter = 0.0;
	for(auto i = 0; i < nPerm; i++)
	  {
	    double singlePerm = 0;
	    for(auto j = 0; j < nInd; j++)
	      singlePerm += rankGenotype[j] * perms[i * nInd + j];
	    corPvalue(singlePerm, nInd);
	    if (singlePerm < cor[2])
	      ++countBetter;
	    if (singlePerm < minPvalues[i])
	      minPvalues[i] = singlePerm;
	  }
	outFile.writeln("\t", countBetter/ nPerm);
      } catch(VarianceException e){
	writeError("NaN", outFile, 4);
      } catch(InputException e){
	writeError("NA", outFile, 4);
      } catch(ConvException e){
	writeError("Idiot", outFile, 4);
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
	  adjusted = sortMin.lowerBound!(SearchPolicy.gallopBackwards)(pVal).length / len;
	  newFile.writeln(line, "\t", adjusted);
	}
    }

  if ("o" in options)
    std.file.remove(options["o"] ~ "temp");
  else
    std.file.remove("temp");
}
