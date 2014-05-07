import std.string, std.conv, std.stdio, std.algorithm, std.range;
import arg_parse, calculation;

class InputException : Exception {
  this(string s) {super(s);}
}

void writeError(string error, File outFile, int count){
  for (auto j = 0; j < count - 1; j++)
    outFile.write(error, "\t");
  outFile.writeln(error);
}

double[] readGenotype(char[] line, File outFile, Opts opts, ulong indCount){
  string[] splitLine;
  double[] genotype;
  double[] rankGenotype;

  splitLine = to!(string[])(split(line));
    
  if (opts.skip > 0)
    outFile.write(join(splitLine[0..opts.skip], "\t"), "\t");

  if (splitLine.length != indCount + opts.skip)
    throw new InputException("");
  else
    {
      genotype = to!(double[])(splitLine[opts.skip..$]);
      rankGenotype = transform(rank(genotype));
    }
  return(rankGenotype);
}

void noPerm(File phenFile, File genFile, File outFile, Opts opts, immutable(double[]) rankPhenotype){
  double[] rankGenotype;
  double[3] cor;

  
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, outFile, opts, rankPhenotype.length);
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


void simplePerm(File phenFile, File genFile, File outFile, Opts opts, immutable(double[]) rankPhenotype){
  double[] rankGenotype;
  double singlePerm;
  double[3] cor;

  immutable(double[])[] perms = cast(immutable)getPerm(opts, rankPhenotype);
  
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, outFile, opts, rankPhenotype.length);
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"));
	outFile.write("\t");
	foreach(i, e; perms)
	  {
	    singlePerm = corPvalue(rankGenotype, e);
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

void pvalPerm(File phenFile, File genFile, File outFile, Opts opts, immutable(double[]) rankPhenotype){
  double[] rankGenotype;
  double singlePerm;
  double[3] cor;

  immutable(double[])[] perms = cast(immutable)getPerm(opts, rankPhenotype);
  
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, outFile, opts, rankPhenotype.length);
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"));
	double countBetter = 1.0;
	foreach(i, e; perms)
	  {
	    singlePerm = corPvalue(rankGenotype, e);
	    if (singlePerm < cor[2])
	      ++countBetter;
	  }
	outFile.writeln("\t", countBetter/(perms.length + 1));
      } catch(VarianceException e){
	writeError("NaN", outFile, 4);
      } catch(InputException e){
	writeError("NA", outFile, 4);
      } catch(ConvException e){
	writeError("Idiot", outFile, 4);
      }
    }
}


double[] minPerm(File phenFile, File genFile, File outFile, Opts opts, immutable(double[]) rankPhenotype){
  double[] minPvalues = new double[opts.number];
  double[] rankGenotype;
  double singlePerm;
  double[3] cor;

  immutable(double[])[] perms = cast(immutable)getPerm(opts, rankPhenotype);

  minPvalues[] = 1.0;
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, outFile, opts, rankPhenotype.length);
	cor = correlation(rankGenotype, rankPhenotype);
	outFile.write(join(to!(string[])(cor), "\t"));
	double countBetter = 1.0;
	foreach(i, e; perms)
	  {
	    singlePerm = corPvalue(rankGenotype, e);
	    if (singlePerm < cor[2])
	      ++countBetter;
	    if (singlePerm < minPvalues[i])
	      minPvalues[i] = singlePerm;
	  }
	outFile.writeln("\t", countBetter/(perms.length + 1));
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

void writeFWER(string[string] options, double[] minPvalues){

  File oldFile;
  if ("o" in options)
    oldFile = File(options["o"] ~ "temp", "r");
  else
    oldFile = File("temp", "r");

  File newFile;
  if ("o" in options) 
    newFile = File(options["o"], "w");
  else
    newFile = stdout;

  auto sortMin = sort!()(minPvalues);
  double len = sortMin.length + 1.0;

  newFile.write(oldFile.readln());

  double pVal;
  double adjusted;
  foreach(line; oldFile.byLine()){
    auto splitLine = split(chomp(line));
    auto lastEntry = splitLine[(splitLine.length - 2)];
    if (lastEntry=="NaN" || lastEntry=="NA" || lastEntry=="Idiot")
      newFile.writeln(line, "\t", lastEntry);
    else
      {
    	pVal = to!double(lastEntry);
    	adjusted = (sortMin.lowerBound!(SearchPolicy.gallopBackwards)(pVal).length + 1) / len;
    	newFile.writeln(line, "\t", adjusted);
      }
  }
  if ("o" in options)
    std.file.remove(options["o"] ~ "temp");
  else
    std.file.remove("temp");
}
