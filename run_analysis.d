import std.string, std.conv, std.stdio, std.algorithm, std.math, std.c.stdlib, std.file, std.random, std.range;
import arg_parse, calculation;

class InputException : Exception {
  this(string s) {super(s);}
}

void writeError(string error, int count){
  for (auto j = 0; j < count; ++j)
    std.stdio.write(error,"\t");
  std.stdio.writeln(error);
}

double[] readGenotype(char[] line, Opts opts, ulong indCount){
  string[] splitLine;
  double[] genotype;
  double[] rankGenotype;

  splitLine = to!(string[])(split(line));
    
  if (opts.skip > 0)
    std.stdio.write(join(splitLine[0..opts.skip], "\t"), "\t");

  if (splitLine.length != indCount + opts.skip)
    throw new InputException("");
  else
    {
      genotype = to!(double[])(splitLine[opts.skip..$]);
      rankGenotype = transform(rank(genotype));
    }
  return(rankGenotype);
}

void noPerm(File phenFile, File genFile, Opts opts, immutable(double[]) rankPhenotype){

  string[] splitLine;
  double[] genotype;
  double[] rankGenotype;

  double[3] cor;

  
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, opts, rankPhenotype.length);
	cor = correlation(rankGenotype, rankPhenotype);
	std.stdio.writeln(join(to!(string[])(cor), "\t"));
      } catch(VarianceException e){
	writeError("NaN", 3);
      } catch(InputException e){
	writeError("NA", 3);
      }
    }
}


void simplePerm(File phenFile, File genFile, Opts opts, immutable(double[]) rankPhenotype){
  double[] rankGenotype;

  double singlePerm;
  double[3] cor;

  immutable(double[])[] perms = cast(immutable)getPerm(opts, rankPhenotype);
  
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, opts, rankPhenotype.length);
	cor = correlation(rankGenotype, rankPhenotype);
	std.stdio.write(join(to!(string[])(cor), "\t"));
	write("\t");
	foreach(i, e; perms)
	  {
	    singlePerm = corPvalue(rankGenotype, e);
	    write(singlePerm,"\t");
	  }
	write("\n");
      } catch(VarianceException e){
	writeError("NaN", 3 + opts.number);
      } catch(InputException e){
	writeError("NA", opts.number);
      }
    }
}

void pvalPerm(File phenFile, File genFile, Opts opts, immutable(double[]) rankPhenotype){
  double[] rankGenotype;

  double singlePerm;
  double[3] cor;

  immutable(double[])[] perms = cast(immutable)getPerm(opts, rankPhenotype);
  
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, opts, rankPhenotype.length);
	cor = correlation(rankGenotype, rankPhenotype);
	std.stdio.write(join(to!(string[])(cor), "\t"));
	float countBetter = 0.0;
	foreach(i, e; perms)
	  {
	    singlePerm = corPvalue(rankGenotype, e);
	    if (singlePerm < cor[2])
	      ++countBetter;
	  }
	writeln("\t", countBetter/perms.length);
      } catch(VarianceException e){
	writeError("NaN", 4);
      } catch(InputException e){
	writeError("NA", 4);
      }
    }
}


void minPerm(File phenFile, File genFile, Opts opts, immutable(double[]) rankPhenotype){
  double[] minPvalues = new double[opts.number];
  double[] rankGenotype;

  double singlePerm;
  double[] cor = new double[3];
  immutable(double[])[] perms;

  perms = cast(immutable)getPerm(opts, rankPhenotype);
  minPvalues[] = 1.0;
  
  foreach(line; genFile.byLine())
    {
      try {
	rankGenotype = readGenotype(line, opts, rankPhenotype.length);
	cor = correlation(rankGenotype, rankPhenotype);
	std.stdio.write(join(to!(string[])(cor), "\t"));
	float countBetter = 0.0;
	foreach(i, e; perms)
	  {
	    singlePerm = corPvalue(rankGenotype, e);
	    if (singlePerm < cor[2])
	      ++countBetter;
	  }
	writeln("\t", countBetter/perms.length);
      } catch(VarianceException e){
	writeError("NaN", 4);
      } catch(InputException e){
	writeError("NA", 4);
      }
    }
}
