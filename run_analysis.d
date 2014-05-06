import std.string, std.conv, std.stdio, std.algorithm, std.math, std.c.stdlib, std.file, std.random, std.range;
import arg_parse, calculation;

class InputException : Exception {
  this(string s) {super(s);}
}

// double[] readGenotype(string line){
//   string[] splitLine;
//   double[] genotype;
//   double[] rankGenotype;


void noPerm(File phenFile, File genFile, Opts opts, immutable(double[]) rankPhenotype){

  string[] splitLine;
  double[] genotype;
  double[] rankGenotype;

  double[3] cor;

  
  foreach(line; genFile.byLine())
    {
      splitLine = to!(string[])(split(line));
    
      if (opts.skip > 0)
    	std.stdio.write(join(splitLine[0..opts.skip], "\t"), "\t");

      if (splitLine.length != rankPhenotype.length + opts.skip)
	writeln("NA\tNA\tNA");
      else
    	{
    	  genotype = to!(double[])(splitLine[opts.skip..$]);
    	  try {
    	    rankGenotype = transform(rank(genotype));
    	    cor = correlation(rankGenotype, rankPhenotype);
    	    std.stdio.writeln(join(to!(string[])(cor), "\t"));
    	  } 
	  catch(VarianceException e){
	    writeln("NaN\tNaN\tNaN");
	  }
	}
    }
}


// void pvalPerm(File phenFile, File genFile, Opts opts, immutable(double[]) rankPhenotype){
//   string[] splitLine;
//   double[] genotype;
//   double[] rankGenotype;

//   double singlePerm;
//   double[3] cor;

//   immutable(double[])[] perms = cast(immutable)getPerm(opts, rankPhenotype);
  
//   foreach(line; genFile.byLine())
//     {
//       splitLine = to!(string[])(split(line));
    
//       if (opts.skip > 0)
//     	std.stdio.write(join(splitLine[0..opts.skip], "\t"), "\t");
//       if (splitLine.length != rankPhenotype.length + opts.skip)
//     	{
//     	  for (auto j=0; j < opts.number + 2; j++)
//     	    std.stdio.write("NA\t");
//     	  writeln("NA");
//     	}
//     	  genotype = to!(double[])(splitLine[opts.skip..$]);
//     	  try {
//     	    rankGenotype = transform(rank(genotype));
//     	    cor = correlation(rankGenotype, rankPhenotype);
//     	    std.stdio.write(join(to!(string[])(cor), "\t"));
//     	    	float countBetter = 0.0;
//     	    	foreach(i, e; perms)
//     	    	  {
//     	    	    singlePerm = corPvalue(rankGenotype, e);
//     	    	    if (singlePerm < cor[2])
//     	    	      ++countBetter;
//     	    	  }
//     	    	writeln("\t", countBetter/perms.length);
//     		foreach(i, e; perms)
//     		  {
//     		    singlePerm = corPvalue(rankGenotype, e);
//     		    if (singlePerm < minPvalues[i])
//     		      minPvalues[i] = singlePerm;
//     		    std.stdio.write("\t", singlePerm);
//     		  }
//     		std.stdio.write("\n");
//     	      }
//     	  } catch(VarianceException e){
//     	    if (opts.pval)
//     	      writeln("NaN\tNaN\tNaN\tNaN");
//     	    else
//     	      {
//     		for (auto j=0; j < opts.number + 2; j++)
//     		  std.stdio.write("NaN\t");
//     		writeln("NaN");
//     	      }
// 	  }
	 
//     	}
//     }

// }
void fuck(File phenFile, File genFile, Opts opts, immutable(double[]) rankPhenotype){
  string[] splitLine;
  double[] minPvalues = new double[opts.number];
  double[] genotype;
  double[] rankGenotype;

  double singlePerm;
  double[] cor = new double[3];
  immutable(double[])[] perms;

  if (opts.run)
    {
      perms = cast(immutable)getPerm(opts, rankPhenotype);
      if (opts.min)
  	{
  	  minPvalues[] = 1.0;
  	}
    }
  
  foreach(line; genFile.byLine())
    {
      splitLine = to!(string[])(split(line));
    
      if (opts.skip > 0)
    	std.stdio.write(join(splitLine[0..opts.skip], "\t"), "\t");
      if (splitLine.length != rankPhenotype.length + opts.skip)
    	{
    	  for (auto j=0; j < opts.number + 2; j++)
    	    std.stdio.write("NA\t");
    	  writeln("NA");
    	}
      else
    	{
    	  genotype = to!(double[])(splitLine[opts.skip..$]);
    	  try {
    	    rankGenotype = transform(rank(genotype));
    	    cor = correlation(rankGenotype, rankPhenotype);
    	    std.stdio.write(join(to!(string[])(cor), "\t"));
    	    if (opts.pval)
    	      {
    	    	float countBetter = 0.0;
    	    	foreach(i, e; perms)
    	    	  {
    	    	    singlePerm = corPvalue(rankGenotype, e);
    	    	    if (singlePerm < minPvalues[i])
    	    	      minPvalues[i] = singlePerm;
    	    	    if (singlePerm < cor[2])
    	    	      ++countBetter;
    	    	  }
    	    	writeln("\t", countBetter/perms.length);
    	      }
    	    else	  
    	      {
    		foreach(i, e; perms)
    		  {
    		    singlePerm = corPvalue(rankGenotype, e);
    		    if (singlePerm < minPvalues[i])
    		      minPvalues[i] = singlePerm;
    		    std.stdio.write("\t", singlePerm);
    		  }
    		std.stdio.write("\n");
    	      }
    	  } catch(VarianceException e){
    	    if (opts.pval)
    	      writeln("NaN\tNaN\tNaN\tNaN");
    	    else
    	      {
    		for (auto j=0; j < opts.number + 2; j++)
    		  std.stdio.write("NaN\t");
    		writeln("NaN");
    	      }
	  }
	 
    	}
    }

}
