/* The MIT License

   Copyright (C) 2014 Genome Research Ltd.
   #
   # Author: Andrew Brown <ab25@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

import std.string, std.conv, std.stdio, std.algorithm, std.math, std.c.stdlib, std.file, std.random, std.range;
import arg_parse, calculation;

void main(string[] args){
  string[string] options;
  double[] phenotype;
  double[] genotype;
  double singlePerm;
  string[] splitLine;
  double[] cor = new double[3];
  int skip = 0;
  int phenColumn = 0;
  bool pvalCalc = false;
  bool permRun = false;
  string[] genId;
  string[] phenId;
  PermOpts permOptions;
  immutable(double[])[] perms;

  auto phenFile = File();
  auto genFile = File();

  if (args.length == 1)
    giveHelp();  
  
  options = getOpts(args[1..$]);

  phenFile = File(options["p"]);
  skip = getGenotypeSkip(options);
  phenColumn = getPhenColumn(options);
  permOptions = getPermOptions(options);
  if (permOptions.run)
    permRun = true;

  if ("pval" in options)
    pvalCalc = true;

  if ("g" in options)
    genFile = File(options["g"]);
  else
    genFile = stdin;

  foreach(line; phenFile.byLine())
    {
      auto phenLine = split(chomp(line));
      phenotype ~= to!double(phenLine[phenColumn]);
      if ("pi" in options)
	phenId ~= phenLine[0].idup;
    }

  string headerLine;
  if ("gi" in options)
    {
      splitLine = split(chomp(genFile.readln()));
      genId = splitLine[skip..$];
      headerLine ~= join(splitLine[0..skip], "\t") ~ "\tCor\tT_stat\tP";
      if (permRun)
	{
	  if(pvalCalc)
	    headerLine ~= "PermP";
	  else
	    {
	      for (auto j = 1; j < permOptions.number + 1; j++)
		headerLine ~= "\tP" ~ to!string(j);
	    }
	}
    }

  if ("pi" in options && "gi" in options && genId!=phenId)
    {
      writeln("Mismatched ID");
      exit(0);
    }

  double[] rankGenotype = new double[phenotype.length];;
  double[] rankTemp;
  try {
    rankTemp = transform(rank(phenotype));
  } catch(VarianceException e) {
    writeln("Phenotype is constant");
    exit(0);
  }

  writeln(headerLine);
  
  immutable(double[]) rankPhenotype = cast(immutable)rankTemp;

  double[] minPvalues = new double[permOptions.number];
  
  if (permOptions.run)
    {
      perms = cast(immutable)getPerm(permOptions, rankPhenotype);
      if (permOptions.min)
	{
	  minPvalues[] = 1.0;
	}
    }
  
  foreach(line; genFile.byLine())
    {
      splitLine = to!(string[])(split(line));
      if (skip > 0)
	std.stdio.write(join(splitLine[0..skip], "\t"), "\t");
      if (splitLine.length != phenotype.length + skip)
	{
	  for (auto j=0; j < permOptions.number + 2; j++)
	    std.stdio.write("NA\t");
	  writeln("NA");
	}
      else
	{
	  genotype = to!(double[])(splitLine[skip..$]);
	  try {
	    rankGenotype = transform(rank(genotype));
	    cor = correlation(rankGenotype, rankPhenotype);
	    std.stdio.write(join(to!(string[])(cor), "\t"));
	    if (pvalCalc)
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
	    if (pvalCalc)
	      writeln("NaN\tNaN\tNaN\tNaN");
	    else
	      {
		for (auto j=0; j < permOptions.number + 2; j++)
		  std.stdio.write("NaN\t");
		writeln("NaN");
	      }
	  }
	 
	}
    }
}
