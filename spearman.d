/* The GPL v3 License

   Copyright (C) 2014 Genome Research Ltd.
   #
   # Author: Andrew Brown <ab25@sanger.ac.uk>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
