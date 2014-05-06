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
import arg_parse, calculation, noPerms;

void main(string[] args){
  string[string] options;
  Opts opts;
  double[] phenotype;
  double[] genotype;
  double singlePerm;
  string[] splitLine;
  double[] cor = new double[3];
  string[] genId;
  string[] phenId;

  auto phenFile = File();
  auto genFile = File();
  auto outFile = File();

  immutable(double[])[] perms;


  if (args.length == 1)
    giveHelp();  
  
  options = getOpts(args[1..$]);
  opts = getOptions(options);

  phenFile = File(options["p"]);

  if ("g" in options)
    genFile = File(options["g"]);
  else
    genFile = stdin;

  foreach(line; phenFile.byLine())
    {
      auto phenLine = split(chomp(line));
      phenotype ~= to!double(phenLine[opts.phenC]);
      if ("pi" in options)
	phenId ~= phenLine[0].idup;
    }

  string headerLine;
  if ("gi" in options)
    {
      splitLine = split(chomp(genFile.readln()));
      genId = splitLine[opts.skip..$];
      headerLine ~= join(splitLine[0..opts.skip], "\t") ~ "\tCor\tT_stat\tP";
      if (opts.run)
	{
	  if(opts.pval)
	    headerLine ~= "\tPermP";
	  else
	    {
	      for (auto j = 1; j < opts.number + 1; j++)
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

  double[] minPvalues = new double[opts.number];
  
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
      if (splitLine.length != phenotype.length + opts.skip)
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
