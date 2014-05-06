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
import arg_parse, calculation, run_analysis;

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
  if (!opts.run)
    noPerm(phenFile, genFile, opts, rankPhenotype);
  else
    fuck(phenFile, genFile, opts, rankPhenotype);

}
