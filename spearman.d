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

import std.string, std.conv, std.stdio, std.c.stdlib;
import arg_parse, calculation, run_analysis;

void main(string[] args){

  if (args.length == 1)
    giveHelp();  
  
  string[string] options = getOpts(args[1..$]);
  Opts tempOpts = getOptions(options);
  immutable(Opts) opts = cast(immutable) tempOpts;
 
  File phenFile = File(options["p"]);

  File genFile;
  if ("g" in options)
    genFile = File(options["g"]);
  else
    genFile = stdin;

  File outFile;
  if ("o" in options)
    {
      if (opts.min)
	outFile = File(options["o"] ~ "temp", "w");
      else
	outFile = File(options["o"], "w");
    } 
  else
    {
      if (opts.min)
	outFile = File("temp","w");
      else
	outFile = stdout;
    }

  double[] phenotype;
  string[] phenId;

  foreach(line; phenFile.byLine())
    {
      auto phenLine = split(chomp(line));
      phenotype ~= to!double(phenLine[opts.phenC]);
      if ("pi" in options)
	phenId ~= phenLine[0].idup;
    }

  string headerLine;
  string[] genId;
  string[] splitLine;

  if ("gi" in options)
    {
      splitLine = split(chomp(genFile.readln()));
      genId = splitLine[opts.skip..$];
      headerLine ~= join(splitLine[0..opts.skip], "\t");
      headerLine ~= "\t";
    } 
  else if(opts.skip > 0)
    {
    for (auto j = 1; j < opts.skip + 1; j++)
      headerLine ~= ("F" ~ to!string(j) ~ "\t");
    }

  headerLine ~= "Cor\tT_stat\tP";

  if (opts.run)
    {
      if(opts.pval)
	headerLine ~= "\tPermP";
      else if(opts.min)
	headerLine ~= "\tPermP\tFWER";
      else
	{
	  for (auto j = 1; j < opts.number + 1; j++)
	    headerLine ~= ("\tP" ~ to!string(j));
	}
    }
    

  if ("pi" in options && "gi" in options && genId!=phenId)
    {
      writeln("Failed to run analysis: Mismatched IDs");
      exit(0);
    }

  double[] rankGenotype = new double[phenotype.length];;
  double[] rankTemp;
  try {
    rankTemp = transform(rank(phenotype));
  } catch(VarianceException e) {
    writeln("Failed to run analysis: Phenotype is constant");
    exit(0);
  }


  outFile.writeln(headerLine);
  

  immutable(double[]) rankPhenotype = cast(immutable)rankTemp;

  if (!opts.run)
    noPerm(phenFile, genFile, outFile, opts, rankPhenotype);
  else if (!opts.pval && !opts.min)
    simplePerm(phenFile, genFile, outFile, opts, rankPhenotype);
  else if (!opts.min)
    pvalPerm(phenFile, genFile, outFile, opts, rankPhenotype);
  else
    {
      double[] minPvalues = minPerm(phenFile, genFile, outFile, opts, rankPhenotype);
      outFile.close();
      File oldFile;
      if ("o" in options)
      	oldFile = File(options["o"] ~ "temp", "r");
      else
      	oldFile = File("temp", "r");
      writeFWER(oldFile, options, minPvalues);
    }
}
