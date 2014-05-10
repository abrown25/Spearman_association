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

import std.string, std.conv, std.stdio, std.c.stdlib, std.file;
import arg_parse, run_analysis;
import calculation : rank, transform, VarianceException;

void main(string[] args){

  if (args.length == 1)
    giveHelp();  
  
  string[string] options = getOpts(args[1..$]);
  immutable opts = cast(immutable) getOptions(options);
 
  File phenFile;
  try{
    phenFile = File(options["p"]);
  } catch(Exception e){
    writeln(e.msg);
    exit(0);
  }

  File genFile;

  try{
    if ("g" in options)
      genFile = File(options["g"]);
    else
      genFile = stdin;
  } catch(Exception e){
    writeln(e.msg);
    exit(0);
  }

  File outFile;
  try{
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
	  outFile = File("temp", "w");
	else
	  outFile = stdout;
      }
  } catch(Exception e){
    writeln(e.msg);
    exit(0);
  }

  double[] phenotype;
  string[] phenId;

  foreach(line; phenFile.byLine())
    {
      auto phenLine = split(chomp(line));
      try{
	phenotype ~= to!double(phenLine[opts.phenC]);
      } catch(ConvException e){
	writeln("Failed to run analysis: Non-numeric data in phenotype");
	exit(0);
      }
      if (opts.pid)
	phenId ~= phenLine[0].idup;
    }

  string headerLine;
  string[] genId;
  string[] splitLine;

  if (opts.gid)
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
    

  if (opts.pid && opts.gid && genId!=phenId)
    {
      writeln("Failed to run analysis: Mismatched IDs");
      exit(0);
    }

  // double[] rankTemp;
  try {
    rank(phenotype);
    transform(phenotype);
  } catch(VarianceException e){
    writeln("Failed to run analysis: Phenotype is constant");
    exit(0);
  }
  
  immutable rankPhenotype = cast(immutable) phenotype;

  outFile.writeln(headerLine);

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
      writeFWER(options, minPvalues);
    }
}
