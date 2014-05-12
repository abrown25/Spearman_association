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

import std.stdio; 
import calculation : rank, transform, VarianceException;
import std.algorithm : reduce;
import std.range : iota;
import std.conv : to, ConvException;

import arg_parse, run_analysis;

void main(in string[] args){

  if (args.length == 1)
    giveHelp();

  string[string] options = getOpts(args[1..$]);
  auto opts = new Opts(options);

  File phenFile;
  try{
    phenFile = File(options["p"]);
  } catch(Exception e){
    writeln(e.msg);
    exit(0);
  }

  File genFile;
  auto pGen = "g" in options;
  if (pGen)
    try{
      genFile = File(*pGen);
    } catch(Exception e){
      writeln(e.msg);
      exit(0);
    }
  else
    genFile = stdin;

  File outFile;
  auto pOut = "o" in options;
  if (!pOut && !opts.min)
      outFile = stdout;
    else 
      {
  	try{
  	  if (pOut && opts.min)
  	    outFile = File(*pOut ~ "temp", "w");
  	  else if (pOut)
  	    outFile = File(*pOut, "w");
  	  else
  	    outFile = File("temp", "w");
  	} catch(Exception e){
  	  writeln(e.msg);
  	  exit(0);
  	}
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
    headerLine ~= "".reduce!((a, b) => a ~ "F" ~ to!string(b) ~ "\t")(iota(1, opts.skip + 1));

  headerLine ~= "Cor\tT_stat\tP";

  if (opts.run)
    {
      headerLine ~= opts.pval ? "\tPermP"
	: opts.min ? "\tPermP\tFWER"
	: "".reduce!((a, b) => a ~ "\tP" ~ to!string(b))(iota(1, opts.number + 1));
    }
    
  if (opts.pid && opts.gid && !opts.nocheck && genId!=phenId)
    {
      writeln("Failed to run analysis: Mismatched IDs");
      exit(0);
    }

  try {
    transform(rank(phenotype));
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
