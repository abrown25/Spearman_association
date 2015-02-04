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

import std.file : exists, File, remove;
import std.stdio : stdin, writeln;
import std.c.stdlib : exit;
import core.sys.posix.signal;

import arg_parse : Opts, giveHelp, helpString;
import calculation : rank, transform, VarianceException, covariates;
import run_analysis : noPerm, simplePerm, pvalPerm, minPerm, writeFWER, fdrCalc;
import setup_all : fileSetup, setup, F;

extern (C)
{
  void del_temp(int value);
}


version(unittest) void main() {writeln("All unit tests completed successfully.");}
 else void main(string[] args)
 {
   //set precision to either double or real (very costly in time and memory)
   alias precision = double;
   //print help string and quit if no options are given
   if (args.length == 1)
     giveHelp(helpString);

   auto opts = new Opts(cast(string[]) args);

   File[3] fileArray;

   fileSetup(fileArray, opts);

   immutable(precision[]) rankPhenotype = cast(immutable)setup!(precision)(fileArray, opts);

   if (!opts.run)
     //simple analysis with no permutations
     noPerm(fileArray, opts, rankPhenotype);
   else if (!opts.pval && !opts.min && !opts.fdr)
     //print analysis and p values for permuted datasets
     simplePerm(fileArray, opts, rankPhenotype);
   else if (!opts.min && !opts.fdr)
     //calculates permutation p values
     pvalPerm(fileArray, opts, rankPhenotype);
   else
     {
       //delete temp file when program finishes
       scope(exit)
	 {
	   if (("AndrewWantsATempFile").exists)
	     remove("AndrewWantsATempFile");
	 }
       //delete temp file if pipe closes or on ctrl c signal
       sigset(SIGPIPE, &del_temp);
       sigset(SIGINT, &del_temp);
       if (!opts.fdr)
	 {
	   //calculates family wise error rate
	   precision[] minPvalues = minPerm(fileArray, opts, rankPhenotype);
	   writeFWER(fileArray, opts, minPvalues);
	 }
       else
	 //calculates FDR
	 fdrCalc(fileArray, opts, rankPhenotype);
     }
 }
