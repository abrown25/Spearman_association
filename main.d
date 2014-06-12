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

import arg_parse : Opts, giveHelp, getOpts, helpString;
import calculation : rank, transform, VarianceException, covariates;
import run_analysis : noPerm, simplePerm, pvalPerm, minPerm, writeFWER, fdrCalc;
import setup_all : fileSetup, setup, F;

extern (C) {
  void del_temp(int value);
}


version(unittest) void main() {writeln("All unit tests completed successfully.");}
 else void main(in string[] args)
 {
   alias precision = double;

   if (args.length == 1)
     giveHelp(helpString);

   string[string] options = getOpts(args[1..$]);
   auto opts = new Opts(options);

   File[3] fileArray;

   scope(exit){
     if ((opts.min || opts.fdr) && (opts.output ~ "temp").exists)
       remove((opts.output ~ "temp"));
   }

   fileSetup(fileArray, opts);

   if ((opts.min || opts.fdr) && opts.output == "")
     sigset(SIGPIPE, &del_temp);

   immutable(precision[]) rankPhenotype = cast(immutable)setup!(precision)(fileArray, opts);

   if (!opts.run)
     noPerm(fileArray, opts.skip, rankPhenotype);
   else if (!opts.pval && !opts.min && !opts.fdr)
     simplePerm(fileArray, opts, rankPhenotype);
   else if (!opts.min && !opts.fdr)
     pvalPerm(fileArray, opts, rankPhenotype);
   else if (!opts.fdr)
     {
       precision[] minPvalues = minPerm(fileArray, opts, rankPhenotype);
       fileArray[F.out_].close();
       writeFWER(opts, minPvalues);
     }
   else
     fdrCalc(fileArray, opts, rankPhenotype);
 }
