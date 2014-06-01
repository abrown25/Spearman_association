
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

import std.algorithm : reduce;
import std.file : exists;
import std.stdio : stdin; 

import arg_parse;
import calculation : rank, transform, VarianceException, covariates;
import run_analysis;
import setup_all : fileSetup, setup;

void main(in string[] args){
  if (args.length == 1)
    giveHelp(helpString);

  string[string] options = getOpts(args[1..$]);
  auto opts = new Opts(options);

  File[3] fileArray;

  scope(failure){
    fileArray[outF].close();
    if (opts.min && (opts.output ~ "temp").exists)
      remove((opts.output ~ "temp"));
   }

  scope(exit){
    fileArray[phenF].close();
    fileArray[genF].close();
  }

  fileSetup(fileArray, opts);

  immutable(double[]) rankPhenotype = cast(immutable)setup(fileArray, opts);

  if (!opts.run)
    noPerm(fileArray, opts.skip, rankPhenotype);
  else if (!opts.pval && !opts.min)
    simplePerm(fileArray, opts, rankPhenotype);
  else if (!opts.min)
    pvalPerm(fileArray, opts, rankPhenotype);
  else
    {
      double[] minPvalues = minPerm(fileArray, opts, rankPhenotype);
      fileArray[outF].close();
      writeFWER(opts, minPvalues);
    }
}
