/* The MIT License

Copyright (C) 2013-2014 Genome Research Ltd.
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

import std.string, std.conv, std.stdio, std.algorithm, std.math, std.c.stdlib, std.file;

extern(C) {
  double gsl_cdf_tdist_P (double x, double nu);
}

auto help_string = "Usage: spearman [options]:
    Options:
        -p    :  phenotype file [mandatory]
        -g    :  genotype file [default stdin]
        -pi   :  phenotype IDs are in the first column, if genotype IDs are also present then we check for mismatches
        -gi   :  genotype IDs are in the first row, if phenotype IDs are also present then we check for mismatches
        -pc   :  column for phenotype values, default is 1 if phenotype IDs are not present, 2 otherwise
        -gs   :  column at which genotype values start, preceding columns are printed
        -perm :  calculated permuted p values, one following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed";

// to do: if headers present for both check IDs, skip first few columns of genotype file, handle permutations, write help file, handle errors

string[string] arg_parse(string[] args){
  // options: -p phenotype, -g genotype, -pi phen ids, -gi gen ids, -perm generate permuations -pc phenotype column, -gs genotype skip
  string[string] opts;
  string prefix;

  foreach(i, arg; args){
    if (arg.startsWith("-")){
      prefix = chompPrefix(arg.idup,"-");
      switch (prefix){
      case "p":
	opts["p"] = args[i+1].idup;
	break;
      case "g":
	opts["g"] = args[i+1].idup;
	break;
      case "gi":
	opts["gi"] = "T";
	break;
      case "pi":
	opts["pi"] = "T";
	break;
      case "pc":
	opts["pc"] = args[i+1].idup;
	break;
      case "gs":
	opts["gs"] = args[i+1].idup;
	break;
      default:
	writeln("Unknown option: ", prefix);
	break;
      }
    }
  }
  if ("p" in opts && !opts["p"].exists){
    writeln("Phenotype file missing");
    exit(0);
  }
  if ("g" in opts && !opts["g"].exists){
    writeln("Genotype file missing");
    exit(0);
  }
  return opts;
}

double[] rank(double[] rank_array){

  ulong[] order_index = new ulong[rank_array.length];
  double[] rank_index = new double[rank_array.length];
  ulong sumrank = 0;
  ulong dupcount = 0;
  double avgrank;

  foreach(i, ref x; order_index){
    x = i;
  }

  sort!((a,b) { return rank_array[a] < rank_array[b] ; } )(order_index);
  foreach(i, e; order_index){
    sumrank += i;
    dupcount++;
    if (i==(order_index.length - 1) 
	|| rank_array[e] != rank_array[order_index[i+1]])
      {
      avgrank = to!double(sumrank)/dupcount +1;
      for (ulong j = i - dupcount + 1; j < i + 1; j++){
	rank_index[order_index[j]] = avgrank;
      }
      sumrank = 0;
      dupcount = 0;
    }
  }

  return rank_index;
}

double[] transform(double[] vector){
  int n = 0;
  double mean = 0;
  double M2 = 0;
  double delta;
  double[] normalised = new double[vector.length];

  foreach(e; vector){
    n++;
    delta = e - mean;
    mean += delta / n;
    M2 += delta * (e - mean);
  }

  M2 = sqrt(M2);
  foreach(i, e; vector){
    normalised[i] = (e - mean) / M2;
  }

  return normalised;
}

double[] correlation(double[] vector1, double[] vector2){

  double[] results = [0.0, 0.0, 0.0];
  foreach(i, e; vector1){
    results[0] += e*vector2[i];
  }
  results[1] = results[0] * sqrt((vector1.length - 2) / (1 - results[0]*results[0]));
  results[2] = gsl_cdf_tdist_P(-fabs(results[1]), vector1.length - 2) * 2;
  return results;
}

void main(string[] args){
  string[string] options;
  double[] phenotype;
  double[] rank_phenotype;
  double[] rank_genotype;
  double[] cor = new double[3];
  int skip = 0;
  int phen_column = 0;
  string[] gen_id;
  string[] phen_id;

  auto phen_file = File();
  auto gen_file = File();

  if (args.length == 1){
    writeln(help_string);
    exit(0);
  }
  
  options = arg_parse(args[1..$]);

  writeln(options);

  phen_file = File(options["p"]);
  if ("g" in options)
    gen_file = File(options["g"]);
  else
    gen_file = stdin;

  if ("gi" in options)
    gen_id = split(chomp(gen_file.readln()));

  if ("gs" in options)
    skip = to!int(options["gs"]);

  writeln(gen_id);

  foreach(line; phen_file.byLine()){
    auto phen_line = split(chomp(line));
    phenotype ~= to!double(phen_line[phen_column]);
    // if ("pi" in options)
    //   phen_id ~= phen_line[to!int(options["pi"])-1];
  }

  rank_phenotype = transform(rank(phenotype));

  foreach(line; gen_file.byLine()){
    double genotype[];
    foreach(e; split(line)){
      genotype ~= to!double(e);
	}
    rank_genotype = transform(rank(genotype));
    cor = correlation(rank_genotype, rank_phenotype);
    writeln(join(to!(string[])(cor), "\t"));
  }

}
