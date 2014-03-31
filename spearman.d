import std.string, std.conv, std.stdio, std.algorithm, std.math, std.c.stdlib, std.file;

extern(C) {
  double gsl_cdf_tdist_P (double x, double nu);
}

// to do: if headers present for both check IDs
//skip first few columns of genotype file, handle permutations, figure out immutable status

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
      default:
	writeln("Unknown option");
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

  double sum = 0;
  double var = 0;
  double[] normalised = new double[vector.length];

  foreach(e; vector){
    sum += e;
  }
  sum /= vector.length;

  foreach(i, e ; vector){
    normalised[i] = e - sum;
  }

  foreach(e; normalised){
    var += e*e;
  }
  var = sqrt(var);
  foreach(ref e; normalised){
    e /= var;
  }

  return normalised;
}

double[] correlation(double[] vector1, double[] vector2){

  double[] results = [0.0, 0.0];
  foreach(i, e; vector1){
    results[0] += e*vector2[i];
  }
  results[1] = results[0] * sqrt((vector1.length - 2) / (1 - results[0]*results[0]));
  return results;
}

void main(string[] args){
  string[string] options;
  double[] phenotype;
  double[] rank_phenotype;
  double[] rank_genotype;
  double[] cor = new double[2];
  auto phen_file = File();
  auto gen_file = File();

  if (args.length == 1){
    writeln("All the options");
    exit(0);
  }
  
  options = arg_parse(args[1..$]);

  writeln(options);

  phen_file = File(options["p"]);
  if ("g" in options)
    gen_file = File(options["g"]);
  else
    gen_file = stdin;
    //  auto phen_file = File(options["p"]);
  //if (!("p" in options))
  //input_file = File(options["p"]);
  //   writeln(e, options[e]);

      foreach(line; phen_file.byLine()){
    phenotype ~= to!double(chomp(line));
  }

  rank_phenotype = transform(rank(phenotype));

  foreach(line; File("genotype.txt").byLine()){
    double genotype[];
    foreach(e; split(line)){
      genotype ~= to!double(e);
	}
    rank_genotype = transform(rank(genotype));
    cor = correlation(rank_genotype, rank_phenotype);
    double p_value = gsl_cdf_tdist_P(-fabs(cor[1]), rank_genotype.length - 2) * 2;
    writeln(cor[0], "\t", cor[1], "\t", p_value);
  }

}

