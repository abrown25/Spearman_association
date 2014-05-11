import std.algorithm, std.conv, std.math, std.random, std.range, std.stdio;
import arg_parse : Opts;

extern(C) {
  double gsl_cdf_tdist_P (double x, double nu);
}

class VarianceException : Exception {
  this(string s) {super(s);}
}

ref double[] rank(ref double[] rankArray){

  size_t len = rankArray.length;
  auto orderIndex = new size_t[len];
  makeIndex!("a < b")(rankArray, orderIndex);

  double sumrank = 0.0;
  size_t dupcount = 0;
  double avgrank;

  foreach(i, e; orderIndex)
    {
      sumrank += i;
      dupcount++;
      if (i == (len - 1) || rankArray[e] != rankArray[orderIndex[i + 1]])
	{
	  avgrank = sumrank / dupcount + 1;
	  for (auto j = i - dupcount + 1; j < i + 1; j++)
	    rankArray[orderIndex[j]] = avgrank;
	  sumrank = 0;
	  dupcount = 0;
	}
    }
  return rankArray;
}

void transform(ref double[] vector){
  int n = 0;
  double mean = 0;
  double M2 = 0;
  double delta;

  foreach(ref e; vector)
    {
      n++;
      delta = e - mean;
      mean += delta / n;
      M2 += delta * (e - mean);
    }

  if (M2==0)
    throw new VarianceException("");

  M2 = sqrt(M2);

  foreach(i, ref e; vector)
    e = (e - mean) / M2;

}

double[3] correlation(in double[] vector1, immutable(double[]) vector2){
  double[3] results = 0;

  foreach(i, ref e; vector1)
    results[0] += e * vector2[i];

  results[1] = results[0] * sqrt((vector1.length - 2) / (1 - results[0] * results[0]));
  results[2] = gsl_cdf_tdist_P(-fabs(results[1]), vector1.length - 2) * 2;
  return results;
}

double corPvalue(in double[] vector1, immutable(double[]) vector2){
  double results = 0;
  foreach(i, ref e; vector1)
    results += e * vector2[i];

  results = results * sqrt((vector1.length - 2) / (1 - results * results));
  results = gsl_cdf_tdist_P(-fabs(results), vector1.length - 2) * 2;
  return results;
}

double[][] getPerm(in Opts permOpts, immutable(double[]) vector){
  double[][] outPerm;

  if (permOpts.give_seed)
    rndGen.seed(permOpts.seed);
  outPerm = new double[][permOpts.number];

  foreach(ref e; outPerm)
    {
      e = vector.dup;
      randomShuffle(e);
    }

  return outPerm;
}
