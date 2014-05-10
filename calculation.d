import std.algorithm, std.conv, std.math, std.random, std.range, std.stdio;
import arg_parse : Opts;

extern(C) {
  double gsl_cdf_tdist_P (double x, double nu);
}

class VarianceException : Exception {
  this(string s) {super(s);}
}

void rank(double[] rankArray){

  auto orderIndex = new size_t[rankArray.length];
  makeIndex!("a < b")(rankArray, orderIndex);

  double sumrank = 0.0;
  size_t dupcount = 0;
  double avgrank;

  foreach(i, e; orderIndex)
    {
      sumrank += i;
      dupcount++;
      if (i==(orderIndex.length - 1)
	  || rankArray[e] != rankArray[orderIndex[i + 1]])
	{
	  avgrank = sumrank / dupcount + 1;
	  for (auto j = i - dupcount + 1; j < i + 1; j++)
	    rankArray[orderIndex[j]] = avgrank;
	  sumrank = 0;
	  dupcount = 0;
	}
    }
}

void transform(double[] vector){
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

double[3] correlation(double[] vector1, immutable(double[]) vector2){
  double[3] results;

  results[0] = reduce!((a,b) => a + b[0] * b[1])(0.0, zip(vector1, vector2));
  results[1] = results[0] * sqrt((vector1.length - 2) / (1 - results[0] * results[0]));
  results[2] = gsl_cdf_tdist_P(-fabs(results[1]), vector1.length - 2) * 2;
  return results;
}

double corPvalue(double[] vector1, immutable(double[]) vector2){

  double results = reduce!((a,b) => a + b[0] * b[1])(0.0, zip(vector1, vector2));
  results = results * sqrt((vector1.length - 2) / (1 - results * results));
  results = gsl_cdf_tdist_P(-fabs(results), vector1.length - 2) * 2;
  return results;
}

double[][] getPerm(Opts permOpts, immutable(double[]) vector){
  double[][] outPerm;

  if (permOpts.give_seed)
    rndGen.seed(permOpts.seed);

  outPerm = new double[][permOpts.number];
  for (int i = 0; i < permOpts.number; i++)
    {
      outPerm[][i] = vector.dup;
      randomShuffle(outPerm[][i]);
    }
  return outPerm;
}



