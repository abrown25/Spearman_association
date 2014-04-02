import std.algorithm, std.conv, std.math, std.random;
import arg_parse;
extern(C) {
  double gsl_cdf_tdist_P (double x, double nu);
}

double[] rank(ref double[] rankArray){
  ulong[] orderIndex = new ulong[rankArray.length];
  double[] rankIndex = new double[rankArray.length];
  ulong sumrank = 0;
  ulong dupcount = 0;
  double avgrank;

  foreach(i, ref x; orderIndex)
    x = i;

  sort!((a,b) { return rankArray[a] < rankArray[b] ; } )(orderIndex);

  foreach(i, e; orderIndex)
    {
      sumrank += i;
      dupcount++;
      if (i==(orderIndex.length - 1)
	  || rankArray[e] != rankArray[orderIndex[i+1]])
	{
	  avgrank = to!double(sumrank)/dupcount +1;
	  for (ulong j = i - dupcount + 1; j < i + 1; j++)
	    rankIndex[orderIndex[j]] = avgrank;
	  sumrank = 0;
	  dupcount = 0;
	}
    }
  return rankIndex;
}

double[] transform(double[] vector){
  int n = 0;
  double mean = 0;
  double M2 = 0;
  double delta;
  double[] normalised = new double[vector.length];

  foreach(e; vector)
    {
      n++;
      delta = e - mean;
      mean += delta / n;
      M2 += delta * (e - mean);
    }

  M2 = sqrt(M2);

  foreach(i, e; vector)
    normalised[i] = (e - mean) / M2;

  return normalised;
}

double[] correlation(double[] vector1, immutable(double[]) vector2){
  double[] results = [0.0, 0.0, 0.0];

  foreach(i, e; vector1)
    results[0] += e*vector2[i];

  results[1] = results[0] * sqrt((vector1.length - 2) / (1 - results[0]*results[0]));
  results[2] = gsl_cdf_tdist_P(-fabs(results[1]), vector1.length - 2) * 2;
  return results;
}

double corPvalue(double[] vector1, immutable(double[]) vector2){
  double results = 0.0;
  
  foreach(i, e; vector1)
    results += e*vector2[i];

  results = results * sqrt((vector1.length - 2) / (1 - results*results));
  results = gsl_cdf_tdist_P(-fabs(results), vector1.length - 2) * 2;
  return results;
}

double[][] getPerm(PermOpts permOpts, immutable(double[]) vector1){
  double[][] outPerm;

  if (permOpts.give_seed)
    rndGen.seed(permOpts.seed);
  outPerm = new double[][permOpts.number];
  for (int i=0; i<permOpts.number; i++)
    {
      outPerm[][i] = vector1.dup;
      randomShuffle(outPerm[][i]);
    }
  return outPerm;
}
