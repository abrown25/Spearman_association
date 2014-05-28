module calculation;

import std.algorithm : makeIndex;
import std.array : split;
import std.conv : to, ConvException;
import std.math : fabs, sqrt;
import std.random : rndGen, randomShuffle;
import std.range : lockstep;
import std.stdio : File, writeln;
import std.numeric: dotProduct;

import arg_parse : Opts;
import run_analysis : InputException;

extern(C) {
  double gsl_cdf_tdist_P(double x, double nu);
}

extern(C) {
  void regress(size_t nInd, size_t nCov, double *x, double *y, double *rOut);
}

class VarianceException : Exception {
  this(string s) {super(s);}
}

void covariates(string covF, ref double[] phen){
  double[] covOut;
  auto covFile = File(covF);
  size_t nInd = 1;
  covOut ~= 1;

  auto firstLine = split(covFile.readln());
  size_t nCov = firstLine.length;
  covOut ~= to!(double[])(firstLine);

  foreach(line; covFile.byLine())
    {
      covOut ~= 1;
      auto splitLine = split(line);
	{
	  if(splitLine.length != nCov)
	    throw new InputException("Failed to run analysis, covariate measurements missing");
	  else
	    {
	      covOut ~= to!(double[])(splitLine);
	      nInd++;
	    }
	}
    }

  if (nInd != phen.length)
    throw new InputException("Failed to run analysis, covariates file does not have the same number of individuals as phenotype file");

  regress(nInd, nCov + 1, covOut.ptr, phen.dup.ptr, phen.ptr);
}


ref double[] rank(ref double[] rankArray){

  immutable size_t len = rankArray.length;
  auto orderIndex = new size_t[len];
  makeIndex!("a < b")(rankArray, orderIndex);

  double sumrank = 0.0;
  size_t dupcount = 0;
  double avgrank;

  foreach(i, ref e; orderIndex)
    {
      sumrank += i;
      dupcount++;
      if (i == (len - 1) || rankArray[e] != rankArray[orderIndex[i + 1]])
	{
	  avgrank = sumrank / dupcount + 1;
	  for(auto j = i - dupcount + 1; j < i + 1; j++)
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

  foreach(ref e; vector)
    e = (e - mean) / M2;
}

double[3] correlation(in double[] vector1, immutable(double[]) vector2){
  double[3] results;
  results[0] = dotProduct(vector1, vector2);
  results[1] = results[0] * sqrt((vector1.length - 2) / (1 - results[0] * results[0]));
  results[2] = gsl_cdf_tdist_P(-fabs(results[1]), vector1.length - 2) * 2;
  return results;
}

ref double corPvalue(ref double results, in size_t len){
  results = results * sqrt((len - 2) / (1 - results * results));
  results = gsl_cdf_tdist_P(-fabs(results), len - 2) * 2;
  return results;
}

double[] getPerm(in Opts permOpts, immutable(double[]) vector){
  if (permOpts.give_seed)
    rndGen.seed(permOpts.seed);

  double[] outPerm = new double[vector.length * permOpts.number];
  for(auto i = 0; i < permOpts.number; i++)
    {
      outPerm[(i * vector.length)..(i + 1) * vector.length] = vector.dup;
      randomShuffle(outPerm[(i * vector.length)..(i + 1) * vector.length]);
    }
  return outPerm;
}
