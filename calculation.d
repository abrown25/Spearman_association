module calculation;

import std.exception : enforce;
import std.math : fabs, sqrt;
import arg_parse : Opts;

enum double EPSILON = 0.00000001; //comparison for X>=Y is done X > Y - epsilon 

version(unittest)
{
  import std.math : approxEqual;
}

class VarianceException : Exception
{
  //thrown if variable is constant
  pure nothrow this(string s) {super(s);}
}

pure nothrow extern(C)
{
  //call GSL to calculate P values from T statistics
  double gsl_cdf_tdist_P(double x, double nu);
}

unittest
{
  // Checks GSL gives right p value for t statistic
  assert(approxEqual(gsl_cdf_tdist_P(-1.6, 7), 0.07681585));
}

pure nothrow extern(C)
{
  //performs regression and returns residuals in rOut 
  void regress(size_t nInd, size_t nCov, double *x, double *y, double *rOut);
}

unittest
{
  //Runs a sample linear regression and checks values against R residuals
  double [] residualsFromR = [-1.00027816411683, -1.72461752433936, -0.275938803894297, 1.51766342141864, 1.48317107093185];

  size_t nInd = 5;
  size_t nCov = 4;
  double[] x = [1, 2, 1, 5, 1, 9, 5, 1, 1, 2, 9, 8, 1, 8, 7, 1, 1, 4, 1, 5];
  double[] y = [7, 2, 2, 8, 5];
  double[5] residuals;

  regress(nInd, nCov, x.ptr, y.ptr, residuals.ptr);

  foreach(i, ref e; residuals)
    assert(approxEqual(e, residualsFromR[i]));
}

void covariates(string covF, ref double[] phen)
{
  //given name of covariate file, and phenotype, produces residuals
  import std.array : split;
  import std.conv : to, ConvException;
  import std.stdio : File;
  import run_analysis : InputException;

  double[] covOut;
  size_t nInd = 1;
  covOut ~= 1;

  auto covFile = File(covF);
  auto firstLine = split(covFile.readln);
  size_t nCov = firstLine.length;
  covOut ~= to!(double[])(firstLine);

  foreach(line; covFile.byLine)
    {
      covOut ~= 1;
      auto splitLine = split(line);
      enforce(splitLine.length == nCov, new InputException("Failed to run analysis, covariate measurements missing"));
      covOut ~= to!(double[])(splitLine);
      nInd++;
    }

  enforce(nInd == phen.length, new InputException("Failed to run analysis, covariates file does not have the same number of individuals as phenotype file"));

  regress(nInd, nCov + 1, covOut.ptr, phen.dup.ptr, phen.ptr);
}

unittest
{
  //Checks residuals from sample files against residuals calculated by R
  double[] residualsFromR = [-0.0988744898987968, 0.32101118013441, 0.598800504240901, 0.798569051911944, 0.00240916094127165, -0.512816562906186, -0.586676086840575, -0.690613979610441, -0.801131138946856, 0.969322360974328];

  double[] phen = [-1.3853088072, -0.785797093643, 1.14540423638, -0.785797093643, 1.03820492508, -1.25652676836, -0.787662180447, -2.05355237841, -0.245457234103, 1.14277217712];

  covariates("cov.txt", phen);

  foreach(i, ref e; phen)
    assert(approxEqual(e, residualsFromR[i]));
}

pure ref T[] rank(T)(ref T[] rankArray)
{
  //ranks array, giving ties mean rank
  import std.algorithm : makeIndex;

  immutable size_t len = rankArray.length;
  auto orderIndex = new size_t[len];
  makeIndex!("a < b")(rankArray, orderIndex);

  T sumrank = 0.0;
  size_t dupcount = 0;
  T avgrank;

  foreach(i, ref e; orderIndex)
    {
      sumrank += i;
      dupcount++;
      if (i == (len - 1) || rankArray[e] != rankArray[orderIndex[i + 1]])
	{
	  avgrank = sumrank / dupcount + 1;
	  foreach(ref j; orderIndex[(i - dupcount + 1) .. (i + 1)])
	    rankArray[j] = avgrank;
	  sumrank = 0;
	  dupcount = 0;
	}
    }
  return rankArray;
}

unittest
{
  //Simple test of ranking with ties
  double[] vector = [10, 9, 2, 9, 3];

  assert(rank!(double)(vector) == [5, 3.5, 1, 3.5, 2]);
}

pure void transform(T)(ref T[] vector)
{
  //transforms array so mean =0 sum of squares = 1
  int n = 0;
  T mean = 0;
  T M2 = 0;
  T delta;

  foreach(ref e; vector)
    {
      n++;
      delta = e - mean;
      mean += delta / n;
      M2 += delta * (e - mean);
    }

  enforce(M2 != 0, new VarianceException(""));

  M2 = sqrt(M2);

  foreach(ref e; vector)
    e = (e - mean) / M2;
}

unittest
{
  //Checks that transform works on randomly generated vector
  import std.algorithm : reduce;
  import std.random : uniform;

  double[] x = new double[10];
  foreach(ref e; x)
    e = uniform(0.0, 10.0);

  transform!(double)(x);
  auto mean = 0.0.reduce!((a, b) => a + b)(x);

  assert(approxEqual(mean, 0.0));
  assert(approxEqual(0.0.reduce!((a, b) => a + (b - mean) * (b - mean))(x), 1));
}

pure nothrow T[3] correlation(T)(in T[] vector1, immutable(T[]) vector2)
{
  //calculates correlation, t stat and p value for two arrays
  import std.numeric: dotProduct;

  T[3] results;
  results[0] = dotProduct(vector1, vector2);
  results[1] = results[0] * sqrt((vector1.length - 2) / (1 - results[0] * results[0]));
  results[2] = gsl_cdf_tdist_P(-fabs(results[1]), vector1.length - 2) * 2;
  return results;
}

unittest
{
  //Check correlation of phenotype with 3rd row genotype against estimates from R
  double[2] corFromR = [-0.2863051, 0.4225695];

  double[] genotype = [0.115, 2, 0.0964, 1, 1, 1, 0, 1, 0, 0.0563];
  double[] tempPhen = [-1.3853088072, -0.785797093643, 1.14540423638, -0.785797093643, 1.03820492508, -1.25652676836, -0.787662180447, -2.05355237841, -0.245457234103, 1.14277217712];

  transform(rank(tempPhen));
  transform(rank(genotype));
  immutable double[] phen = cast(immutable)tempPhen;
  double[3] cor = correlation(genotype, phen);

  assert(approxEqual(cor[0], corFromR[0]));
  assert(approxEqual(cor[2], corFromR[1]));
}

pure nothrow corPvalue(T)(T results, immutable size_t len)
{
  //just returns p value for 2 arrays (only used on permutations)
  results = results * sqrt((len - 2) / (1 - results * results));
  results = gsl_cdf_tdist_P(-fabs(results), len - 2) * 2;
  return results;
}

T[] getPerm(T)(in Opts permOpts, immutable(T[]) vector)
{
  //takes array and generates a continuous array of permuted versions of this array
  import std.random : rndGen, randomShuffle;
  import std.range : chunks, cycle, take;
  import std.array : array;

  if (permOpts.give_seed)
    rndGen.seed(permOpts.seed);

  T[] outPerm = vector.dup
                      .cycle
                      .take(permOpts.number * vector.length)
                      .array;

  foreach(ref perm; chunks(outPerm, vector.length))
    randomShuffle(perm);

  return outPerm;
}

unittest
{
  //Checks permutation p values (option 4,12) for 3rd phenotype against 5th genotype row against values calculated in R based on those permutations
  //Permutations = matrix(c(2, 9, 5, 3, 1, 8, 6, 7, 10, 4, 2, 7, 5, 6, 9, 10, 4, 1, 3, 8, 7, 10, 8, 3, 1, 5, 2, 6, 9, 4, 3, 5, 7, 9, 10, 2, 8, 4, 1, 6), 10, 4)
  import std.numeric: dotProduct;
  double[] permPvalsR = [0.295293228837933, 0.620081751733492, 0.248565693244932, 0.523263392865814];

  double[] genotype = [0, 2, 0, 0, 0, 2, 0.252, 1, 0.018, 0.367];
  double[] tempPhen = [-1.3853088072, -0.785797093643, 1.14540423638, -0.785797093643, 1.03820492508, -1.25652676836, -0.787662180447, -2.05355237841, -0.245457234103, 1.14277217712];

  string[] options = ["dummy", "--perm" , "4,12", "--p" , ""];
  auto testOpts = new Opts(options);

  transform(rank(tempPhen));
  transform(rank(genotype));
  immutable double[] phen = cast(immutable)tempPhen;

  double[] perms = getPerm(testOpts, phen);
  for (auto i = 0; i < 4; i++)
    {
      auto singlePerm = dotProduct(genotype, perms[i * 10 .. (i + 1) * 10]);
      assert(approxEqual(corPvalue!(double)(singlePerm, 10), permPvalsR[i]));
    }
}
