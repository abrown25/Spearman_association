import std.stdio;
import std.array;
import std.algorithm;
import std.conv;
import std.c.stdlib : exit;
import std.file;
import std.string;
import std.range;
// Usage is vcf_parse DS:PP:GT vcfFile.vcf, options can be altered to change preferential order for extracting fields.
// If you write a new function for a field, add it to the vector below
immutable auto functions = ["PP", "DS", "GT"];

pure nothrow string genFunctionPointer(immutable string[] x){
  string y = "size_t ind;
";
  foreach(ref e; x)
    y ~= "if ((ind = options.countUntil(\"" ~ e ~ "\")) != -1)
    functionArray[ind] = &" ~ e ~ ";
";
  return y;
}

pure nothrow int convInt(const ubyte x)
{
  return x > 57 || x < 48 ? -1 : x - 48;
}

pure nothrow auto convDouble(const char[] x){
  auto y = cast(ubyte[])x;
  if (y.length == 0)
    return -1;

  double value = 0;
  auto z = convInt(y[0]);
  if (z == -1)
    return -1;

  value += z;
  if (y.length == 1)
    return value;
  if (y[1] != cast(ubyte)'.')
    return -1;
  if (y.length > 2)
    {
      double offset = 0.1;
      foreach(ref e; y[2..$])
	{
	  z = convInt(e);
	  if (z == -1)
	    return -1;
	  value += z * offset;
	  offset = offset / 10;
	}
    }
  return value;
}

pure nothrow double GT(const char[] x)
{
  auto y = cast(ubyte[])x;
  if (y.length != 3 || (y[1] != cast(ubyte)'|' && y[1] != cast(ubyte)'/'))
    return -1.0;
  int[2] z = [convInt(y[0]), convInt(y[2])];
  return z[0] == -1 || z[1] == -1 ? -1 : (z[0]!=0) + (z[1]!=0);
}

unittest{
  assert(GT("|||")==-1);
  assert(GT("3\4")==-1);
  assert(GT("0|0")==0);
  assert(GT("3333|4")==-1);
  assert(GT("3|0")==1);
  assert(GT("3/4")==2);
}

pure nothrow auto DS(const char[] x)
{
  return convDouble(x);
}

unittest{
  assert(DS("0.542653")==0.542653);
  assert(DS("40.542653")==-1);
  assert(DS("0l542353")==-1);
  assert(DS("")==-1);
  assert(DS("4")==4);
  assert(DS("0.")==0);
}

pure auto PP(const char[] x)
{
  if (x.count(",")!=2)
    return -1;
  auto  y = x.splitter(',');
  auto  z = y.drop(1).map!convDouble.array;
  return z[0] == -1 || z[1] == -1 ? -1 : z[0] + 2 * z[1];
}

unittest{
  assert(PP("1,3,2")==7);
  assert(PP("1,2.4,3.6")==9.6);
  assert(PP("1,2,A")==-1);
  assert(PP("1,A,3")==-1);
  assert(PP("1")==-1);
}

void main(string[] args)
{
  if (args.length != 3)
    {
      writeln("Need parsing options and input file.");
      exit(0);
    }

  immutable auto options = cast(immutable string[])args[1].split(":")
					        	  .filter!(x => functions.countUntil(x) != -1)
							  .array;

  if (options.length == 0)
    {
      writeln("No suitable options to run.");
      exit(0);
    }

  auto functionArray = new typeof(&PP)[options.length];

  mixin(genFunctionPointer(functions));

  File inFile;
  try{
    inFile = File(args[2]);
  }catch (Exception e){
    stdout.write(e.msg);
    exit(0);
  }

  char[][] infoFields;
  char[][] indVCF;
  long[] mapping;
  double dosage;

  foreach(line; inFile.byLine)
    {
      if (line[0..2]=="##")
	continue;
      else if (line[0..2] =="#C")
	{
	  auto splitLine = line.strip.splitter("\t");
	  stdout.writeln(join(splitLine.take(5), "\t"), "\t", join(splitLine.drop(9), "\t"));
	}
      else
	{
	  auto splitLine = line.strip.splitter("\t");
	  stdout.write(join(splitLine.take(5), "\t"));
	  infoFields = splitLine.drop(8).front.split(":");
	  mapping = options.map!(x => infoFields.countUntil(x)).array;
	  foreach(ref e; splitLine.drop(9))
	    {
	      indVCF = e.split(":");
	      bool done = false;
	      foreach(ref f; zip(mapping, functionArray))
	      	{
	      	  if (f[0] == -1)
		    continue;
		  else
	      	    {
	      	      dosage = f[1](indVCF[f[0]]);
	      	      if (dosage != -1)
	      	      	{
			  stdout.write("\t", dosage);
			  done = true;
	      	      	  break;
	      	      	}
	      	    }
		}
	      if (!done)
		stdout.write("\tNA");
	    }
	  stdout.write("\n");
	}
    }
}
