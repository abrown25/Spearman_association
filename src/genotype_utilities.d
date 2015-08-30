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

immutable auto functions = ["PP", "GP", "DS", "GT"];

alias GP = PP;

pure nothrow string genFunctionPointer(const string[] x)
{
  string y = "size_t ind;
";
  foreach (ref e; x)
    y ~= "if ((ind = options.countUntil(\"" ~ e ~ "\")) != -1)
    functionArray[ind] = &" ~ e ~ ";
";
  return y;
}

static const string helpString = "
genotype_utilities: utilities for processing genotype data files. Contains routines to:

Convert VCF and SNPTEST format data to dosage files.
To reorder and match genotype IDs to a given sample ID file.

Output is sent to stdout.

Usage: genotype_utilities <command> <options>:

Commands:

     --match:           Rearranges genotype file to match file with sample (usually phenotype) IDS. Any individuals only present in phenotype file are ignored. Requires either 2 options, depending on whether the input file is given from the stdin. For example, if IDs are in ID_FILE, genotypes are in GEN_FILE, and 3 columns are used for identifying SNPS, then the following command would work: cat GEN_FILE | ./genotype_utilities --match ID_FILE 3.

       First option:    ID file containing phenotype IDs.
       Second option:   File containing genotype data (Can also be taken from the stdin.).
       Third option:    Number of columns containing SNP identification information (Second option if IDs are given by stdin). For example, this is 9 for vcf files. These columns will be written to stdout without parsing.

     --plink:           Produces dosage file from plink bed file. Takes one option, INPUT, which should be the root of the plink file formats (i.e. INPUT.fam, INPUT.bim and INPUT.bed are the names of the files). Output to stdout is a SNP major dosage file, with extra columns for chromosome, rs ID, base pair location and reference and alternative alleles.

     --vcf:             Reports dosage from vcf files, takes two options by place.

       First option:    FORMAT FIELD, possible values of FORMAT column, separated by :, in order of preference. Available values are currently: GT (genotype threshold), PP (or equivalently GP, posterior probability), and DS (dosage). For example, GT:PP:DS means calculate dosage on GT if available and well-formed, then PP, then DS (reporting NA if none succeed).
       Second option:   VCF FILE NAME. Name of vcf file.
";

@safe pure nothrow int convInt(const ubyte x)
{
  return x > 57 || x < 48 ? -1 : x - 48;
}

pure auto convDouble(const char[] x)
{
  try
  {
    return to!double(x);
  }
  catch (ConvException e)
  {
    return -1.0;
  }
  //     auto y = cast(ubyte[]) x;
  //     if (y.length == 0)
  //         return -1;

  //     double value = 0;
  //     auto z = convInt(y[0]);
  //     if (z == -1)
  //         return -1;

  //     value += z;
  //     if (y.length == 1)
  //         return value;
  //     if (y[1] != cast(ubyte) '.')
  //         return -1;
  //     if (y.length > 2)
  //     {
  //         double offset = 0.1;
  //         foreach (ref e; y[2 .. $])
  //         {
  //             z = convInt(e);
  //             if (z == -1)
  //                 return -1;
  //             value += z * offset;
  //             offset = offset / 10;
  //         }
  //     }
  //     return value;
  // }
}

@safe pure auto GT(const char[] x)
{
  auto y = x.to!(ubyte[]);
  if (y.length != 3 || (y[1] != '|'.to!ubyte && y[1] != '/'.to!ubyte))
    return -1.0;
  int[2] z = [convInt(y[0]), convInt(y[2])];
  return z[0] == -1 || z[1] == -1 ? -1 : (z[0] != 0) + (z[1] != 0);
}

unittest
{
  assert(GT("|||") == -1);
  assert(GT("3\4") == -1);
  assert(GT("0|0") == 0);
  assert(GT("3333|4") == -1);
  assert(GT("3|0") == 1);
  assert(GT("3/4") == 2);
}

pure auto DS(const char[] x)
{
  return convDouble(x);
}

unittest
{
  assert(DS("0.542653") == 0.542653);
  assert(DS("40.542653") == -1);
  assert(DS("0l542353") == -1);
  assert(DS("") == -1);
  assert(DS("4") == 4);
  assert(DS("0.") == 0);
}

pure auto PP(const char[] x)
{
  if (x.count(",") != 2)
    return -1;
  auto y = x.splitter(',');
  auto z = y.drop(1).map!convDouble.array;
  return z[0] == -1 || z[1] == -1 ? -1 : z[0] + 2 * z[1];
}

unittest
{
  assert(PP("1,3,2") == 7);
  assert(PP("1,2.4,3.6") == 9.6);
  assert(PP("1,2,A") == -1);
  assert(PP("1,A,3") == -1);
  assert(PP("1") == -1);
}

void main(string[] args)
{
  if (args.length < 2)
  {
    writeln(helpString);
    exit(0);
  }
  else if (args[1] == "--help")
  {
    writeln(helpString);
    exit(0);
  }
  else if (args[1] == "--vcf")
    vcfFile(args[2 .. $]);
  // else if (args[1] == "snptest")
  //     snpTest(args[2 .. $]);
  else if (args[1] == "--match")
    matchIds(args[2 .. $]);
  else if (args[1] == "--plink")
    plinkConvert(args[2 .. $]);
  else
  {
    writeln(helpString);
    exit(0);
  }
}

void vcfFile(string[] args)
{
  if (args.length == 0 || args.length > 2)
  {
    stderr.writeln(
      "Need parsing options and (optional) input file. Run genotype_utilities --help for help file.");
    exit(0);
  }

  const auto options = args[0].split(":").filter!(x => functions.countUntil(x) != -1).array;

  if (options.length == 0)
  {
    stderr.writeln("No suitable options to run. Run genotype_utilities --help for help file.");
    exit(0);
  }

  auto functionArray = new typeof(&PP)[options.length];

  mixin(genFunctionPointer(functions));

  File inFile;
  try
  {
    if (args.length == 2)
      inFile = File(args[1]);
    else
      inFile = stdin;
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }

  char[][] infoFields;
  char[][] indVCF;
  long[] mapping;
  double dosage;

  foreach (line; inFile.byLine)
  {
    if (line[0 .. 2] == "##")
      continue;
    else if (line[0 .. 2] == "#C")
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
      foreach (ref e; splitLine.drop(9))
      {
        indVCF = e.split(":");
        bool done = false;
        foreach (ref f; zip(mapping, functionArray))
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

// void snpTest(string[] args)
// {
//     File sampleFile, inFile;

//     auto outFile = stdout;

//     try
//     {
//         sampleFile = File(args[0]);
//         inFile = File(args[1]);
//     }
//     catch (Exception e)
//     {
//         writeln(e.msg);
//         exit(0);
//     }

//     auto ids = sampleFile.byLine.drop(2).map!(a => a.split.front).to!string.array;

//     stdout.writeln("#CHR\tLOC\tREF\tALT\t", ids.to!(string[]).join("\t"));

//     foreach (ref line; inFile.byLine)
//     {
//         auto splitLine = line.split;
//         write(splitLine[1].split("-").join("\t"), "\t", splitLine[3], "\t", splitLine[
//             4
//         ], "\t");
//         iota(6, splitLine.length, 3).map!(a => convDouble(
//             splitLine[a]) + 2 * convDouble(splitLine[a + 1])).array.to!(string[]).join(
//             "\t").writeln;
//     }
// }

void matchIds(string[] args)
{
  File idFile, inFile;
  auto outFile = stdout;
  long skip;

  try
  {
    idFile = File(args[0]);
    skip = to!int(args[args.length - 1]);
    if (args.length == 3)
    {
      inFile = File(args[1]);
    }
    else if (args.length == 2)
    {
      inFile = stdin;
    }
    else
    {
      stderr.writeln("Either pass 2 or three arguments to --match");
      exit(0);
    }
  }
  catch (Exception e)
  {
    stderr.write(e.msg);
    exit(0);
  }

  try
  {
  }
  catch (Exception e)
  {
    stderr.write(e.msg);
    exit(0);
  }

  auto ids = idFile.byLine.map!(to!string).array;

  auto line = inFile.readln.split;

  auto places = ids.map!(a => line.countUntil(a)).filter!(a => a != -1).array;
  places = iota(skip).array ~ places;

  line.indexed(places).join("\t").writeln;

  inFile.byLine.map!(a => a.split.indexed(places).joiner("\t")).joiner("\n").writeln;
}

void plinkConvert(string[] args)
{

  if (args.length != 1)
  {
    writeln("Pass only one option to --plink");
    exit(0);
  }

  string[] outputString = ["2", "NA", "1", "0"];

  string input = args[0];

  File bimFile;
  File famFile;

  class InputException : Exception
  {
    //thrown if fam file has too few fields
    pure nothrow this(string s)
    {
      super(s);
    }
  }

  try
  {
    famFile = File(input ~ ".fam");
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }

  string[] id;
  string getID(string line)
  {
    auto splitLine = line.splitter;
    splitLine.popFront;
    if (splitLine.empty)
      throw new InputException("");
    else
      return splitLine.front;
  }

  try
  {
    id = famFile.byLine.map!(a => getID(a.idup)).array;
  }
  catch (InputException e)
  {
    stderr.writeln("Lines in ", input, ".fam have fewer than 2 fields.");
    exit(0);
  }

  size_t nInd = id.length;

  size_t nSnps = 0;

  try
  {
    bimFile = File(input ~ ".bim");
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }

  foreach (line; bimFile.byLine)
  {
    if (line.split.length < 6)
    {
      stderr.writeln("Lines in ", input, ".bim have fewer than 6 fields.");
      exit(0);
    }
    nSnps++;
  }

  bimFile.seek(0);

  size_t nBytes = nInd >> 2;
  size_t offset = nInd & 3;

  if (offset != 0)
    nBytes++;

  ubyte[] bytes;
  try
  {
    bytes = cast(ubyte[]) read(input ~ ".bed");
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }

  if (bytes[0] != 108 || bytes[1] != 27 || bytes[2] != 1)
  {
    stderr.writeln("The magic numbers imply that ", input, ".bed is not a valid bed file.");
    exit(0);
  }
  if ((nBytes * nSnps + 3) != bytes.length)
  {
    stderr.writeln(
      "Number of SNPs in bed file does not match the numbers of individuals and SNPs in bim and fam file.");
    exit(0);
  }

  writeln("CHROM\tID\tPOS\tREF\tALT\t", id.joiner("\t"));

  auto writeByte(ubyte inByte)
  {
    return [outputString[(inByte & 0x03)], outputString[(inByte & 0x0C) >> 2],
      outputString[(inByte & 0x30) >> 4], outputString[(inByte & 0xC0) >> 6]].joiner("\t");
  }

  foreach (snpID, currSNP; zip(bimFile.byLine, bytes[3 .. $].chunks(nBytes)))
  {
    snpID.split.indexed([0, 1, 3, 4, 5]).joiner("\t").write;
    stdout.write("\t");
    if (offset == 0)
      currSNP.map!(writeByte).joiner("\t").writeln;
    else
    {
      currSNP[0 .. ($ - 1)].map!(writeByte).joiner("\t").write;
      stdout.write("\t", outputString[(currSNP[$ - 1] & 0x03)]);
      if (offset > 1)
        stdout.write("\t", outputString[(currSNP[$ - 1] & 0x0C) >> 2]);
      if (offset > 2)
        stdout.write("\t", outputString[(currSNP[$ - 1] & 0x30) >> 4]);
      writeln();
    }
  }

}
