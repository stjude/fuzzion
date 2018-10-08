//------------------------------------------------------------------------------------
//
// fuzzion.cpp - program to find the reads in a BAM file containing two target
//               sequences, or containing one target sequence and not containing a
//               second target sequence; each sequence is matched approximately,
//               allowing a limited number of substitutions
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2018 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include <sstream>
#include "api/BamReader.h"

const std::string VERSION = "fuzzion 2.0";

const int MIN_TARGET_LENGTH = 8; // a target sequence must be at least this long

const int DEFAULT_MAXSUB = 2;    // default maximum substitutions allowed
int maxsub = DEFAULT_MAXSUB;     // maximum substitutions allowed when matching

std::string bam_filename = "";   // name of BAM file (specified on command line)

typedef std::vector<std::string> StringVector;

//------------------------------------------------------------------------------------

class Target // represents one or more target sequences
{
public:
   Target(const std::string& targetString);

   ~Target();

   std::string reverseComplement() const;

   bool findLeftmost (const char *readseq, int readseqlen, int rightpad,
                      int& matchIndex, int& matchStart) const;
   bool findRightmost(const char *readseq, int readseqlen, int leftpad,
                      int& matchIndex, int& matchStart) const;

   bool   want;      // true if we want to find any one of these target sequences
   int    minseqlen; // length of shortest target sequence in this set
   int    maxseqlen; // length of longest  target sequence in this set
   int    seqcount;  // number of target sequences in this set
   int   *seqlen;    // array containing target sequence lengths
   char **seq;       // array containing target sequences
};

//------------------------------------------------------------------------------------

class TargetPair // represents a labeled pair of Target objects
{
public:
   TargetPair(const std::string& inLabel, const std::string& leftTargetString,
              const std::string& rightTargetString);

   ~TargetPair() { delete left; delete right; }

   TargetPair *createReverseComplement() const;

   void findMatch (const std::string& readName, const std::string& readString) const;

   void writeMatch(const std::string& readName, const std::string& readString,
                   int leftIndex,  int leftStart,
		   int rightIndex, int rightStart) const;

   std::string label;
   Target *left, *right;
};

std::vector<TargetPair *> targetPair;
int numTargetPairs;

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stdout

void showUsage(const char *progname)
{
   std::cout << VERSION << std::endl << std::endl;

   std::cout << "Usage: " << progname
             << " [-maxsub=N]"
             << " bam_file"
             << " < target_sequences"
             << " > matching_reads"
             << std::endl << std::endl;

   std::cout << "  -maxsub=N  maximum substitutions allowed, default is "
             << DEFAULT_MAXSUB << std::endl;
}

//------------------------------------------------------------------------------------
// parseArgs() parses the command-line arguments and returns true if all are valid

bool parseArgs(int argc, char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      std::string arg = argv[i];
      int arglen = arg.length();
      if (arglen == 0)
         continue;

      if (arg[0] == '-') // found an option
         if (arglen > 8 && arg.substr(1, 7) == "maxsub=")
	 {
            std::string s = arg.substr(8);
	    std::stringstream stream(s);
	    stream >> maxsub;
	    if (maxsub < 0)
               return false;
	 }
         else
            return false; // unrecognized option
      else
         if (bam_filename == "")
            bam_filename = arg;
         else
            return false; // extraneous argument
   }

   if (bam_filename == "")
      return false; // missing argument

   return true; // all command-line arguments are valid
}

//------------------------------------------------------------------------------------
// toupperSequence() converts all letters in a sequence to uppercase

std::string toupperSequence(const std::string& sequence)
{
   std::string upperSeq = "";

   int len = sequence.length();

   for (int i = 0; i < len; i++)
   {
      char ch   = std::toupper(sequence[i]);
      upperSeq += ch;
   }

   return upperSeq;
}

//------------------------------------------------------------------------------------
// getDelimitedStrings() extracts delimited string values from a string and appends
// them to a string vector

void getDelimitedStrings(const std::string& s, char delimiter, StringVector& v)
{
   int i = 0, len = s.length();

   do
   {
      std::string value = "";

      while (i < len && s[i] != delimiter)
	 value += s[i++];

      v.push_back(value);
   }
   while (i++ < len);
}

//------------------------------------------------------------------------------------
// isACGT() returns true if the given character is A, C, G or T

bool isACGT(char ch)
{
   ch = std::toupper(ch);

   return (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T');
}

//------------------------------------------------------------------------------------
// isAllACGT() returns true if all characters in the sequence are A, C, G or T

bool isAllACGT(const std::string& sequence)
{
   int len = sequence.length();

   for (int i = 0; i < len; i++)
      if (!isACGT(sequence[i]))
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// Target::Target() parses the given string to obtain one or more target sequences and
// saves them in the new object it is constructing

Target::Target(const std::string& targetString)
{
   std::string s = toupperSequence(targetString);

   if (s.length() > 0 && s[0] == '-')
   {
      want = false;
      s = s.substr(1);
   }
   else
      want = true;

   StringVector t;
   getDelimitedStrings(s, '|', t);

   seqcount  = t.size();
   minseqlen = t[0].length();
   maxseqlen = t[0].length();

   for (int i = 1; i < seqcount; i++)
      if (t[i].length() < minseqlen)
         minseqlen = t[i].length();
      else if (t[i].length() > maxseqlen)
         maxseqlen = t[i].length();

   if (minseqlen < MIN_TARGET_LENGTH)
      throw std::runtime_error("invalid sequence length in " + targetString);

   for (int i = 0; i < seqcount; i++)
      if (!isAllACGT(t[i]))
         throw std::runtime_error("invalid character in " + t[i]);

   seqlen = new int[seqcount];
   seq    = new char *[seqcount];

   for (int i = 0; i < seqcount; i++)
   {
      seqlen[i] = t[i].length();
      seq[i]    = new char[seqlen[i] + 1];
      std::strcpy(seq[i], t[i].c_str());
   }
}

//------------------------------------------------------------------------------------
// Target::~Target() de-allocates the arrays allocated by the constructor

Target::~Target()
{
   for (int i = 0; i < seqcount; i++)
      delete[] seq[i];

   delete[] seq;
   delete[] seqlen;
}

//------------------------------------------------------------------------------------
// reverseSequence() reverses the order of the characters in a sequence

std::string reverseSequence(const std::string& sequence)
{
   std::string reverseSeq = "";

   int len = sequence.length();

   for (int i = len - 1; i >= 0; i--)
      reverseSeq += sequence[i];

   return reverseSeq;
}

//------------------------------------------------------------------------------------
// invertSequence() swaps A and T, and C and G in the sequence

std::string invertSequence(const std::string& sequence)
{
   std::string inverseSeq = "";

   int len = sequence.length();

   for (int i = 0; i < len; i++)
   {
      char ch;

      switch (sequence[i])
      {
         case 'A': ch = 'T'; break;
	 case 'a': ch = 't'; break;

	 case 'T': ch = 'A'; break;
	 case 't': ch = 'a'; break;

         case 'C': ch = 'G'; break;
	 case 'c': ch = 'g'; break;

	 case 'G': ch = 'C'; break;
	 case 'g': ch = 'c'; break;

	 default : ch = sequence[i];
      }

      inverseSeq += ch;
   }

   return inverseSeq;
}

//------------------------------------------------------------------------------------
// Target::reverseComplement() returns a string representing the reverse complement of
// the target sequence(s)

std::string Target::reverseComplement() const
{
   std::string revcomp = (want ? "" : "-");

   for (int i = 0; i < seqcount; i++)
   {
      std::string temp = seq[i];
      temp = reverseSequence(temp);
      temp = invertSequence(temp);

      revcomp += (i > 0 ? "|" + temp : temp);
   }

   return revcomp;
}

//------------------------------------------------------------------------------------
// isMatch() performs a fuzzy match of two sequences; it returns true if the target
// sequence matches the read sequence with no more than maxsub substitutions

inline bool isMatch(const char *readseq, const char *target, int targetlen)
{
   int numsubs = 0;

   for (int i = 0; i < targetlen; i++)
      if (readseq[i] != target[i] && ++numsubs > maxsub)
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// Target::findLeftmost() identifies the target sequence that matches the given read
// sequence with its last base farthest left in the read sequence (leaving the most
// opportunity to find a match to the right); if found, true is returned, matchIndex
// is set to the subscript identifying the matching target sequence, and matchStart is
// set to the start index of the match within the read sequence

bool Target::findLeftmost(const char *readseq, int readseqlen, int rightpad,
                          int& matchIndex, int& matchStart) const
{
   matchIndex = -1; // no match found yet

   int lastMatchEnd = readseqlen - rightpad; // exclusive end point

   for (int i = 0; i < seqcount; i++)
   {
      int lastStart = lastMatchEnd - seqlen[i];

      for (int start = 0; start <= lastStart; start++)
         if (isMatch(&readseq[start], seq[i], seqlen[i]))
	 {
            matchIndex   = i;
	    matchStart   = start;
	    lastMatchEnd = start + seqlen[i] - 1;
	    break;
	 }
   }

   return (matchIndex >= 0);
}

//------------------------------------------------------------------------------------
// Target::findRightmost() identifies the target sequence that matches the given read
// sequence with its first base farthest right in the read sequence (leaving the most
// opportunity to find a match to the left); if found, true is returned, matchIndex
// is set to the subscript identifying the matching target sequence, and matchStart is
// set to the start index of the match within the read sequence

bool Target::findRightmost(const char *readseq, int readseqlen, int leftpad,
                           int& matchIndex, int& matchStart) const
{
   matchIndex = -1; // no match found yet

   int lastMatchStart = leftpad; // inclusive starting point

   for (int i = 0; i < seqcount; i++)
   {
      int firstStart = readseqlen - seqlen[i];

      for (int start = firstStart; start >= lastMatchStart; start--)
         if (isMatch(&readseq[start], seq[i], seqlen[i]))
	 {
            matchIndex     = i;
	    matchStart     = start;
	    lastMatchStart = start + 1;
	    break;
	 }
   }

   return (matchIndex >= 0);
}

//------------------------------------------------------------------------------------
// TargetPair::TargetPair() gets the left and right target sequences from the input
// strings; an exception is thrown if there is something wrong

TargetPair::TargetPair(const std::string& inLabel,
                       const std::string& leftTargetString,
                       const std::string& rightTargetString)
   : label(inLabel), left(new Target(leftTargetString)),
     right(new Target(rightTargetString))
{
   if (label.length() == 0)
      throw std::runtime_error("missing label before " + leftTargetString);

   if (!left->want && !right->want)
      throw std::runtime_error("double negative specified for " + label);
}

//------------------------------------------------------------------------------------
// TargetPair::createReverseComplement() returns a TargetPair object that represents
// the reverse complement of this one

TargetPair *TargetPair::createReverseComplement() const
{
   std::string leftTargetString  = right->reverseComplement();
   std::string rightTargetString = left ->reverseComplement();

   return new TargetPair(label, leftTargetString, rightTargetString);
}

//------------------------------------------------------------------------------------
// TargetPair::findMatch() determines whether this target pair can be found in the
// given read sequence; if so, the read sequence is written with the matches
// highlighted

void TargetPair::findMatch(const std::string& readName,
                           const std::string& readString) const
{
   const char *readseq = readString.c_str();
   int readseqlen      = readString.length();

   int leftIndex, leftStart, rightIndex, rightStart;

   if (left->want &&
       left->findLeftmost(readseq, readseqlen,
                          (right->want ? right->minseqlen : right->maxseqlen),
                          leftIndex, leftStart) &&
       right->findRightmost(readseq, readseqlen, leftStart + left->seqlen[leftIndex],
                            rightIndex, rightStart) == right->want ||
       !left->want &&
       right->findRightmost(readseq, readseqlen, left->maxseqlen,
                            rightIndex, rightStart) &&
       !left->findLeftmost (readseq, readseqlen, readseqlen - rightStart,
                           leftIndex, leftStart))
   {
      writeMatch(readName, readString, leftIndex, leftStart, rightIndex, rightStart);
   }
}

//------------------------------------------------------------------------------------
// highlight() highlights a match as it is written to stdout

void highlight(std::string readString, const char *target)
{
   std::cout << "[";

   int len = readString.length();

   for (int i = 0; i < len; i++)
      std::cout << (char)(readString[i] == target[i] ?
                          readString[i] : std::tolower(readString[i]));

   std::cout << "]";
}

//------------------------------------------------------------------------------------
// TargetPair::writeMatch() writes a hit to stdout

void TargetPair::writeMatch(const std::string& readName,
                            const std::string& readString,
                            int leftIndex,  int leftStart,
			    int rightIndex, int rightStart) const
{
   std::cout << readName << "\t";

   int initialBases = (left->want ? leftStart : rightStart);

   if (initialBases > 0)
      std::cout << readString.substr(0, initialBases);

   if (left->want)
   {
      int leftlen = left->seqlen[leftIndex];
      
      highlight(readString.substr(leftStart, leftlen), left->seq[leftIndex]);

      int nextlen =
         (right->want ? rightStart : readString.length()) - leftStart - leftlen;

      if (nextlen > 0)
         std::cout << readString.substr(leftStart + leftlen, nextlen);
   }

   if (right->want)
   {
      int rightlen = right->seqlen[rightIndex];

      highlight(readString.substr(rightStart, rightlen), right->seq[rightIndex]);

      int nextlen = readString.length() - rightStart - rightlen;

      if (nextlen > 0)
         std::cout << readString.substr(rightStart + rightlen, nextlen);
   }

   std::cout << "\t" << label << "\n";
}

//------------------------------------------------------------------------------------
// readTargetPairs() reads a list of target pairs from stdin and stores each pair and
// its reverse complement in a vector of target pairs

void readTargetPairs()
{
   std::string line;

   while (std::getline(std::cin, line))
   {
      StringVector column;
      getDelimitedStrings(line, '\t', column);

      if (column.size() != 3)
         throw std::runtime_error("unexpected #columns in " + line);

      TargetPair *tp = new TargetPair(column[0], column[1], column[2]);

      targetPair.push_back(tp);
      targetPair.push_back(tp->createReverseComplement());
   }

   numTargetPairs = targetPair.size();

   if (numTargetPairs == 0)
      throw std::runtime_error("no input targets");
}

//------------------------------------------------------------------------------------
// readBamFile() reads a BAM file and writes hits to stdout

void readBamFile()
{
   BamTools::BamReader bamReader;

   if (!bamReader.Open(bam_filename))
      throw std::runtime_error("unable to open " + bam_filename);

   readTargetPairs();

   BamTools::BamAlignment alignment;

   while (bamReader.GetNextAlignment(alignment))
      for (int i = 0; i < numTargetPairs; i++)
         targetPair[i]->findMatch(alignment.Name, alignment.QueryBases);

   bamReader.Close();
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   if (!parseArgs(argc, argv))
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      readBamFile();
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
