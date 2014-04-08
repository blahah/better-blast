#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/store.h>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/alignment_free.h>
#include <seqan/alignment_free/kmer_functions.h>
#include <iostream>
#include <sstream>

using namespace seqan;

// --------------------------------------------------------------------------
// Class Options
// --------------------------------------------------------------------------

struct Options
{
  std::string query;
  std::string target;
};

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
  // Setup ArgumentParser.
  seqan::ArgumentParser parser("bblast");

  // Short desc, version and date
  setShortDescription(parser, "Prototype for better BLAST");
  setVersion(parser, "0.1");
  setDate(parser, "November 2013");

  // Query
  addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "query"));
  setValidValues(parser, 0, getFileFormatExtensions(seqan::AutoSeqFormat()));
  setHelpText(parser, 0, "Query (multi)FASTA file.");

  // Target
  addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "target", true));
  setValidValues(parser, 1, getFileFormatExtensions(seqan::AutoSeqFormat()));
  setHelpText(parser, 1, "Target multiFASTA file.");

  // Usage
  addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIQUERY\\fP> <\\fITARGET\\fP>");

  // Long description
  addDescription(parser, "Better BLAST (bblast) is a protoype of a massively improved BLAST");
  addDescription(parser, "Currently implemented: nucleotide query -> nucleotide target");
  addDescription(parser, "Input to bblast is a (multi)FASTA query file and a multiFASTA target file");

  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  // Only extract  options if the program will continue after parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res;

  getArgumentValue(options.query, parser, 0);
  getArgumentValue(options.target, parser, 1);
  
  return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function d2CompareCounts()
// --------------------------------------------------------------------------
template <typename TValue>
void
d2CompareCounts(TValue & result,
                String<unsigned> const & kmerCounts1,
                String<unsigned> const & kmerCounts2,
                AFScore<D2> const & score)
{
  typedef typename Iterator<String<unsigned> const>::Type TIteratorInt;

  TIteratorInt it1 = begin(kmerCounts1);
  TIteratorInt it2 = begin(kmerCounts2);

  result = 0;
  for (; it1 < end(kmerCounts1); ++it1) {
    result += value(it1) * value(it2);
    ++it2;
  }
}

// --------------------------------------------------------------------------
// Function computeD2DistanceMatrix()
// --------------------------------------------------------------------------
template <typename TStringSet, typename TValue>
void computeD2DistanceMatrix(Matrix<TValue, 2> & scoreMatrix,
                           TStringSet const & querySet,
			                     TStringSet const & targetSet,
                           AFScore<D2> const & score)
{

  typedef typename Iterator<TStringSet const>::Type               TIteratorSet;
  typedef typename Iterator<StringSet<String<unsigned> > >::Type  TIteratorSetInt;

  unsigned queryNumber = length(querySet);
  unsigned targetNumber = length(targetSet);

  StringSet<String<unsigned> > qKmerCounts;
  StringSet<String<unsigned> > tKmerCounts;

  resize(qKmerCounts, queryNumber);
  resize(tKmerCounts, targetNumber);

  // Count query kmers
  TIteratorSetInt itQKmerCounts = begin(qKmerCounts);
  TIteratorSet itQSeqSet = begin(querySet);

  for (; itQSeqSet < end(querySet); ++itQSeqSet) {
    countKmers(value(itQKmerCounts), value(itQSeqSet), score.kmerSize);
    ++itQKmerCounts;
  }
  // Count target kmers
  TIteratorSetInt itTKmerCounts = begin(tKmerCounts);
  TIteratorSet itTSeqSet = begin(targetSet);

  for (; itTSeqSet < end(targetSet); ++itTSeqSet) {
    countKmers(value(itTKmerCounts), value(itTSeqSet), score.kmerSize);
    ++itTKmerCounts;
  }
  // Resize the scoreMatrix
  setLength(scoreMatrix, 0, queryNumber);
  setLength(scoreMatrix, 1, targetNumber);
  resize(scoreMatrix, (TValue) 0);

  // Calculate all pairwise scores and store them in scoreMatrix
  for (unsigned rowIndex = 0; rowIndex < queryNumber; ++rowIndex) {
    if(score.verbose) {
	    std::cout << "\nSequence number " << rowIndex;
    }
    for (unsigned colIndex = rowIndex; colIndex < targetNumber; ++colIndex) {
      d2CompareCounts(value(scoreMatrix, rowIndex, colIndex), qKmerCounts[rowIndex], tKmerCounts[colIndex], score);
	    value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
    }
  }
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, const char *argv[]) {
  // Parse the command line.
  Options options;

  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

  // If there was an error parsing or built-in argument parser functionality
  // was triggered then we exit the program.  The return code is 1 if there
  // were errors and 0 if there were none.
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  // Load the query sequences
  std::cerr << "Loading query sequences from " << options.query << std::endl;
  seqan::StringSet<seqan::CharString> queryIds;
  seqan::StringSet<seqan::Dna5String> querySeqs;
  seqan::SequenceStream queryStream(options.query.c_str());
  readAll(queryIds, querySeqs, queryStream);

  // Load the target sequences
  std::cerr << "Loading target sequences from " << options.target << std::endl;
  seqan::StringSet<seqan::CharString> targetIds;
  seqan::StringSet<seqan::Dna5String> targetSeqs;
  seqan::SequenceStream targetStream(options.target.c_str());
  readAll(targetIds, targetSeqs, targetStream);
  
  std::cerr << "Sequences loaded" << std::endl;

  // Setup distance matrix
  Matrix<double, 2> distMatrix;
  unsigned kmerSize = 5;
  bool verbose = false;
  AFScore<D2> myScoreD2(kmerSize, verbose);

  // Compute distances
  std::cerr << "Computing D2 distances" << std::endl;
  computeD2DistanceMatrix(distMatrix, querySeqs, targetSeqs, myScoreD2);

  std::cout << distMatrix;
}
