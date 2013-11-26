#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/store.h>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

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
  getArgumentValue(options.target, parser, 0);
  
  return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function loadSequences()
// --------------------------------------------------------------------------
template <typename TSeqSet, typename TNameSet>
bool loadSequences(TSeqSet& sequences, 
                    TNameSet& fastaIDs,
                    const char *fileName)
{
  MultiFasta multiFasta;
  if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
  AutoSeqFormat format;
  guessFormat(multiFasta.concat, format); 
  split(multiFasta, format);
  unsigned seqCount = length(multiFasta);
  resize(sequences, seqCount, Exact());
  resize(fastaIDs, seqCount, Exact());
  for(unsigned i = 0; i < seqCount; ++i) 
    {
      assignSeqId(fastaIDs[i], multiFasta[i], format);
      assignSeq(sequences[i], multiFasta[i], format);
    }
  return (seqCount > 0);
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, const char *argv[])
{
  // Parse the command line.
  Options options;

  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

  // If there was an error parsing or built-in argument parser functionality
  // was triggered then we exit the program.  The return code is 1 if there
  // were errors and 0 if there were none.
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  // Load the query and target sequences
  std::cerr << "Loading query sequences from " << options.query << std::endl;
  seqan::StringSet<seqan::CharString> queryIds;
  seqan::StringSet<seqan::Dna5String> querySeqs;
  seqan::SequenceStream queryStream(options.query.c_str());
  readAll(queryIds, querySeqs, queryStream);

  std::cerr << "Loading target sequences from " << options.target << std::endl;
  seqan::StringSet<seqan::CharString> targetIds;
  seqan::StringSet<seqan::Dna5String> targetSeqs;
  seqan::SequenceStream targetStream(options.target.c_str());
  readAll(targetIds, targetSeqs, targetStream);
  
  std::cerr << "Sequences loaded" << std::endl;
}
