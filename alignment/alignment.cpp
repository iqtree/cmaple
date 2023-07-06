#include "alignment.h"

using namespace std;
using namespace cmaple;

cmaple::Alignment::Alignment(const std::string& aln_filename, const std::string& ref_seq, const bool overwrite, const std::string& format, const std::string& seqtype):aln_base(nullptr)
{
    ASSERT(aln_filename.length() && "Please specify an alignnment file");
    
    // Init alignment base
    aln_base = new AlignmentBase();
    
    // set the aln_filename
    aln_base->aln_filename = aln_filename;
    
    // if users specify the format -> parse it
    if (format.length())
    {
        aln_base->aln_format = aln_base->getAlignmentFormat(format);
        
        // validate the format
        if (aln_base->aln_format == IN_UNKNOWN)
            outError("Unsupported alignment format " + format + ". Please use MAPLE, FASTA, or PHYLIP");
    }
    // Otherwise, detect the format from the alignment
    else
    {
        aln_base->aln_format = detectInputFile(aln_filename.c_str());
        
        // validate the format
        if (aln_base->aln_format == IN_UNKNOWN)
            outError("Failed to detect the format from the alignment!");
    }
    
    // If users specify the sequence type -> parse it. Otherwise, we'll dectect it later when reading the alignment
    if (seqtype.length())
    {
        // parse the sequence type
        aln_base->setSeqType(aln_base->getSeqType(seqtype));
        
        // validate the sequence type
        if (aln_base->getSeqType() == SEQ_UNKNOWN)
            outError("Unknown sequence type " + seqtype + ", please use DNA or AA");
    }
    
    // if alignment is in PHYLIP or FASTA format -> convert to MAPLE format
    if (aln_base->aln_format != IN_MAPLE)
    {
        // record the starting time
        auto start = getRealTime();
        
        // prepare output (MAPLE) file
        std::string maple_path = aln_filename + ".maple";
        
        aln_base->extractMapleFile(maple_path, ref_seq, overwrite);
        
        // record the end time and show the runtime
        auto end = getRealTime();
        cout << "The input alignment is converted into MAPLE format at " << maple_path << endl;
        cout << " - Converting time: " << end-start << endl;
    }
    // otherwise, alignment is in MAPLE format -> read it
    else
        aln_base->readMapleFile();
}

cmaple::Alignment::~Alignment()
{
    if (aln_base)
    {
        delete aln_base;
        aln_base = nullptr;
    }
}
