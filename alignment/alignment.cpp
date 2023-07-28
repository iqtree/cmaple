#include "alignment.h"

using namespace std;
using namespace cmaple;

cmaple::Alignment::Alignment():aln_base(new AlignmentBase()) {};

cmaple::Alignment::Alignment(std::istream& aln_stream, const std::string& ref_seq, const InputType format, const SeqType seqtype):aln_base(new AlignmentBase())
{
    // Initialize an alignment instance from the input stream
    read(aln_stream, ref_seq, format, seqtype);
}

cmaple::Alignment::Alignment(const std::string& aln_filename, const std::string& ref_seq, const InputType format, const SeqType seqtype):aln_base(new AlignmentBase())
{
    // Initialize an alignment instance from the input file
    read(aln_filename, ref_seq, format, seqtype);
}

void cmaple::Alignment::read(std::istream& aln_stream, const std::string& ref_seq, const InputType format, const SeqType seqtype)
{
    // Make sure aln_base is initialized
    ASSERT(aln_base && "Null pointer to aln_base");
    
    if (cmaple::verbose_mode >= cmaple::VB_MED)
        std::cout << "Reading an alignment" << std::endl;
    
    // Reset aln_base
    aln_base->reset();
    
    // Set format (is specified)
    if (format != IN_UNKNOWN)
        aln_base->aln_format = format;
    // Otherwise, detect the format from the alignment
    else
    {
        aln_base->aln_format = detectInputFile(aln_stream);
        
        // validate the format
        if (aln_base->aln_format == IN_UNKNOWN)
            outError("Failed to detect the format from the alignment!");
    }
    
    // Set seqtype. If it's unknown (not specified), we'll dectect it later when reading the alignment
    aln_base->setSeqType(seqtype);
    
    // Read the input alignment
    // in FASTA or PHYLIP format
    if (aln_base->aln_format != IN_MAPLE)
        aln_base->readFastaOrPhylip(aln_stream, ref_seq);
    // in MAPLE format
    else
        aln_base->readMaple(aln_stream);
}

void cmaple::Alignment::read(const std::string& aln_filename, const std::string& ref_seq, const InputType format, const SeqType seqtype)
{
    ASSERT(aln_filename.length() && "Please specify an alignnment file");
    
    // Create a stream from the input alignment
    ifstream aln_stream;
    try {
        aln_stream.exceptions(ios::failbit | ios::badbit);
        aln_stream.open(aln_filename);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, aln_filename);
    }
    
    // Initialize an alignment instance from the input stream
    read(aln_stream, ref_seq, format, seqtype);
    
    // close aln_stream
    aln_stream.close();
}

cmaple::Alignment::~Alignment()
{
    if (aln_base)
    {
        delete aln_base;
        aln_base = nullptr;
    }
}

void cmaple::Alignment::write(std::ostream& aln_stream, const InputType format) const
{
    ASSERT(aln_base);
    
    // Handle empty alignment
    if (!aln_base->data.size())
        outError("Alignment is empty. Please call read(...) first!");
        
    // Validate the format
    if (aln_base->aln_format == IN_UNKNOWN)
        outError("Unknown format! Please use IN_MAPLE, IN_FASTA, or IN_PHYLIP");
    
    // Write the alignment
    aln_base->write(aln_stream, format);
}

void cmaple::Alignment::write(const std::string& aln_filename, const InputType format, const bool overwrite) const
{
    // Validate the input
    if (!aln_filename.length())
        outError("Please specify a filename to output the alignment");
    
    // Check whether the output file already exists
    if (!overwrite && fileExists(aln_filename))
        outError("File " + aln_filename + " already exists. Please set overwrite = true to overwrite it.");
    
    // Open a stream to write the output
    std::ofstream aln_stream = ofstream(aln_filename);
    
    // Write alignment to the stream
    write(aln_stream, format);
    
    // Close the stream
    aln_stream.close();
}

std::ostream& cmaple::operator<<(std::ostream& out_stream, const cmaple::Alignment& aln)
{
    // write the alignment to the stream in MAPLE format
    aln.write(out_stream);
    
    // return the stream
    return out_stream;
}

std::istream& cmaple::operator>>(std::istream& in_stream, cmaple::Alignment& aln)
{
    // read the alignment from a stream
    aln.read(in_stream);
    
    // return the stream
    return in_stream;
}
