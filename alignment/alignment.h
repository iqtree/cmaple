#pragma once
#include "alignmentbase.h"

namespace cmaple
{
    /** Class represents the input alignment */
    class Alignment {
    public:
        
        /*! \brief Default constructor
         */
        Alignment();
        
        /*! \brief Constructor from a stream of an alignment in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq A feference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_UNKNOWN (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_UNKNOWN (auto detection)
         * @throw std::invalid\_argument if the alignment is empty or in an incorrect format
         */
        Alignment(std::istream& aln_stream, const std::string& ref_seq = "", const InputType format = IN_UNKNOWN, const SeqType seqtype = SEQ_UNKNOWN);
        
        /*! \brief Constructor from an alignment file in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_filename Name of an alignment file
         * @param[in] ref_seq A reference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_UNKNOWN (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_UNKNOWN (auto detection)
         * @throw std::invalid\_argument if the alignment is empty or in an incorrect format
         * @throw ios::failure if the alignment file is not found
         */
        Alignment(const std::string& aln_filename, const std::string& ref_seq = "", const InputType format = IN_UNKNOWN, const SeqType seqtype = SEQ_UNKNOWN);
        
        /*! \brief Destructor
         */
        ~Alignment();
        
        /*! \brief Read an alignment from a stream in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq A reference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_UNKNOWN (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_UNKNOWN (auto detection)
         * @throw std::invalid\_argument if the alignment is empty or in an incorrect format
         */
        void read(std::istream& aln_stream, const std::string& ref_seq = "", const InputType format = IN_UNKNOWN, const SeqType seqtype = SEQ_UNKNOWN);
        
        /*! \brief Read an alignment from a file in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_filename Name of an alignment file
         * @param[in] ref_seq A reference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_UNKNOWN (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_UNKNOWN (auto detection)
         * @throw std::invalid\_argument if the alignment is empty or in an incorrect format
         * @throw ios::failure if the alignment file is not found
         */
        void read(const std::string& aln_filename, const std::string& ref_seq = "", const InputType format = IN_UNKNOWN, const SeqType seqtype = SEQ_UNKNOWN);
        
        /** \brief Write the alignment to a stream in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of the output alignment file
         * @param[in] format Format of the output alignment (optional): IN_MAPLE, IN_FASTA, or IN_PHYLIP
         * @throw std::invalid\_argument if format is unknown
         * @throw std::logic\_error if the alignment is empty (i.e., nothing to write)
         */
        void write(std::ostream& aln_stream, const InputType format = IN_MAPLE) const;
        
        /** \brief Write the alignment to a file in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_filename Name of the output alignment file
         * @param[in] format Format of the output alignment (optional): IN_MAPLE, IN_FASTA, or IN_PHYLIP
         * @param[in] overwrite TRUE to overwrite the existing output file (optional)
         * @throw std::invalid\_argument if any of the following situations occur.
         * - aln_filename is empty
         * - format is unknown
         *
         * @throw std::logic\_error if the alignment is empty (i.e., nothing to write)
         * @throw ios::failure if aln\_filename already exists and overwrite = false
         */
        void write(const std::string& aln_filename, const InputType format = IN_MAPLE, const bool overwrite = false) const;
        
        // declare Tree as a friend class
        friend class Tree;
        
    private:
        /**
         A base instance of alignment
         */
        AlignmentBase* aln_base;
    };

    /** \brief Customized << operator to output the alignment to a stream in MAPLE format
     */
    std::ostream& operator<<(std::ostream& out_stream, const cmaple::Alignment& aln);

    /** \brief Customized >> operator to read the alignment from a stream
     */
    std::istream& operator>>(std::istream& in_stream, cmaple::Alignment& aln);
}
