#pragma once
#include "alignmentbase.h"

namespace cmaple
{
    /** Class represents the input alignment */
    class Alignment {
    public:
        /*! \brief Constructor from a stream of an alignment in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq A feference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): "" (auto detect), "MAPLE", "PHYLIP", "FASTA"
         * @param[in] seqtype Data type of sequences (optional): "" (auto detect), "DNA" (nucleotide data), "AA" (amino acid data)
         */
        Alignment(std::istream& aln_stream, const std::string& ref_seq = "", const std::string& format = "", const std::string& seqtype = "");
        
        /*! \brief Constructor from an alignment file in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_filename Name of an alignment file
         * @param[in] ref_seq A reference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): "" (auto detect), "MAPLE", "PHYLIP", "FASTA"
         * @param[in] seqtype Data type of sequences (optional): "" (auto detect), "DNA" (nucleotide data), "AA" (amino acid data)
         */
        Alignment(const std::string& aln_filename, const std::string& ref_seq = "", const std::string& format = "", const std::string& seqtype = "");
        
        /*! \brief Destructor
         */
        ~Alignment();
        
        /** \brief Write the alignment to a stream in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of the output alignment file
         * @param[in] format Format of the output alignment (optional): "MAPLE", "PHYLIP", "FASTA"
         */
        void write(std::ostream& aln_stream, const std::string& format = "MAPLE");
        
        /** \brief Write the alignment to a file in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_filename Name of the output alignment file
         * @param[in] format Format of the output alignment (optional): "MAPLE", "PHYLIP", "FASTA"
         * @param[in] overwrite TRUE to overwrite the existing output file
         */
        void write(const std::string& aln_filename, const std::string& format = "MAPLE", const bool overwrite = false);
        
        // declare Tree as a friend class
        friend class Tree;
        
    private:
        /**
         A base instance of alignment
         */
        AlignmentBase* aln_base;
        
        /*! \brief Initialize an alignment instance
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq A reference sequence
         * @param[in] format Format of the alignment: "" (auto detect), "MAPLE", "PHYLIP", "FASTA"
         * @param[in] seqtype Data type of sequences: "" (auto detect), "DNA" (nucleotide data), "AA" (amino acid data)
         */
        void init(std::istream& aln_stream, const std::string& ref_seq, const std::string& format, const std::string& seqtype);
        
    };
}
