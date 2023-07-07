#pragma once
#include "alignmentbase.h"

namespace cmaple
{
    /** Class presents the input alignment */
    class Alignment {
    public:
        /*! \brief Alignment constructor
         *
         * Alignment constructor from an alignment file in FASTA, PHYLIP, or MAPLE format
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq Reference sequence
         * @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
         * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
         */
        Alignment(std::istream& aln_stream, const std::string& ref_seq = "", const std::string& format = "", const std::string& seqtype = "");
        
        /*! \brief Alignment constructor
         *
         * Alignment constructor from an alignment file in FASTA, PHYLIP, or MAPLE format
         * @param[in] aln_filename Name of an alignment file
         * @param[in] ref_seq Reference sequence
         * @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
         * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
         */
        Alignment(const std::string& aln_filename, const std::string& ref_seq = "", const std::string& format = "", const std::string& seqtype = "");
        
        /*! \brief Alignment destructor
         *
         * Alignment destructor
         */
        ~Alignment();
        
        // TODO: write the alignment in FASTA, PHYLIP, or MAPLE format
        /**
         Write alignment to a stream
         @param[in] ostream A stream of the output alignment file
         @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
         */
        void write(std::ostream& aln_stream, const std::string& format = "MAPLE");
        
        /**
         Write alignment to a file
         @param[in] ostream A stream of the output alignment file
         @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
         @param[in] overwrite TRUE to overwrite the existing output file
         */
        void write(const std::string& aln_filename, const std::string& format = "MAPLE", const bool overwrite = false);
        
        // declare Tree as a friend class
        friend class Tree;
        
    private:
        /**
         A base instance of alignment
         */
        AlignmentBase* aln_base;
        
        /*! \brief Init an alignment instance
         *
         * Init an alignment instance
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq Reference sequence
         * @param[in] format Alignment format: "", "MAPLE", "PHYLIP", "FASTA"
         * @param[in] seqtype Data type of sequences: "", "DNA", "AA"
         */
        void init(std::istream& aln_stream, const std::string& ref_seq, const std::string& format, const std::string& seqtype);
        
    };
}
