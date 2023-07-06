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
         * @param[in] aln_filename Name of an alignment file
         * @param[in] ref_seq Reference sequence
         * @param[in] overwrite TRUE to overwrite the existing output file
         * @param[in] format Alignment format (optional): "", "MAPLE", "PHYLIP", "FASTA"
         * @param[in] seqtype Data type of sequences (optional): "", "DNA", "AA"
         */
        Alignment(const std::string& aln_filename, const std::string& ref_seq = "", const bool overwrite = false, const std::string& format = "", const std::string& seqtype = "");
        
        /*! \brief Alignment destructor
         *
         * Alignment destructor
         */
        ~Alignment();
        
        // TODO: write the alignment in FASTA, PHYLIP, or MAPLE format
        
    private:
        /**
         A base instance of alignment
         */
        AlignmentBase* aln_base;
        
    };
}
