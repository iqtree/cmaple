#include "sequence.h"
#include "../utils/timeutil.h"

#ifndef CMAPLE_ALIGNMENT_H
#define CMAPLE_ALIGNMENT_H

namespace cmaple
{
    /** Class presents the input alignment */
    class Alignment {
    public:
        /*!
            Alignment format
         */
        enum InputType {
            IN_FASTA, /*!< FASTA format */
            IN_PHYLIP, /*!< PHYLIP format */
            IN_MAPLE, /*!< [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format */
            IN_AUTO, /*!< Auto detect */
            IN_UNKNOWN, /*!< Unknown format */
        };
        
        // ----------------- BEGIN OF PUBLIC APIs ------------------------------------ //
        /*! \brief Default constructor
         */
        Alignment();
        
        /*! \brief Constructor from a stream of an alignment in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq A feference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_AUTO (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_AUTO (auto detection)
         * @throw std::invalid\_argument if any of the following situations occur.
         * - the alignment is empty or in an incorrect format
         * - the sequences contains invalid states
         */
        Alignment(std::istream& aln_stream, const std::string& ref_seq = "", const InputType format = IN_AUTO, const cmaple::SeqRegion::SeqType seqtype = cmaple::SeqRegion::SEQ_AUTO);
        
        /*! \brief Constructor from an alignment file in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_filename Name of an alignment file
         * @param[in] ref_seq A reference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_AUTO (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_AUTO (auto detection)
         * @throw std::invalid\_argument if any of the following situations occur.
         * - the alignment is empty or in an incorrect format
         * - the sequences contains invalid states
         *
         * @throw ios::failure if the alignment file is not found
         */
        Alignment(const std::string& aln_filename, const std::string& ref_seq = "", const InputType format = IN_AUTO, const cmaple::SeqRegion::SeqType seqtype = cmaple::SeqRegion::SEQ_AUTO);
        
        /*! \brief Destructor
         */
        ~Alignment();
        
        /*! \brief Read an alignment from a stream in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of an alignment file
         * @param[in] ref_seq A reference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_AUTO (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_AUTO (auto detection)
         * @throw std::invalid\_argument if any of the following situations occur.
         * - the alignment is empty or in an incorrect format
         * - the sequences contains invalid states
         */
        void read(std::istream& aln_stream, const std::string& ref_seq = "", const InputType format = IN_AUTO, const cmaple::SeqRegion::SeqType seqtype = cmaple::SeqRegion::SEQ_AUTO);
        
        /*! \brief Read an alignment from a file in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_filename Name of an alignment file
         * @param[in] ref_seq A reference sequence (optional). If not specified, it will be read from the alignment (in MAPLE format) or automatically generated from the alignment (in FASTA or PHYLIP format)
         * @param[in] format Format of the alignment (optional): IN_MAPLE, IN_FASTA, IN_PHYLIP, or IN_AUTO (auto detection)
         * @param[in] seqtype Data type of sequences (optional): SEQ_DNA (nucleotide data), SEQ_PROTEIN (amino acid data), or SEQ_AUTO (auto detection)
         * @throw std::invalid\_argument if any of the following situations occur.
         * - the alignment is empty or in an incorrect format
         * - the sequences contains invalid states
         *
         * @throw ios::failure if the alignment file is not found
         */
        void read(const std::string& aln_filename, const std::string& ref_seq = "", const InputType format = IN_AUTO, const cmaple::SeqRegion::SeqType seqtype = cmaple::SeqRegion::SEQ_AUTO);
        
        /** \brief Write the alignment to a stream in FASTA, PHYLIP, or [MAPLE](https://www.nature.com/articles/s41588-023-01368-0) format
         * @param[in] aln_stream A stream of the output alignment file
         * @param[in] format Format of the output alignment (optional): IN_MAPLE, IN_FASTA, or IN_PHYLIP
         * @throw std::invalid\_argument if format is unknown
         * @throw std::logic\_error if the alignment is empty (i.e., nothing to write)
         */
        void write(std::ostream& aln_stream, const InputType& format = IN_MAPLE);
        
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
        void write(const std::string& aln_filename, const InputType& format = IN_MAPLE, const bool overwrite = false);
        
        // ----------------- END OF PUBLIC APIs ------------------------------------ //
        
        /**
         * @private
         * Get seq_type
         */
        inline cmaple::SeqRegion::SeqType getSeqType() const { return seq_type_; };
        
        /**
         * @private
         * Set seq_type
         */
        inline void setSeqType(cmaple::SeqRegion::SeqType seq_type) {
            seq_type_ = seq_type;
            updateNumStates();
        };
        
        /**
         @private
         Read a reference genome from an alignment file in FASTA or PHYLIP format
         @param ref_filename Name of an alignment file
         @param seq_name Name of the reference sequence
         @return a reference genome
         @throw std::invalid\_argument if seq\_name is empty
         @throw ios::failure if ref\_filename is not found
         @throw std::logic\_error if the alignment is empty or in an incorrect format
         */
         std::string readRefSeq(const std::string& ref_filename, const std::string& seq_name);
        
        /**
         @private
         Convert a state ID, indexed from 0, to a raw character
         @param state ID input a state ID
         @param seqtype The sequence type
         @return a raw state
         */
        static char convertState2Char(const cmaple::StateType& state, const cmaple::SeqRegion::SeqType& seqtype);
        
        /**
         * @private
         * Parse  alignment format from a string
         * @param n_format a format in string
         * @return an InputType
         */
        static InputType parseAlnFormat(const std::string& n_format);
        
        /**
         @private
         A vector stores all sequences
         */
        std::vector<Sequence> data; // note: this is inefficient, but only used briefly
        
        /**
         @private
         The reference sequence
         */
        std::vector<cmaple::StateType> ref_seq;
        
        /**
         @private
         The number of states
         */
        cmaple::StateType num_states;
        
        /**
         @private
         Alignment format
         */
        InputType aln_format = IN_AUTO;
        
        /**
         @private
         A set of trees that this alignment attached to
         */
        std::unordered_set<void*> attached_trees;
    
    private:
        /**
         Type of sequences
         */
        cmaple::SeqRegion::SeqType seq_type_ = cmaple::SeqRegion::SEQ_AUTO;
        
        /**
         @private
         Reset all members
         */
        void reset();
        
        /**
         @private
         update num_states according to the seq_type
         */
        void updateNumStates();
        
        /**
         @private
         detect the data type of the input sequences
         @param sequences vector of strings
         @return the data type of the input sequences
         */
        cmaple::SeqRegion::SeqType detectSequenceType(cmaple::StrVector& sequences);
        
        /**
         @private
         Compute the distance between a sequence and the ref sequence
         distance = num_differents * hamming_weight + num_ambiguities
         @param hamming_weight weight to calculate the hamming distance
         
         @throw std::logic\_error if the sequence contains an invalid type (R)
         */
        cmaple::PositionType computeSeqDistance(Sequence& sequence, cmaple::RealNumType hamming_weight);
        
        /**
         @private
         Sort sequences by their distances to the reference genome
         distance = num_differents * hamming_weight + num_ambiguities
         
         @throw std::logic\_error if the sequences contain an invalid type (R)
         */
        void sortSeqsByDistances();
        
        /**
         @private
         Convert a raw character state into ID, indexed from 0
         @param state input raw state
         @return state ID
         @throw std::invalid\_argument if state is invalid
         */
        cmaple::StateType convertChar2State(char state);
        
        /**
         @private
         Extract Mutation from sequences regarding the reference sequence
         @param sequences a vector of sequences
         @param seq_names a vector of sequence names
         @param ref_sequence the reference sequence
         @throw std::logic\_error if any of the following situations occur.
         - the lengths of sequences are different from that of the reference genome
         - sequences contains invalid states
         */
        void extractMutations(const cmaple::StrVector &sequences, const cmaple::StrVector &seq_names, const std::string& ref_sequence);
        
        /**
         @private
         Read an alignment in MAPLE format from a stream
         @param aln_stream A stream of an alignment file
         @throw std::logic\_error if any of the following situations occur.
         - the alignment is empty or in an incorrect format
         - the sequences contains invalid states
         - unexpected values/behaviors found during the operations
         */
        void readMaple(std::istream& aln_stream);
        
        /**
         @private
         Read an alignment in FASTA or PHYLIP format from a stream
         @param aln_stream A stream of an alignment file
         @param[in] ref_seq The reference sequence
         @throw std::logic\_error if any of the following situations occur.
         - the alignment is empty or in an incorrect format
         - the alignment format is unknown
         */
        void readFastaOrPhylip(std::istream& aln_stream, const std::string& ref_seq = "");
        
        /**
         @private
         Parse the reference sequence into vector of state
         @param ref_sequence reference genome in string
         @throw std::logic\_error if ref\_sequence contains invalid states
         */
        void parseRefSeq(std::string& ref_sequence);
        
        /**
         @private
         Read alignment file in FASTA format
         @param aln_stream A stream of the alignment;
         @param check_min_seqs check the minimum number of input sequences
         @return sequences, seq_names
         
         @throw std::logic\_error if the alignment is in an incorrect format
         */
        void readFasta(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, bool check_min_seqs = true);
        
        /**
         @private
         Read alignment file in PHYLIP format
         @param aln_stream A stream of the alignment;
         @param check_min_seqs: check the minimum number of input sequences
         @return sequences, seq_names
         
         @throw std::logic\_error if the alignment is in an incorrect format
         */
        void readPhylip(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, bool check_min_seqs = true);
        
        /**
         @private
         Read alignment file
         @param aln_stream A stream of the alignment;
         @param[in,out] aln_format the format of the alignment
         @param check_min_seqs: check the minimum number of input sequences
         @return sequences, seq_names
         
         @throw std::logic\_error if the alignment is in an incorrect format
         */
        void readSequences(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, InputType aln_format = IN_AUTO, bool check_min_seqs = true);
        
        /**
         @private
         Generate a reference genome from input_sequences
         @param sequences the input sequences; only_extract_diff: TRUE to only extract MAPLE file without running inference
         @return a reference genome
         
         @throw std::logic\_error if the input sequences are empty
         */
        std::string generateRef(cmaple::StrVector &sequences);
        
        /**
         Read sequence from a string line
         @throw std::logic\_error if the sequence is an incorrect format
         */
        void processSeq(std::string &sequence, std::string &line, cmaple::PositionType line_num);
        
        /**
         Add a mutation to a sequence
         @throw std::logic\_error if sequence contains invalid states
         */
        void addMutation(Sequence* sequence, char state_char, cmaple::PositionType pos, cmaple::PositionType length = -1);
        
        /**
         Write alignment in MAPLE format
         @param[in] aln_stream A stream of the output alignment file
         */
        void writeMAPLE(std::ostream& aln_stream);
        
        /**
         Write alignment in FASTA format
         @param[in] aln_stream A stream of the output alignment file
         */
        void writeFASTA(std::ostream& aln_stream);
        
        /**
         Write alignment in PHYLIP format
         @param[in] aln_stream A stream of the output alignment file
         */
        void writePHYLIP(std::ostream& aln_stream);
        
        /**
         Get reference sequence in string
         */
        std::string getRefSeqStr();
        
        /**
         Get a sequence in string
         */
        std::string getSeqString(const std::string& ref_seq_str, Sequence* sequence);
        
        /**
        Detect the format of input file in MAPLE or FASTA format
        @param aln_stream A stream of the alignment file
        @return
            IN_FASTA if in fasta format,
            IN_MAPLE if in MAPLE format
         */
        InputType detectMAPLEorFASTA(std::istream& aln_stream);
        
        /**
        Detect the format of input file
        @param aln_stream A stream of the alignment file
        @return
            IN_FASTA if in fasta format,
            IN_PHYLIP if in phylip format,
            IN_MAPLE if in MAPLE format,
            IN_UNKNOWN if file format unknown.
         */
        InputType detectInputFile(std::istream& aln_stream);
    };

    /*!
     *  \addtogroup cmaple
     *  @{
     */

    /** \brief Customized << operator to output the alignment (in MAPLE format) to a stream
     */
    std::ostream& operator<<(std::ostream& out_stream, cmaple::Alignment& aln);

    /** \brief Customized >> operator to read the alignment from a stream
     */
    std::istream& operator>>(std::istream& in_stream, cmaple::Alignment& aln);

    /*! @} End of Doxygen Groups*/

    /** Sets of valid characters/states of different data types
     */
    extern char symbols_protein[];
    extern char symbols_dna[];
    extern char symbols_rna[];
    extern char symbols_morph[];
}
#endif
