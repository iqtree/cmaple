#include "sequence.h"
#include "../utils/timeutil.h"
#include "../utils/gzstream.h"

#ifndef CMAPLE_ALIGNMENT_H
#define CMAPLE_ALIGNMENT_H

namespace cmaple
{
    extern char symbols_protein[];
    extern char symbols_dna[];
    extern char symbols_rna[];
    extern char symbols_morph[];

    /** Base class presents the input alignment */
    class AlignmentBase {
    private:
        /**
         Type of sequences
         */
        cmaple::SeqType seq_type_ = cmaple::SEQ_UNKNOWN;
        
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
        std::string getRefSeq();
        
        /**
         Get a sequence in string
         */
        std::string getSeqString(const std::string& ref_seq_str, Sequence* sequence);
        
    public:
        
        std::vector<Sequence> data; // note: this is inefficient, but only used briefly
        
        /**
         reference sequence
         */
        std::vector<cmaple::StateType> ref_seq;
        
        /**
         The number of states
         */
        cmaple::StateType num_states;
        
        /**
         Alignment format
         */
        InputType aln_format = IN_UNKNOWN;
        
        /**
         A set of trees that this alignment attached to
         */
        std::unordered_set<void*> attached_trees;
        
        /**
         *  Alignment constructor
         */
        AlignmentBase();
        
        /**
         *  Alignment deconstructor
         */
        ~AlignmentBase();
        
        /**
         *  Get seq_type
         */
        inline cmaple::SeqType getSeqType() const { return seq_type_; };
        
        /**
         *  Set seq_type
         */
        inline void setSeqType(cmaple::SeqType seq_type) {
            seq_type_ = seq_type;
            updateNumStates();
        };
        
        /**
         Read alignment file in FASTA format
         @param aln_stream A stream of the alignment;
         @param check_min_seqs check the minimum number of input sequences
         @return sequences, seq_names
         
         @throw std::logic\_error if the alignment is in an incorrect format
         */
        void readFasta(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, bool check_min_seqs = true);
        
        /**
         Read alignment file in PHYLIP format
         @param aln_stream A stream of the alignment;
         @param check_min_seqs: check the minimum number of input sequences
         @return sequences, seq_names
         
         @throw std::logic\_error if the alignment is in an incorrect format
         */
        void readPhylip(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, bool check_min_seqs = true);
        
        /**
         Read alignment file
         @param aln_stream A stream of the alignment;
         @param[in,out] aln_format the format of the alignment
         @param check_min_seqs: check the minimum number of input sequences
         @return sequences, seq_names
         
         @throw std::logic\_error if the alignment is in an incorrect format
         */
        void readSequences(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, InputType aln_format = IN_UNKNOWN, bool check_min_seqs = true);
        
        /**
         Generate a reference genome from input_sequences
         @param sequences the input sequences; only_extract_diff: TRUE to only extract MAPLE file without running inference
         @return a reference genome
         
         @throw std::logic\_error if the input sequences are empty
         */
        std::string generateRef(cmaple::StrVector &sequences);
        
        /**
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
         Read an alignment in MAPLE format from a stream
         @param aln_stream A stream of an alignment file
         @throw std::logic\_error if the alignment is empty or in an incorrect format
         */
        void readMaple(std::istream& aln_stream);
        
        /**
         Read an alignment in FASTA or PHYLIP format from a stream
         @param aln_stream A stream of an alignment file
         @param[in] ref_seq The reference sequence
         @throw std::logic\_error if any of the following situations occur.
         - the alignment is empty or in an incorrect format
         - the alignment format is unknown
         */
        void readFastaOrPhylip(std::istream& aln_stream, const std::string& ref_seq = "");
        
        /**
         Write alignment to a stream
         @param[in] ostream A stream of the output alignment file
         @param[in] format the format of the output alignment file: IN_FASTA, IN_PHYLIP, IN_MAPLE
         @throw std::invalid\_argument if format is unknown/unsupported
         */
        void write(std::ostream& aln_stream, const InputType& format);
        
        /**
         Parse the reference sequence into vector of state
         @param ref_sequence reference genome in string
         @throw std::logic\_error if ref\_sequence contains invalid states
         */
        void parseRefSeq(std::string& ref_sequence);
        
        /**
         Convert a state ID, indexed from 0, to a raw character
         @param state ID input a state ID
         @return a raw state
         */
        char convertState2Char(cmaple::StateType state);
        
        /**
         Convert a raw character state into ID, indexed from 0
         @param state input raw state
         @return state ID
         @throw std::invalid\_argument if state is invalid
         */
        cmaple::StateType convertChar2State(char state);
        
        /**
         Sort sequences by their distances to the reference genome
         distance = num_differents * hamming_weight + num_ambiguities
         
         @throw std::logic\_error if the sequences contain an invalid type (R)
         */
        void sortSeqsByDistances();
        
        /**
         Compute the distance between a sequence and the ref sequence
         distance = num_differents * hamming_weight + num_ambiguities
         @param hamming_weight weight to calculate the hamming distance
         
         @throw std::logic\_error if the sequence contains an invalid type (R)
         */
        cmaple::PositionType computeSeqDistance(Sequence& sequence, cmaple::RealNumType hamming_weight);
        
        /**
         detect the data type of the input sequences
         @param sequences vector of strings
         @return the data type of the input sequences
         */
        cmaple::SeqType detectSequenceType(cmaple::StrVector& sequences);
        
        /**
         update num_states according to the seq_type
         */
        void updateNumStates();
        
        /**
         Reset all members
         */
        void reset();
    };
}
#endif
