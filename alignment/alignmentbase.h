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
         */
        void processSeq(std::string &sequence, std::string &line, cmaple::PositionType line_num);
        
        /**
         Add a mutation to a sequence
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
         */
        void readFasta(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, bool check_min_seqs = true);
        
        /**
         Read alignment file in PHYLIP format
         @param aln_stream A stream of the alignment;
         @param check_min_seqs: check the minimum number of input sequences
         @return sequences, seq_names
         */
        void readPhylip(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, bool check_min_seqs = true);
        
        /**
         Read alignment file
         @param aln_stream A stream of the alignment;
         @param check_min_seqs: check the minimum number of input sequences
         @return sequences, seq_names
         */
        void readSequences(std::istream& aln_stream, cmaple::StrVector &sequences, cmaple::StrVector &seq_names, bool check_min_seqs = true);
        
        /**
         Generate a reference genome from input_sequences
         @param sequences the input sequences; only_extract_diff: TRUE to only extract MAPLE file without running inference
         @return a reference genome
         */
        std::string generateRef(cmaple::StrVector &sequences);
        
        // DISABLE due to new implementation
        /**
         Read a reference genome from file
         @param ref_path an alignment contains the reference sequence at the beginning of the file
         @return a reference genome
         */
        // std::string readRef(const std::string& ref_path);
        
        /**
         Extract Mutation from sequences regarding the reference sequence
         @param sequences a vector of sequences
         @param seq_names a vector of sequence names
         @param ref_sequence the reference sequence
         */
        void extractMutations(const cmaple::StrVector &sequences, const cmaple::StrVector &seq_names, const std::string& ref_sequence);
        
        /**
         Read an alignment in MAPLE format from a stream
         @param aln_stream A stream of an alignment file
         */
        void readMaple(std::istream& aln_stream);
        
        /**
         Read an alignment in FASTA or PHYLIP format from a stream
         @param aln_stream A stream of an alignment file
         @param[in] ref_seq The reference sequence
         */
        void readFastaOrPhylip(std::istream& aln_stream, const std::string& ref_seq = "");
        
        /**
         Write alignment to a stream
         @param[in] ostream A stream of the output alignment file
         @param[in] format the format of the output alignment file: IN_FASTA, IN_PHYLIP, IN_MAPLE
         */
        void write(std::ostream& aln_stream, const InputType& format);
        
        /**
         Parse the reference sequence into vector of state
         @param ref_sequence reference genome in string
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
         */
        cmaple::StateType convertChar2State(char state);
        
        /**
         Sort sequences by their distances to the reference genome
         distance = num_differents * hamming_weight + num_ambiguities
         @param hamming_weight weight to calculate the hamming distance
         */
        void sortSeqsByDistances(cmaple::RealNumType hamming_weight);
        
        /**
         Compute the distance between a sequence and the ref sequence
         distance = num_differents * hamming_weight + num_ambiguities
         @param hamming_weight weight to calculate the hamming distance
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
         Get alignment format from a string
         */
        InputType getAlignmentFormat(const std::string& n_format);
        
        /**
         Get alignment format from a string
         */
        SeqType getSeqType(const std::string& n_seqtype_str);
    };
}
#endif
