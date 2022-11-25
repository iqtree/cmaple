//
//  alignment.h
//  alignment
//
//  Created by Nhan Ly-Trong on 24/01/2022.
//

#include "sequence.h"
#include "utils/timeutil.h"
#include "utils/gzstream.h"

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

/** Class presents the input alignment */
class Alignment {
private:
    /**
        Read sequence from a string line
     */
    void processSeq(std::string &sequence, std::string &line, PositionType line_num);
    
    /**
        Output a mutation into Diff file
     */
    void outputMutation(std::ofstream &out, Sequence* sequence, char state_char, PositionType pos, PositionType length = -1);
    
public:

    std::vector<Sequence*> data; // note: this is inefficient, but only used briefly

    /**
        reference sequence
     */
    std::vector<StateType> ref_seq;
    
    /**
        Type of sequences
     */
    SeqType seq_type = SEQ_DNA;
    
    /**
        The number of states
     */
    StateType num_states = 4;
    
    /**
    *  Alignment constructor
    */
    Alignment();
    
    /**
    *  Alignment deconstructor
    */
    ~Alignment();
    
    /**
        Read alignment file in FASTA format
        @param aln_path path to the alignment; check_min_seqs: check the minimum number of input sequences
        @return sequences, seq_names
     */
    void readFasta(char *aln_path, StrVector &sequences, StrVector &seq_names, bool check_min_seqs = true);
    
    /**
        Read alignment file in PHYLIP format
        @param aln_path path to the alignment; check_min_seqs: check the minimum number of input sequences
        @return sequences, seq_names
     */
    void readPhylip(char *aln_path, StrVector &sequences, StrVector &seq_names, bool check_min_seqs = true);
    
    /**
        Read alignment file
        @param aln_path path to the alignment; check_min_seqs: check the minimum number of input sequences
        @return sequences, seq_names
     */
    void readSequences(char* aln_path, StrVector &sequences, StrVector &seq_names, bool check_min_seqs = true);
    
    /**
        Generate a reference genome from input_sequences
        @param sequences the input sequences; only_extract_diff: TRUE to only extract Diff file without running inference
        @return a reference genome
     */
    std::string generateRef(StrVector &sequences, bool only_extract_diff);
    
    /**
        Read a reference genome from file
        @param ref_path; only_extract_diff: TRUE to only extract Diff file without running inference
        @return a reference genome
     */
    std::string readRef(char* ref_path, bool only_extract_diff);
    
    /**
        Extract Mutation from sequences regarding the reference sequence
        @param sequences, seq_names: the input sequences,  ref_sequence; ref_sequence, out: output stream to write the Diff file; only_extract_diff: TRUE to only extract Diff file without running inference
     */
    void extractMutations(StrVector &sequences, StrVector &seq_names, std::string ref_sequence, std::ofstream &out, bool only_extract_diff);
    
    /**
        Read Diff file to load reference sequence and vector of Sequence (represented by vector of Mutations)
        @param diff_path path to the Diff file; ref_path path to the reference sequence
     */
    void readDiff(char* diff_path, char* ref_path);
    
    /**
        Reconstruct an alignment file from a Diff file
        @param diff_path path to the Diff file; output_file path to the output aln
     */
    void reconstructAln(char* diff_path, char* output_file);
    
    /**
        Extract Diff file from and alignment file
     */
    void extractDiffFile(Params& params);
    
    /**
        Parse the reference sequence into vector of state
        @param ref_sequence reference genome in string
     */
    void parseRefSeq(std::string ref_sequence);
    
    /**
        Convert a state ID, indexed from 0, to a raw character
        @param state ID input a state ID
        @return a raw state
    */
    char convertState2Char(StateType state);
    
    /**
        Convert a raw character state into ID, indexed from 0
        @param state input raw state
        @return state ID
    */
    StateType convertChar2State(char state);
    
    /**
        Sort sequences by their distances to the reference genome
        distance = num_differents * hamming_weight + num_ambiguities
        @param hamming_weight weight to calculate the hamming distance
    */
    void sortSeqsByDistances(RealNumType hamming_weight);
};
#endif
