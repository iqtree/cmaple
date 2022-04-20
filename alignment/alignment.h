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

class Alignment:public vector<Sequence*> {
private:
    
    // read sequence from a string line
    void processSeq(string &sequence, string &line, PositionType line_num);
    
    // output a mutation into Diff file
    void outputMutation(ofstream &out, Sequence* sequence, char state_char, PositionType pos, PositionType length = -1);
    
    /**
        move to the next mutation in the vector of mutations
        @param sequence: a vector of mutations; mutation_index: the index of the next mutation
        @return mutation: a mutation; current_pos: the current position; end_pos: the end position
     */
    void move2NextMutation(vector<Mutation*> sequence, PositionType mutation_index, Mutation* &mutation, PositionType &current_pos, PositionType &end_pos);
    
    /**
        get the shared segment between the next mutations of two sequences
        @param current_pos: current site posisition; sequence1, sequence2: vectors of mutations; seq1_index, seq2_index: the indexes of the current mutations; seq1_end_pos, seq2_end_pos: the end positions of the current mutations
        @return seq1_mutation, seq2_mutation: the mutations contains the shared segment; length: length of the shared segment
     */
    void getNextSharedSegment(PositionType current_pos, vector<Mutation*> sequence1, vector<Mutation*> sequence2, PositionType &seq1_index, PositionType &seq2_index, Mutation* &seq1_mutation, Mutation* &seq2_mutation, PositionType &seq1_end_pos, PositionType &seq2_end_pos, PositionType &length);
    
public:
    // reference sequence
    vector<StateType> ref_seq;
    
    // type of sequences
    SeqType seq_type;
    
    // the number of states
    StateType num_states;
    
    /**
    *  Alignment constructor
    */
    Alignment();
    
    /**
    *  Alignment constructor
    */
    Alignment(vector<StateType> ref_seq, vector<Sequence*> sequences);
    
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
        @param sequences, seq_names: the input sequences,  seq_names; only_extract_diff: TRUE to only extract Diff file without running inference
        @return a reference genome
     */
    string generateRef(StrVector sequences, StrVector seq_names, bool only_extract_diff);
    
    /**
        Read a reference genome from file
        @param ref_path; only_extract_diff: TRUE to only extract Diff file without running inference
        @return a reference genome
     */
    string readRef(char* ref_path, bool only_extract_diff);
    
    /**
        extract Mutation from sequences regarding the reference sequence
        @param sequences, seq_names: the input sequences,  ref_sequence; ref_sequence, out: output stream to write the Diff file; only_extract_diff: TRUE to only extract Diff file without running inference
     */
    void extractMutations(StrVector sequences, StrVector seq_names, string ref_sequence, ofstream &out, bool only_extract_diff);
    
    /**
            Read Diff file to load reference sequence and vector of Sequence (represented by vector of Mutations)
            @param diff_path path to the Diff file; ref_path path to the reference sequence
     */
    void readDiff(char* diff_path, char* ref_path);
    
    /**
            parse the reference sequence into vector of state
            @param ref_sequence reference genome in string
     */
    void parseRefSeq(string ref_sequence);
    
    /**
        convert a state ID, indexed from 0, to a raw characer
        @param state ID input a state ID
        @return a raw state
    */
    char convertState2Char(StateType state);
    
    /**
        convert a raw characer state into ID, indexed from 0
        @param state input raw state
        @return state ID
    */
    StateType convertChar2State(char state);
    
    /**
    *  Sort sequences by their distances to the reference genome
     distance = num_differents * hamming_weight + num_ambiguities
     @param hamming_weight weight to calculate the hamming distance
    */
    void sortSeqsByDistances(double hamming_weight);
    
    /**
        compare two sequences regarding the amount of information
        @param sequence1, sequence2
        @return 0: if the two sequences are incomparable; 1: if sequence1 is more or equally informative than/to sequence2; -1; if sequence1 is less informative than sequence2
     */
    int compareSequences(vector<Mutation*> sequence1, vector<Mutation*> sequence2);
};
#endif
