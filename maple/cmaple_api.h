#pragma once
#include <ostream>
#include "cmaple.h"

/** The namespace of CMaple library */
namespace cmaple {

/** Print the copyright of CMaple */
void printCMapleCopyright(std::ostream &out);

/**
 Specify an input aligment
 @param aln path to an input alignment file
 */
void setAln(char* aln);

/**
 Specify a substitution model
 @param model name of a substitution/evolutionary model
 */
void setModel(std::string& model);

/**
 Run the inference
 */
void runInference();

/**
 Infer a phylogenetic tree from an input (alignment/DIFF) file using all default settings
 @param input path to an input (alignment/DIFF) file
 */
void inferTree(char* input);

/**
 Infer a phylogenetic tree from an input (alignment/DIFF) file using a user-specified model and all other default settings
 @param input path to an input (alignment/DIFF) file
 @param model name of a substitution/evolutionary model
 */
void inferTree(char* input, std::string& model);

/**
 Infer a phylogenetic tree from an input (alignment/DIFF) file using user settings
 @param input path to an input (alignment/DIFF) file
 @param params an instance of user settings
 */
void inferTree(char* input, Params& params);

/**
 Convert an alignment in FASTA/PHYLIP format into MAPLE format (.DIFF)
 @param input_aln path to an input alignment
 */
void extractDiff(char* input_aln);

/**
 Convert an alignment in FASTA/PHYLIP format into MAPLE format (.DIFF)
 @param input_aln path to an input alignment
 @param ref path to an reference sequence
 */
void extractDiff(char* input_aln, char* ref);

/**
 Convert a DIFF file to an alignment in FASTA format
 @param input_diff path to an input DIFF file
 */
void reconstructAln(char* input_diff);

/**
 Extend an existing tree using all default settings
 @param tree path to a tree file
 @param input path to an input (alignment/DIFF) file
 */
void extendTree(char* tree, char* input);

/**
 Extend an existing tree using a user-specified model and all other default settings
 @param tree path to a tree file
 @param input path to an input (alignment/DIFF) file
 @param model name of a substitution/evolutionary model
 */
void extendTree(char* tree, char* input, std::string& model);

/**
 Extend an existing tree using users’ settings
 @param tree path to a tree file
 @param input path to an input (alignment/DIFF) file
 @param params an instance of user settings
 */
void extendTree(char* tree, char* input, Params& params);

/**
 Infer a phylogenetic tree from an input (alignment/DIFF) file and compute aLRT-SH using all default settings
 @param input path to an input (alignment/DIFF) file
 @param num_threads number of threads (optional)
 @param random_seed an seed number for random generators (optional)
 */
void computeBranchSupports(char* input, uint32_t num_threads = 1, int64_t random_seed = -1);

/**
 Compute aLRT-SH for an existing tree using all default settings
 @param tree path to a tree file
 @param input path to an input (alignment/DIFF) file
 @param num_threads number of threads (optional)
 @param random_seed an seed number for random generators (optional)
 */
void computeBranchSupports(char* tree, char* input, uint32_t num_threads = 1, int64_t random_seed = -1);

/**
 Compute aLRT-SH for an existing tree using a user-specified model and all other default settings
 @param tree path to a tree file
 @param input path to an input (alignment/DIFF) file
 @param model name of a substitution/evolutionary model
 @param num_threads number of threads (optional)
 @param random_seed an seed number for random generators (optional)
 */
void computeBranchSupports(char* tree, char* input, std::string& model, uint32_t num_threads = 1, int64_t random_seed = -1);

/**
 Compute aLRT-SH for an existing tree using users’ settings
 @param tree path to a tree file
 @param input path to an input (alignment/DIFF) file
 @param params an instance of user settings
 @param num_threads number of threads (optional)
 @param random_seed an seed number for random generators (optional)
 */
void computeBranchSupports(char* tree, char* input, Params& params, uint32_t num_threads = 1, int64_t random_seed = -1);
    
}
