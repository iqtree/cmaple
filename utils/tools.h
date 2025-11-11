/****************************************************************************
 *   Copyright (C) 2022 by
 *   Nhan Ly-Trong <trongnhan.uit@gmail.com>
 *   Chris Bielow <chris.bielow@fu-berlin.de>
 *   Nicola De Maio <demaio@ebi.ac.uk>
 *   BUI Quang Minh <m.bui@anu.edu.au>
 *
 *
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the
 *   Free Software Foundation, Inc.,
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 ***************************************************************************/

#pragma once

#include <cmaple_config.h>
#include <sys/stat.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <vector>
#include <cassert>
#include "operatingsystem.h"
#ifdef _OPENMP
#include <omp.h> /* OpenMP */
#endif

/*#ifdef NDEBUG
#define ASSERT(EXPRESSION) ((void)0)
#else
#if defined(__GNUC__) || defined(__clang__)
#define ASSERT(EXPRESSION)                                             \
  ((EXPRESSION) ? (void)0                                              \
                : cmaple::_my_assert(#EXPRESSION, __PRETTY_FUNCTION__, \
                                     __FILE__, __LINE__))
#else
#define ASSERT(EXPRESSION) \
  ((EXPRESSION)            \
       ? (void)0           \
       : cmaple::_my_assert(#EXPRESSION, __func__, __FILE__, __LINE__))
#endif
#endif*/

#define USE_HASH_MAP

#if defined(__GNUC__) && !defined(GCC_VERSION)
#define GCC_VERSION \
  (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
// #else
// #define GCC_VERSION 0
#endif

// for MSVC
#ifndef __func__
#define __func__ __FUNCTION__
#endif

#if defined(USE_HASH_MAP)
//    #include <unordered_map>
//    #include <unordered_set>

#if defined(_MSC_VER)
#include <unordered_map>
#include <unordered_set>
#elif defined(__clang__)
// libc++ detected:     _LIBCPP_VERSION
// libstdc++ detected:  __GLIBCXX__
#if __has_include(<unordered_map>)  // defines _LIBCPP_VERSION
#include <unordered_map>
#include <unordered_set>
#else
#include <tr1/unordered_map>
#include <tr1/unordered_set>
using namespace std::tr1;
#endif
#elif !defined(__GNUC__)
#include <hash_map>
#include <hash_set>
using namespace stdext;
#elif GCC_VERSION < 40300
#include <ext/hash_map>
#include <ext/hash_set>
using namespace __gnu_cxx;
#define unordered_map hash_map
#define unordered_set hash_set
/*
NHANLT: resolve ambiguous reference
std::tr1 is Technical Report 1, a proposal of extensions for C++0X when waiting
for C++11 to get approved. Most of features (including
std::tr1::unordered_map/set) in Technical Report 1 had been merged into C++11.
C++11 had been fully implemented since GCC 4.8.1 -> we only need to use std::tr1
in GCC older than 4.8.1
*/
#elif GCC_VERSION < 40800
#include <tr1/unordered_map>
#include <tr1/unordered_set>
using namespace std::tr1;
#else
#include <unordered_map>
#include <unordered_set>
#endif

#else
#include <map>
#include <set>
#endif

#if defined(USE_HASH_MAP) && GCC_VERSION < 40300 && !defined(_MSC_VER) && \
    !defined(__clang__)
/*
        Define the hash function of Split
 */
#if !defined(__GNUC__)
namespace stdext
#else
namespace __gnu_cxx
#endif
{
template <>
struct hash<string> {
  size_t operator()(string str) const {
    hash<const char*> hash_str;
    return hash_str(str.c_str());
  }
};
}  // namespace stdext
#endif  // USE_HASH_MAP

namespace cmaple {
// redefine assertion
inline void _my_assert(const char* expression,
                       const char* func,
                       const char* file,
                       int line) {
  char* sfile = const_cast<char*>(strrchr(file, '/'));
  if (!sfile)
    sfile = const_cast<char*>(file);
  else
    sfile++;
  std::cerr << sfile << ":" << line << ": " << func << ": Assertion `"
            << expression << "' failed." << std::endl;
  abort();
}

/**
    Type of site states
 */
typedef uint16_t StateType;

/**
    Type of site positions
 */
typedef int32_t PositionType;

/**
    Type of the number of sequences
 */
typedef uint32_t NumSeqsType;

/**
    Type of segment lengths
 */
typedef int16_t LengthType;
/**
    Type of segment lengths (to check if 16 bit is sufficient)
 */
typedef int32_t LengthTypeLarge;

/**
    Type of real numbers
 */
typedef double RealNumType;

/**
    vector of real number number
 */
typedef std::vector<RealNumType> RealNumberVector;

/**
    vector of int
 */
typedef std::vector<int> IntList;

/**
    vector of int
 */
typedef std::vector<int> IntVector;

/**
    vector of bool
 */
typedef std::vector<bool> BoolVector;

/**
    vector of char
 */
typedef std::vector<char> CharVector;

/**
    vector of string
 */
typedef std::vector<std::string> StrVector;

/**
    Unsigned integers
 */
typedef unsigned int UINT;

/**
 *  Types of sequence regions
 */
enum RegionType : StateType {
  TYPE_R = 250,
  TYPE_O,
  TYPE_N,
  TYPE_DEL,
  TYPE_INVALID
};

/**
    verbose mode, determine how verbose should the screen be printed.
 */
enum VerboseMode { VB_QUIET, VB_MIN, VB_MED, VB_MAX, VB_DEBUG };

/**
    verbose level on the screen
 */
extern VerboseMode verbose_mode;

/**
 Index for a mininode of a phylonode
 */
enum MiniIndex : NumSeqsType {
  TOP,
  LEFT,
  RIGHT,
  UNDEFINED  // to handle case return node = null
};

/*--------------------------- NODE's INDEX -----------------------------------*/
/*! \cond PRIVATE */
/**
 Holds a space efficient index into a vector<PhyloNode> and subindex for the
 MiniNode inside the Phylonode
 */
struct Index {
  /**
      constructor
      it's required by the constructor of LeafNode
   */
    Index() : vector_index_(0), mini_index_(UNDEFINED){};

  /**
      constructor
   */
  Index(NumSeqsType vec_index, MiniIndex mini_index)
      : vector_index_(vec_index), mini_index_(mini_index){};

  /**
      return the mini-index
   */
  MiniIndex getMiniIndex() const { return mini_index_; }

  /**
      return the mini-index of the other next node (i.e., getFlipMiniIndex(LEFT)
     = RIGHT; getFlipMiniIndex(RIGHT) = LEFT)
   */
  MiniIndex getFlipMiniIndex() const { return MiniIndex(3 - mini_index_); }

  /**
      set the mini-index
   */
  void setMiniIndex(MiniIndex mini_index) { mini_index_ = mini_index; }

  /**
      return the index of a phylonode
   */
  NumSeqsType getVectorIndex() const { return vector_index_; }

  /**
      set the index of a phylonode
   */
  void setVectorIndex(NumSeqsType vec_index) { vector_index_ = vec_index; }

  /**
      Check two indexes are equally
   */
  bool operator==(const Index& rhs) const {
    return (vector_index_ == rhs.getVectorIndex() &&
            mini_index_ == rhs.getMiniIndex());
  }

 private:
  NumSeqsType vector_index_ : 30;
  MiniIndex mini_index_ : 2;
};
/*! \endcond */
/*--------------------------- NODE's INDEX -----------------------------------*/

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
const char REF_NAME[] = "REF";
const int MIN_NUM_TAXA = 3;
const RealNumType MIN_NEGATIVE =
    std::numeric_limits<RealNumType>::lowest();  // -FLT_MAX;
const RealNumType MIN_POSITIVE =
    (std::numeric_limits<RealNumType>::min)();  // FLT_MIN;
// const RealNumType INVERSE_MIN_POSITIVE = 1.0 / MIN_POSITIVE;
// const RealNumType LOG_MIN_POSITIVE = log(MIN_POSITIVE);
const RealNumType MAX_POSITIVE = (std::numeric_limits<RealNumType>::max)();
const RealNumType LOG_MAX_POSITIVE = log(MAX_POSITIVE);

template <class T>
constexpr T getMinCarryOver();
template <>
constexpr double getMinCarryOver<double>() {
  return 1e-250;
};  // of -308
template <>
constexpr float getMinCarryOver<float>() {
  return 1e-36f;
};  // of -38

const RealNumType MIN_CARRY_OVER = getMinCarryOver<RealNumType>();
const RealNumType MEAN_SUBS_PER_SITE = 0.02;
const RealNumType MAX_SUBS_PER_SITE = 0.067;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*! \cond PRIVATE */
/**
 Program parameters, everything is specified here */
class Params {
 private:
  /*! \brief Default constructor - Initializing all program parameters using the
   * default values.
   */
  Params();

  // declare ParamsBuilder as a friend class to allow it access to the default
  // contructor Params();
  friend class ParamsBuilder;

 public:
  /**
   *  Path to input sequences
   */
  std::string aln_path;

  /**
   *  Name of an alignment that contains the reference sequence
   */
  std::string ref_path;

  /**
   *  Name of the reference sequence
   */
  std::string ref_seqname;

  /**
   *  Alignment format
   */
  std::string aln_format_str;

  /**
   * Substitution model (e.g., HKY, GTR, JC, etc.)
   */
  std::string sub_model_str;

  /**
    TRUE to override the output files
   */
  bool overwrite_output;

  /**
      TRUE to fixed the branch lengths in the input tree
   */
  bool fixed_blengths;
    
    /**
        TRUE to allow CMAPLE to reroot the tree
     */
    bool allow_rerooting;

  /**
   * A relative probability threshold, which is used to ignore possible states
   * with very low probabilities. Default: 1e-8
   */
  RealNumType threshold_prob;

  /**
   * threshold_prob ^ 2
   */
  RealNumType threshold_prob2;

  /**
    The number to limit the attempts of seeking a placement for a sample
   */
  int failure_limit_sample;

  /**
    The number to limit the attempts of seeking a placement for a subtree
  */
  int failure_limit_subtree;

  /**
    The number to limit the attempts of seeking a placement for a subtree
    (short range search)
  */
  int failure_limit_subtree_short_search;

  /**
  *       TRUE to apply strict stop rules when seeking placement for a new
  sample
   */
  bool strict_stop_seeking_placement_sample;

  /**
  *       TRUE to apply strict stop rules when seeking placement for a subtree
   */
  bool strict_stop_seeking_placement_subtree;

  /**
  *       TRUE to apply strict stop rules when seeking placement for a subtree
  (short range search)
   */
  bool strict_stop_seeking_placement_subtree_short_search;
    
    /**
    *  Factor to compute the threshold of loglh to continue explore the subtree
     to seek a placement for a sample
    */
    RealNumType thresh_log_lh_sample_factor;

  /**
  *  Threshold of loglh to continue explore the subtree to seek a placement for
  a sample
  */
  RealNumType thresh_log_lh_sample;
    
    /**
    * Factor to compute the threshold of loglh to continue explore the subtree
     to seek a placement for a subtree
    */
    RealNumType thresh_log_lh_subtree_factor;

  /**
  *  Threshold of loglh to continue explore the subtree to seek a placement for
  a subtree
  */
  RealNumType thresh_log_lh_subtree;
    
    /**
    *  Factor to compute the threshold of loglh to continue explore the subtree
     to seek a placement for a subtree (short range search)
    */
    RealNumType thresh_log_lh_subtree_short_search_factor;

  /**
  *  Threshold of loglh to continue explore the subtree to seek a placement for
  a subtree (short range search)
  */
  RealNumType thresh_log_lh_subtree_short_search;

  /**
  *  Threshold of loglh to count failure
  */
  RealNumType thresh_log_lh_failure;

  /**
   * A minimum value of the branch lengths. Default: -1 (UNSPECIFIED), the
   * minimum branch length is computed from 'min_blength_factor'
   */
  RealNumType fixed_min_blength;

  /**
   *  A factor to compute the minimum branch length (if fixed_min_blength is
   * specified, this factor is ignored). Default: 0.2 <br>\<Minimum branch
   * length> = \<min_blength_factor> * \<default branch length> <br> where
   * \<default branch length> is one mutation per site.
   */
  RealNumType min_blength_factor;

  /**
   * A factor to compute the maximum branch length. Default: 40
   * <br> \<Maximum branch length> = \<max_blength_factor> * \<default branch
   * length> <br> where \<default branch length> is one mutation per site.
   */
  RealNumType max_blength_factor;

  /**
  *  The minium branch length (for mid-branch point) to try for placement
  */
  RealNumType min_blength_mid_factor;

  /**
  *  Threshold to determine whether a changed partial is different from its
  former value
  */
  RealNumType thresh_diff_update;

  /**
  *  Threshold to determine whether a changed partial is different from its
  former value by folds
  */
  RealNumType thresh_diff_fold_update;

  /**
   * The period (in term of the number of sample placements) to update the
   * substitution rate matrix. Default: 25
   */
  PositionType mutation_update_period;
    
    /**
     * The minimum number of taxa existing on the tree required to enable
     * parallel placement search. Default: 1000
     */
    NumSeqsType min_taxa_parallel_placement;
    
    /**
     * The number of samples are processed per threads
     * when searching placements in parallel. Default: 5
     */
    NumSeqsType num_samples_per_thread;
    
    /**
     * The number of steps moving upwards to extend the second (sequential)
     * search of placement starting from the placement found from
     * the first (parallel) search. Default: 2
     */
    NumSeqsType upward_search_extension;

  /**
  *  Name of the output alignment
  */
  std::string output_aln;

  /**
  *  Format of the output alignment
  */
  std::string output_aln_format_str;

  /**
  *  Path to an input tree
  */
  std::string input_treefile;

  /**
   * The number of times we traverse the tree looking for topological
   * improvements (applying SPR moves). Default: 1
   */
  int32_t num_tree_improvement;

  /**
   *  Do not try to apply SPR moves on nodes that have the placement cost (i.e.
   * the likelihood contribution by placing a node on the tree) exceeds this
   * threshold. Default: -1e-5
   */
  RealNumType thresh_placement_cost;

  /**
   *  Threshold to stop the tree search. If the total log likelihood improvement
   * obtained by an iteration of tree search is lower than this threshold,
   * CMAPLE stops doing tree search . Default: 1
   */
  RealNumType thresh_entire_tree_improvement;

  /**
  *  Don't try to re-place nodes (during short range topology search) that have
  the placement cost exceeds this threshold
  */
  RealNumType thresh_placement_cost_short_search;

  /**
  *  format of the output tree
  */
  std::string tree_format_str;

  /**
  *  TRUE to run an additional short range search for tree topology improvement
  */
  bool shallow_tree_search;

  /**
  *  path to output testing codes
  */
  char* output_testing;

  /**
  * TRUE to compute aLRT-SH
  */
  bool compute_aLRT_SH;

  /**
  * Number of replicates to compute aLRT-SH
  */
  PositionType aLRT_SH_replicates;

  /**
  * (Half) epsilon value when computing aLRT-SH
  */
  RealNumType aLRT_SH_half_epsilon;

  /**
  * number of threads
  */
  uint32_t num_threads;

  /**
  * A seed number for random generators. Default: the clock of the PC. Be
  careful! To make the results reproducible, users should specify the seed
  number.
  */
  uint64_t ran_seed;

  /**
  *  prefix output
  */
  std::string output_prefix;

  /**
  *  TRUE to allow replace the input tree by its NNI neighbor (with a higher lh)
  when computing aLRT-SH
  */
  bool allow_replace_input_tree;

  /**
  *  type of sequences
  */
  std::string seq_type_str;

  /**
  *  type of tree search
  */
  std::string tree_search_type_str;

  /**
   * TRUE to make the processes of outputting->re-inputting a tree consistent
   */
  bool make_consistent;
    
    /**
     * TRUE to ignore annotations from the input tree
     */
    bool ignore_input_annotations;
    
    /**
     * TRUE to print ids of internal nodes in the newick tree
     */
    bool print_internal_ids;
    
    /**
     * TRUE to also output tree in NEXUS format
     */
    bool output_NEXUS;
    
    /**
     * TRUE to compute the SPRTA branch supports
     */
    bool compute_SPRTA;
    
    /**
     * TRUE to compute the SPRTA for zero-length branches
     */
    bool compute_SPRTA_zero_length_branches;
    
    /**
     * TRUE to print SPRTA supports for less-infomative seqs
     */
    bool print_SPRTA_less_info_seqs;
    
    /**
     * TRUE to output the alternative SPRs with their supports in the tree
     */
    bool output_alternative_spr;
    
    /**
     * The minimum SPRTA support to be considered as
     * an alternative branches (for outputting a network)
     */
    RealNumType min_support_alt_branches;
    
    /**
     * A threshold to consider a "too short" branch as a zero-length branch
     * when computing SPRTA.
     */
    RealNumType thresh_zero_blength;
    
    /**
     * A factor determines SPRs
     * are close to the optimal one when computing SPRTA.
     * This factor is relative to the reference length.
     * Default: 1.0
     */
    RealNumType thresh_loglh_optimal_diff_fac;
    
    /**
     * A loglh threshold determines SPRs
     * are close to the optimal one when computing SPRTA.
     * It is computed from thresh_loglh_optimal_diff_fac.
     */
    RealNumType thresh_loglh_optimal_diff;
    
  /**
   * the maximum number of substitution per sites that CMAPLE is effective
  */
  RealNumType max_subs_per_site;
    
  /**
   * the mean number of substitution per sites that CMAPLE is effective
  */
  RealNumType mean_subs_per_site;
    
    /**
     * The maximum number of positive-branch descendants for
     * a subtree with a local reference
     * Default: 50
     */
    NumSeqsType max_desc_ref;
    
    /**
     * Minimum number of mutations required for a local reference
     * Default: 2
     */
    int min_mut_ref;
    
    /**
     * FALSE to turn off local references
     * Default: TRUE
     */
    bool local_refs;

  /*
      TRUE to log debugging
   */
  // bool debug = false;
};
/*! \endcond */

/**
 A PatternBuilder to build an instance of Params which stores all program
 parameters. <br> One can build an instance of Params and specify several
 parameters by:
 ParamsBuilder().withRandomSeed().withFixedMinBlength().with<...>().build();
 */
class ParamsBuilder {
 public:
  /*! \brief Default constructor - Initializing a builder using the default
   * values.
   */
  ParamsBuilder();

  /*! \brief Specify a seed number for random generators. Default: the clock of
   * the PC. Be careful! To make the results reproducible, users should specify
   * the seed number.
   * @param[in] seed A non-negative number
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if seed is negative
   */
  ParamsBuilder& withRandomSeed(const uint64_t& seed);

  /*! \brief Specify a relative probability threshold, which is used to ignore
   * possible states with very low probabilities. Default: 1e-8
   * @param[in] thresh_prob A positive relative probability threshold
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if thresh\_prob is non-positive
   */
  ParamsBuilder& withThreshProb(const double& thresh_prob);

  /*! \brief Specify a factor to compute the minimum branch length (if
   * fixed\_min\_blength is specified, this factor is ignored). Default: 0.2
   * <br>\<Minimum branch length> = \<min_blength_factor> * \<default branch
   * length> <br> where \<default branch length> is one mutation per site.
   *
   * @param[in] min_blength_factor A positive factor to compute the minimum
   * branch length
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if min\_blength\_factor is non-positive
   */
  ParamsBuilder& withMinBlengthFactor(const double& min_blength_factor);

  /*! \brief Specify a factor to compute the maximum branch length. Default: 40
   * <br> \<Maximum branch length> = \<max_blength_factor> * \<default branch
   * length> <br> where \<default branch length> is one mutation per site.
   *
   * @param[in] max_blength_factor A positive factor to compute the maximum
   * branch length
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if max\_blength\_factor is non-positive
   */
  ParamsBuilder& withMaxBlengthFactor(const double& max_blength_factor);

  /*! \brief Specify a minimum value of the branch lengths. Default: the minimum
   * branch length is computed from 'min\_blength\_factor'
   * @param[in] fixed_min_blength A positive value for the minimum branch
   * lengths
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if fixed\_min\_blength is non-positive
   */
  ParamsBuilder& withFixedMinBlength(const double& fixed_min_blength);

  /*! \brief Specify the period (in term of the number of sample placements) to
   * update the substitution rate matrix. Default: 25
   * @param[in] mutation_update_period A positive value for the period (in term
   * of the number of sample placements)
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if mutation\_update\_period is non-positive
   */
  ParamsBuilder& withMutationUpdatePeriod(
      const int32_t& mutation_update_period);

  /*! \brief Specify the number of times we traverse the tree looking for
   * topological improvements (applying SPR moves). Default: 1
   * @param[in] num_tree_traversal A positive number of tree traversals
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if num\_tree\_traversal is non-positive
   */
  ParamsBuilder& withNumTreeTraversal(const int32_t& num_tree_traversal);

  /*! \brief Specify a threshold that avoids CMAPLE trying to apply SPR moves on
   * nodes that have the placement cost (i.e. the likelihood contribution by
   * placing a node on the tree) exceeds this threshold. Default: -1e-5
   * @param[in] SPR_thresh A positive threshold to consider applying SPR moves
   * at nodes
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if SPR\_thresh is non-positive
   */
  ParamsBuilder& withSPRThresh(const double& SPR_thresh);

  /*! \brief Specify a threshold to stop the tree search. If the total log
   * likelihood improvement obtained by an iteration of tree search is lower
   * than this threshold, CMAPLE stops doing tree search . Default: 1
   * @param[in] stop_search_thresh A positive value for the threshold to stop
   * the tree search
   * @return A reference to the ParamsBuilder instance
   * @throw std::invalid\_argument if stop\_search\_thresh is non-positive
   */
  ParamsBuilder& withStopTreeSearchThresh(const double& stop_search_thresh);

  /*! \brief Build the Params object after initializing parameters
   * @return a unique pointer to an instance of Params
   */
  std::unique_ptr<cmaple::Params> build();

 private:
  /**
   * Pointer to a Params object
   */
  std::unique_ptr<Params> params_ptr{};
};

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
    Print error message then exit program
 */
void outError(const char* error, bool quit = true);

/**
    Print error message then exit program
 */
void outError(const std::string& error, bool quit = true);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
    Print error messages then exit program
 */
void outError(const char* error, const char* msg, bool quit = true);

/**
    Print error messages then exit program
 */
void outError(const char* error, const std::string& msg, bool quit = true);

/**
    Output a warning message to screen
    @param error warning message
 */
void outWarning(const char* warn);
void outWarning(const std::string& warn);

/** safe version of std::getline to deal with files from different platforms */
std::istream& safeGetline(std::istream& is, std::string& t);

/**
        @return TRUE of ch is a control character (ascii <= 32)
 */
inline bool controlchar(char& ch) {
  return ch <= 32;
}

inline bool is_newick_token(char& ch) {
  return ch == ':' || ch == ';' || ch == ',' || ch == ')' || ch == '(' ||
         ch == '[' || ch == ']';
}

/**
    change unusual character in names into underscore (_)
    @param[in/out] name string name
    @return true if renamed, false otherwise
 */
bool renameString(std::string& name);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*
        Error messages
 */
const char ERR_NO_TAXON[] = "Find no taxon with name ";
const char ERR_NO_AREA[] = "Find no area with name ";
const char ERR_NO_MEMORY[] = "Not enough memory!";
const char ERR_NEG_BRANCH[] = "Negative branch length not allowed.";
const char ERR_READ_INPUT[] =
    "File not found or incorrect input, pls check it again.";
const char ERR_READ_ANY[] =
    "Unidentified error while reading file, pls check it carefully again.";

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
 * Convert int to string
 * @param int
 * @return string
 */
std::string convertPosTypeToString(PositionType number);
std::string convertIntToString(int number);
std::string convertInt64ToString(int64_t number);

std::string convertDoubleToString(RealNumType number);

std::string convertDoubleToString(RealNumType number, uint8_t precision);

/**
 Replace the first subtring, found in a string, by a new one
 */
void replaceSubStr(std::string& input_str,
                   const std::string& old_sub_str,
                   const std::string& new_sub_str);

/**
 Case-insensitive comparison between two strings
 @return true if two strings are equal.
 */
bool iEquals(const std::string& a, const std::string& b);

/**
 *
 * @param SRC
 * @param DEST
 * @return bool
 */
bool copyFile(const char SRC[], const char DEST[]);

/**
 * Check if the file exists
 * @param strFilename
 * @return
 */
bool fileExists(const std::string& strFilename);

/**
    Check that path is a directory
 */
int isDirectory(const char* path);

/**
    Convert string to int, with error checking
    @param str original string
    @return the number
    @throw std::invalid\_argument if the input str is invalid
 */
int convert_int(const char* str);

/**
    Convert string to int64, with error checking
    @param str original string
    @return the number
    @throw std::invalid\_argument if the input str is invalid
 */
int64_t convert_int64(const char* str);

/**
    Convert string to int, with error checking
    @param str original string
    @param end_pos end position
    @return the number
    @throw std::invalid\_argument if the input str is invalid
 */
int convert_int(const char* str, int& end_pos);

/**
    Convert comma-separated string to integer vector, with error checking
    @param str original string with integers separated by comma
    @param vec (OUT) integer vector
    @throw std::invalid\_argument if the input str is invalid
 */
void convert_int_vec(const char* str, IntVector& vec);

/**
    Convert string to int64_t, with error checking
    @param str original string
    @return the number
    @throw std::invalid\_argument if the input str is invalid
 */
int64_t convert_int64(const char* str);

/**
    Convert string to int64_t, with error checking
    @param str original string
    @param end_pos end position
    @return the number
    @throw std::invalid\_argument if the input str is invalid
 */
int64_t convert_int64(const char* str, int& end_pos);

/**
    Convert string to a real number, with error checking
    @param str original string
    @return the RealNumType
    @throw std::invalid\_argument if the input str is invalid
 */
RealNumType convert_real_number(const char* str);

/**
    Convert string to real number, with error checking
    @param str original string
    @param end_pos end position
    @return the RealNumType
    @throw std::invalid\_argument if the input str is invalid
 */
RealNumType convert_real_number(const char* str, int& end_pos);

/**
    Parse an array of real numbers from a tring
    @param input_str: a string of  real_numbers; arr: the output array of
   real_numbers
    @throw std::invalid\_argument if the input str is invalid
 */
void convert_real_numbers(RealNumType*& arr, std::string input_str);

/**
    Convert comma-separated string to integer vector, with error checking
    @param str original string with integers separated by comma
    @param vec (OUT) integer vector
    @param separator char separating elements
    @throw std::invalid\_argument if the input str is invalid
 */
void convert_real_number_vec(const char* str,
                             RealNumberVector& vec,
                             char separator = ',');

/**
    Normalize state frequencies so that sum of them is equal to 1
    @param freqs original state frequencies
    @param num_states the number of states
    @param total_freq sum of all original state frequencies
    @throw std::logic\_error if sum of freqs is zero
 */
void normalize_frequencies_from_index(RealNumType* freqs,
                                      int num_states,
                                      int starting_index);

/**
    Normalize entries so that sum of them is equal to 1
    @param entries original entries
    @param num_entries the number of entries
    @param sum_entries Precomputed sum of all original state frequencies
 */
inline void normalize_arr(RealNumType* const entries,
                          const int num_entries,
                          RealNumType sum_entries) {
  assert(num_entries > 0);

#ifndef NDEBUG
  // if (fabs(sum_entries) < 1e-5)
  //   outError("Sum of entries must be greater than zero!");
#endif

  sum_entries = 1.0 / sum_entries;
  for (int i = 0; i < num_entries; ++i)
    entries[i] *= sum_entries;
    // entries[i] /= sum_entries;
}

/**
    Normalize entries so that sum of them is equal to 1
    @param entries original entries
    @param num_entries the number of entries
 */
inline void normalize_arr(RealNumType* const entries, const int num_entries) {
  RealNumType sum_entries = 0;
  for (int i = 0; i < num_entries; ++i)
    sum_entries += entries[i];
  normalize_arr(entries, num_entries, sum_entries);
}
/**
 * Convert seconds to hour, minute, second
 * @param sec
 * @return string represent hour, minute, second
 */
std::string convert_time(const RealNumType sec);

/**
    Convert a string to to range lower:upper:step_size with error checking
    @param str original string
    @param lower (OUT) lower bound of the range
    @param upper (OUT) upper bound of the range
    @param step_size (OUT) step size of the range
    @throw std::invalid\_argument if the input str is invalid
 */
void convert_range(const char* str, int& lower, int& upper, int& step_size);

/**
    Convert a string to to range lower:upper:step_size with error checking
    @param str original string
    @param lower (OUT) lower bound of the range
    @param upper (OUT) upper bound of the range
    @param step_size (OUT) step size of the range
    @throw std::invalid\_argument if the input str is invalid
 */
void convert_range(const char* str,
                   RealNumType& lower,
                   RealNumType& upper,
                   RealNumType& step_size);

/**
    Reinitialize an array of real number (RealNumType*)
    @param arr the input RealNumType*
    @param size the size of the input array
    @param delete_first TRUE to delete the current array before reinitializing
    @param set_zero TRUE to initialize all new items at 0
 */
void reinitDoubleArr(RealNumType*& arr,
                     StateType size,
                     bool delete_first = true,
                     bool set_zero = true);

/**
    Convert char* into vector of strings separated by separator
 */
void convert_string_vec(const char* str,
                        StrVector& str_vec,
                        char separator = ',');

/**
    Convert string to PositionType, with error checking
    @param str original string
    @return the number
    @throw std::invalid\_argument if the input str is invalid
 */
PositionType convert_positiontype(const char* str);

/**
    Check whether a string is a number
    @param s storing a string
 */
bool is_number(const std::string& s);

/**
    Parse program argument into params
    @param argc number of arguments
    @param argv list of arguments
    @param params (OUT) program parameters
 */
void parseArg(int argc, char* argv[], Params& params);

/**
    Show quick start guide
 */
void quickStartGuide();

/**
    Show help
 */
void usage_cmaple();

/**
    Remove white space at the beginning and end of the string
    @param str (IN/OUT) string to be trimmed
*/
void trimString(std::string& str);

/**
    get number of processor cores
*/
int countPhysicalCPUCores();

/**
    Sort an array by quicksort
    @param arr array for comparison
    @param arr2 array for sorting
    @param left the left index
    @param right the right index
*/
template <class T1, class T2>
void quicksort(T1* arr, int left, int right, T2* arr2 = NULL) {
  if (left > right)
    return;
  assert(left <= right);
  int i = left, j = right;
  T1 pivot = arr[(left + right) / 2];

  /* partition */
  while (i <= j) {
    while (arr[i] < pivot)
      ++i;
    while (arr[j] > pivot)
      --j;
    if (i <= j) {
      T1 tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;
      if (arr2) {
        T2 tmp2 = arr2[i];
        arr2[i] = arr2[j];
        arr2[j] = tmp2;
      }
      ++i;
      --j;
    }
  };

  /* recursion */
  if (left < j)
    quicksort(arr, left, j, arr2);
  if (i < right)
    quicksort(arr, i, right, arr2);
}

/**
 * Print usage for CMAPLE
 * @param program arguments list
 * @param full_command TRUE to print all available commands, FALSE to print
 * normal usage dialog
 */
void usage_cmaple(char* argv[], bool full_command);

/**
 * Print copyright
 */
void printCopyright(std::ostream& out);

/**
 * Set number of threads
 * @param num_threads the number of OpenMP threads
 * @throw std::invalid\_argument if num\_threads > the number of CPU cores
 */
void setNumThreads(const int num_threads);

/**
 * make unique pointer, an implementation to make it compatible to c++11
 */
template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/**
 * Reset a stream (to the beginning)
 * @param instream the input stream
 */
void resetStream(std::istream& instream);
}  // namespace cmaple
