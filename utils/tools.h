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

#include <cmaple_config.h>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <stack>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <stdint.h>
#include <string.h>
#include <sstream>
#include <random>
#include <sys/stat.h>
#include <cfloat>
#include <thread>
#include <iomanip>
#include <regex>
#define SPRNG
#include "../libraries/sprng/sprng.h"
#include "operatingsystem.h"
#ifdef _OPENMP
#include <omp.h> /* OpenMP */
#endif

#ifndef CMAPLE_TOOLS_H
#define CMAPLE_TOOLS_H

#ifdef NDEBUG
#define ASSERT(EXPRESSION) ((void)0)
#else
    #if defined(__GNUC__) || defined(__clang__)
        #define ASSERT(EXPRESSION) ((EXPRESSION) ? (void)0 : cmaple::_my_assert(#EXPRESSION, __PRETTY_FUNCTION__, __FILE__, __LINE__))
    #else
        #define ASSERT(EXPRESSION) ((EXPRESSION) ? (void)0 : cmaple::_my_assert(#EXPRESSION, __func__, __FILE__, __LINE__))
    #endif
#endif


#define USE_HASH_MAP

#if defined(__GNUC__) && !defined(GCC_VERSION)
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
//#else
//#define GCC_VERSION 0
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
        #if __has_include(<unordered_map>) // defines _LIBCPP_VERSION
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
    std::tr1 is Technical Report 1, a proposal of extensions for C++0X when waiting for C++11 to get approved.
    Most of features (including std::tr1::unordered_map/set) in Technical Report 1 had been merged into C++11.
    C++11 had been fully implemented since GCC 4.8.1 -> we only need to use std::tr1 in GCC older than 4.8.1
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


#if    defined(USE_HASH_MAP) && GCC_VERSION < 40300 && !defined(_MSC_VER) && !defined(__clang__)
/*
        Define the hash function of Split
 */
#if !defined(__GNUC__)
namespace stdext {
#else
namespace __gnu_cxx {
#endif

    template<>
    struct hash<string> {

        size_t operator()(string str) const {
            hash<const char*> hash_str;
            return hash_str(str.c_str());
        }
    };
} // namespace
#endif // USE_HASH_MAP

namespace cmaple
{
    // redefine assertion
    inline void _my_assert(const char* expression, const char *func, const char* file, int line)
    {
        char *sfile = (char*)strrchr(file, '/');
        if (!sfile) sfile = (char*)file; else sfile++;
        std::cerr << sfile << ":" << line << ": " << func << ": Assertion `" << expression << "' failed." << std::endl;
        abort();
    }

    struct Distribution {
      std::string random_numbers_str;
      int pool_size;
    } ;

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
        input type, tree or splits graph
     */
    enum InputType {
        IN_NEWICK, IN_NEXUS, IN_FASTA, IN_PHYLIP, IN_COUNTS, IN_CLUSTAL, IN_MSF, IN_MAPLE, IN_OTHER
    };

    /**
        verbose mode, determine how verbose should the screen be printed.
     */
    enum VerboseMode {
        VB_QUIET, VB_MIN, VB_MED, VB_MAX, VB_DEBUG
    };

    /**
        types of sequences
     */
    enum SeqType {
        SEQ_DNA, SEQ_PROTEIN, SEQ_BINARY, SEQ_MORPH, SEQ_MULTISTATE, SEQ_CODON, SEQ_POMO, SEQ_UNKNOWN
    };

    /**
        verbose level on the screen
     */
    extern VerboseMode verbose_mode;

    /** Index for a mininode of a phylonode */
    enum MiniIndex : NumSeqsType
    {
        TOP,
        LEFT,
        RIGHT,
        UNDEFINED // to handle case return node = null
    };

    /*--------------------------- NODE's INDEX -----------------------------------*/
    /** Holds a space efficient index into a vector<PhyloNode> and subindex for the MiniNode inside the Phylonode */
    struct Index
    {
        /**
            constructor
            it's required by the constructor of LeafNode
         */
        Index(): vector_index_(0), mini_index_(UNDEFINED) {};
        
        /**
            constructor
         */
        Index(NumSeqsType vec_index, MiniIndex mini_index): vector_index_(vec_index), mini_index_(mini_index) {};

        /**
            return the mini-index
         */
        const MiniIndex getMiniIndex() const
        {
            return mini_index_;
        }
        
        /**
            return the mini-index of the other next node (i.e., getFlipMiniIndex(LEFT) = RIGHT; getFlipMiniIndex(RIGHT) = LEFT)
         */
        const MiniIndex getFlipMiniIndex() const
        {
            return MiniIndex(3 - mini_index_);
        }
        
        /**
            set the mini-index
         */
        void setMiniIndex(MiniIndex mini_index)
        {
            mini_index_ = mini_index;
        }

        /**
            return the index of a phylonode
         */
        const NumSeqsType getVectorIndex() const
        {
            return vector_index_;
        }
        
        /**
            set the index of a phylonode
         */
        void setVectorIndex(NumSeqsType vec_index)
        {
            vector_index_ = vec_index;
        }
        
        /**
            Check two indexes are equally
         */
        bool operator==(Index &rhs) const
        {
            return (vector_index_ == rhs.getVectorIndex() && mini_index_ == rhs.getMiniIndex());
        }
        
    private:
        NumSeqsType vector_index_ : 30 {0};
        MiniIndex mini_index_ : 2 {UNDEFINED};
    };
    /*--------------------------- NODE's INDEX -----------------------------------*/


    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    const char REF_NAME[] = "REF";
    const int MIN_NUM_TAXA = 3;
    const RealNumType MIN_NEGATIVE = std::numeric_limits<RealNumType>::lowest();// -FLT_MAX;
    const RealNumType MIN_POSITIVE = (std::numeric_limits<RealNumType>::min)();// FLT_MIN;
    //const RealNumType INVERSE_MIN_POSITIVE = 1.0 / MIN_POSITIVE;
    //const RealNumType LOG_MIN_POSITIVE = log(MIN_POSITIVE);
    const RealNumType MAX_POSITIVE = (std::numeric_limits<RealNumType>::max)();
    const RealNumType LOG_MAX_POSITIVE = log(MAX_POSITIVE);

    template <class T> constexpr T getMinCarryOver();
    template<> constexpr double getMinCarryOver<double>() { return 1e-250; }; // of -308
    template<> constexpr float getMinCarryOver<float>() { return 1e-36f; };    // of -38

    const RealNumType MIN_CARRY_OVER = getMinCarryOver<RealNumType>();



    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/

    /**
            program parameters, everything is specified here
            Use singleton pattern to avoid using global variable or
            having to pass the params variable around
     */
    class Params {
    public:
        static Params& getInstance();
    private:
        Params () {}; // Disable constructor
    public:
        
        /**
        *  Path to input sequences
        */
        char* aln_path;
        
        /**
        *  Path to a Diff file
        */
        char* diff_path;
        
        /**
        *  Path to the reference sequence
        */
        char* ref_path;
        
        /**
        *  TRUE to only extract Diff file (from alignment) without running inference
        */
        bool only_extract_diff;
        
        /**
        *  Weight to calculate the hamming distance
        */
        RealNumType hamming_weight;
        
        /**
            Name of the substitution model (e.g., HKY, GTR, TN+I+G, JC+G, etc.)
         */
        std::string model_name;
        
        /**
            TRUE to override the output files
         */
        bool overwrite_output;
        
        /**
            TRUE to redo the inference and overwrite output files
         */
        bool redo_inference;
        
        /**
            TRUE to re-estimate the branch lengths
         */
        bool redo_blength;
        
        /**
        *       Threshold to ignore possible states with very low probabilities
         */
        RealNumType threshold_prob;
        
        /**
        *       threshold_prob ^ 2
         */
        RealNumType threshold_prob2;
        
        /**
        *      The number to limit the attempts of seeking a placement for a sample
        */
        int failure_limit_sample;
        
        /**
        *      The number to limit the attempts of seeking a placement for a subtree
        */
        int failure_limit_subtree;
        
        /**
        *      The number to limit the attempts of seeking a placement for a subtree (short range search
        */
        int failure_limit_subtree_short_search;
        
        /**
        *       The period to update the mutation matrix
         */
        PositionType mutation_update_period;
        
        /**
        *       TRUE to apply strict stop rules when seeking placement for a new sample
         */
        bool strict_stop_seeking_placement_sample;
        
        /**
        *       TRUE to apply strict stop rules when seeking placement for a subtree
         */
        bool strict_stop_seeking_placement_subtree;
        
        /**
        *       TRUE to apply strict stop rules when seeking placement for a subtree (short range search)
         */
        bool strict_stop_seeking_placement_subtree_short_search;
        
        /**
        *  Threshold of loglh to continue explore the subtree to seek a placement for a sample
        */
        RealNumType thresh_log_lh_sample;
        
        /**
        *  Threshold of loglh to continue explore the subtree to seek a placement for a subtree
        */
        RealNumType thresh_log_lh_subtree;
        
        /**
        *  Threshold of loglh to continue explore the subtree to seek a placement for a subtree (short range search)
        */
        RealNumType thresh_log_lh_subtree_short_search;
        
        /**
        *  Threshold of loglh to count failure
        */
        RealNumType thresh_log_lh_failure;
        
        /**
        *  The minium branch length to try for placement
        */
        RealNumType min_blength_factor;
        
        /**
        *  The maximum branch length to try for placement
        */
        RealNumType max_blength_factor;
        
        /**
        *  The minium branch length (for mid-branch point) to try for placement
        */
        RealNumType min_blength_mid_factor;
        
        /**
        *  Threshold to determine whether a changed partial is different from its former value
        */
        RealNumType thresh_diff_update;
        
        /**
        *  Threshold to determine whether a changed partial is different from its former value by folds
        */
        RealNumType thresh_diff_fold_update;
        
        /**
        *  Path to output aln
        */
        char* output_aln;
        
        /**
        *  Path to an input tree
        */
        char* input_treefile;
        
        /**
        *       The number of attempts to improve the tree with SPR moves
        */
        int num_tree_improvement;
        
        /**
        *  Threshold for improvement on the entire tree
        */
        RealNumType thresh_entire_tree_improvement;
        
        /**
        *  Don't try to re-place nodes that have the placement cost exceeds this threshold
        */
        RealNumType thresh_placement_cost;
        
        /**
        *  Don't try to re-place nodes (during short range topology search) that have the placement cost exceeds this threshold
        */
        RealNumType thresh_placement_cost_short_search;
        
        /**
        *  export output tree in binary tree format
        */
        bool export_binary_tree;
        
        /**
        *  TRUE to optimize the branch lengths after optimizing the tree topology
        */
        bool optimize_branch_length;
        
        /**
        *  TRUE to run an additional short range search for tree topology improvement
        */
        bool short_range_topo_search;
        
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
        * epsilon value when computing aLRT-SH
        */
        RealNumType aLRT_SH_epsilon;
        
        /**
        * fixed min_blength
        */
        RealNumType fixed_min_blength;
        
        /**
        * number of threads
        */
        uint32_t num_threads;
        
        /**
        * Random seed
        */
        uint64_t ran_seed;
        
        /**
        *  prefix output
        */
        char* output_prefix;
        
        /**
        *  TRUE to allow replace the input tree by its NNI neighbor (with a higher lh) when computing aLRT-SH
        */
        bool allow_replace_input_tree;
        
        /**
        *  type of sequences
        */
        SeqType seq_type;
        
        /*
            TRUE to log debugging
         */
        // bool debug = false;
    };

    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/

    /**
        Print error message then exit program
     */
    void outError(const char *error, bool quit = true);

    /**
        Print error message then exit program
     */
    void outError(const std::string &error, bool quit = true);


    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/

    /**
        Print error messages then exit program
     */
    void outError(const char *error, const char *msg, bool quit = true);

    /**
        Print error messages then exit program
     */
    void outError(const char *error, const std::string &msg, bool quit = true);

    /**
        Output a warning message to screen
        @param error warning message
     */
    void outWarning(const char *warn);
    void outWarning(const std::string &warn);


    /** safe version of std::getline to deal with files from different platforms */
    std::istream& safeGetline(std::istream& is, std::string& t);

    /**
            @return TRUE of ch is a control character (ascii <= 32)
     */
    inline bool controlchar(char& ch) {
        return ch <= 32;
    }

    inline bool is_newick_token(char& ch) {
        return ch == ':' || ch == ';' || ch == ',' || ch == ')' || ch == '(' || ch == '[' || ch == ']';
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
    const char ERR_READ_INPUT[] = "File not found or incorrect input, pls check it again.";
    const char ERR_READ_ANY[] = "Unidentified error while reading file, pls check it carefully again.";

    /*--------------------------------------------------------------*/
    /* random number generator */
    /*--------------------------------------------------------------*/

    extern int *randstream;

    /**
     * initialize the random number generator
     * @param seed seed for generator
     * @param write_info true to write information, false otherwise (default)
     */
    int init_random(int seed, bool write_info = false, int** rstream = NULL);

    /**
     * finalize random number generator (e.g. free memory
     */
    int finish_random(int *rstream = NULL);

    /**
     * returns a random integer in the range [0; n - 1]
     * @param n upper-bound of random number
     */
    int random_int(int n, int *rstream = NULL);

    /**
     *  return a random integer in the range [a,b]
     */
    //int randint(int a, int b);

    /**
     * returns a random integer in the range [0; RAND_MAX - 1]
     * = random_int(RAND_MAX)
     */
    //int random_int(int *rstream = NULL);

    /**
     * returns a random floating-point nuber in the range [0; 1)
     */
    double random_double(int *rstream = NULL);

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
     Case-insensitive comparison between two strings
     @return true if two strings are equal.
     */
    bool iEquals(const std::string &a, const std::string &b);
        
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
    bool fileExists(const std::string &strFilename);

    /**
        Check that path is a directory
     */
    int isDirectory(const char *path);

    /**
        Convert string to int, with error checking
        @param str original string
        @return the number
     */
    int convert_int(const char *str);

    /**
        Convert string to int64, with error checking
        @param str original string
        @return the number
     */
    int64_t convert_int64(const char *str);

    /**
        Convert string to int, with error checking
        @param str original string
        @param end_pos end position
        @return the number
     */
    int convert_int(const char *str, int &end_pos);

    /**
        Convert comma-separated string to integer vector, with error checking
        @param str original string with integers separated by comma
        @param vec (OUT) integer vector
     */
    void convert_int_vec(const char *str, IntVector &vec);

    /**
        Convert string to int64_t, with error checking
        @param str original string
        @return the number
     */
    int64_t convert_int64(const char *str);

    /**
        Convert string to int64_t, with error checking
        @param str original string
        @param end_pos end position
        @return the number
     */
    int64_t convert_int64(const char *str, int &end_pos);

    /**
        Convert string to a real number, with error checking
        @param str original string
        @return the RealNumType
     */
    RealNumType convert_real_number(const char *str);

    /**
        Convert string to real number, with error checking
        @param str original string
        @param end_pos end position
        @return the RealNumType
     */
    RealNumType convert_real_number(const char *str, int &end_pos);

    /**
        Parse an array of real numbers from a tring
        @param input_str: a string of  real_numbers; arr: the output array of real_numbers
     */
    void convert_real_numbers(RealNumType* &arr, std::string input_str);

    /**
        Convert string to real number, or generate it from a distribution
        @param str original string
        @param end_pos end position
        @param separator char separating elements
        @return the real number
     */
    RealNumType convert_real_number_with_distribution(const char *str, int &end_pos, char separator = ',');

    /**
        Convert comma-separated string to integer vector, with error checking
        @param str original string with integers separated by comma
        @param vec (OUT) integer vector
        @param separator char separating elements
     */
    void convert_real_number_vec(const char *str, RealNumberVector &vec, char separator = ',');

    /**
        Convert comma-separated string to real number vector or generate real number vector from distributions
        @param str original string with real numbers separated by comma
        @param vec (OUT) real number vector
        @param separator char separating elements
     */
    void convert_real_number_vec_with_distributions(const char *str, RealNumberVector &vec, char separator = ',');

    /**
        Convert separated string to an array of real number (RealNumType*) or generate them from distributions
        @param tmp_str original string with real numbers separated by separator
        @param array an array of RealNumType number (RealNumType*)
        @param num_items the number of items in the array
        @param separator char separating elements
     */
    void convert_real_number_array_with_distributions(std::string tmp_str, RealNumType* array, int num_items, char separator);

    /**
        Normalize state frequencies so that sum of them is equal to 1
        @param freqs original state frequencies
        @param num_states the number of states
        @param total_freq sum of all original state frequencies
     */
    void normalize_frequencies_from_index(RealNumType* freqs, int num_states, int starting_index);

    /**
        Normalize entries so that sum of them is equal to 1
        @param entries original entries
        @param num_entries the number of entries
        @param sum_entries Precomputed sum of all original state frequencies
     */
    inline void normalize_arr(RealNumType* const entries, const int num_entries, RealNumType sum_entries)
    {
      ASSERT(num_entries > 0);

    #ifndef NDEBUG
      //if (fabs(sum_entries) < 1e-5)
      //  outError("Sum of entries must be greater than zero!");
    #endif

      sum_entries = 1.0 / sum_entries;
      for (int i = 0; i < num_entries; ++i)
        entries[i] *= sum_entries;
    }

    /**
        Normalize entries so that sum of them is equal to 1
        @param entries original entries
        @param num_entries the number of entries
     */
    inline void normalize_arr(RealNumType* const entries, const int num_entries)
    {
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
     */
    void convert_range(const char *str, int &lower, int &upper, int &step_size);

    /**
        Convert a string to to range lower:upper:step_size with error checking
        @param str original string
        @param lower (OUT) lower bound of the range
        @param upper (OUT) upper bound of the range
        @param step_size (OUT) step size of the range
     */
    void convert_range(const char *str, RealNumType &lower, RealNumType &upper, RealNumType &step_size);

    /**
        Reinitialize an array of real number (RealNumType*)
        @param arr the input RealNumType*
        @param size the size of the input array
        @param delete_first TRUE to delete the current array before reinitializing
        @param set_zero TRUE to initialize all new items at 0
     */
    void reinitDoubleArr(RealNumType* &arr, StateType size, bool delete_first = true, bool set_zero = true);

    /**
        Convert char* into vector of strings separated by separator
     */
    void convert_string_vec(const char *str, StrVector &str_vec, char separator = ',');

    /**
        Convert string to PositionType, with error checking
        @param str original string
        @return the number
     */
    PositionType convert_positiontype(const char *str);

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
    void parseArg(int argc, char *argv[], Params &params);

    /**
        Init default parameters
        @param params (OUT) program parameters
     */
    void initDefaultValue(Params &params);

    /**
        Init default parameters with a random seed
        @param params (OUT) program parameters
        @param params (OUT) program parameters
     */
    void initDefaultValue(Params &params, uint64_t n_random_seed);

    /**
        Show quick start guide
     */
    void quickStartGuide();

    /**
        Detect the format of input file
        @param input_file file name
        @return
            IN_NEWICK if file in newick format,
            IN_NEXUS if in nexus format,
            IN_FASTA if in fasta format,
            IN_PHYLIP if in phylip format,
            IN_COUNTSFILE if in counts format (PoMo),
            IN_OTHER if file format unknown.
     */
    InputType detectInputFile(const char *input_file);

    /**
        If file exists, ask user to overwrite it or not
        @param filename file name
        @return TRUE if agree to overwrite an existing file, or simply file does not exist
     */
    bool overwriteFile(char *filename);

    /**
        Remove white space at the beginning and end of the string
        @param str (IN/OUT) string to be trimmed
    */
    void trimString(std::string &str);

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
    template<class T1, class T2>
    void quicksort(T1* arr, int left, int right, T2* arr2 = NULL) {
        if (left > right) return;
        ASSERT(left <= right);
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
     * Print usage for CMaple
     * @param program arguments list
     * @param full_command TRUE to print all available commands, FALSE to print normal usage dialog
     */
    void usage_cmaple(char* argv[], bool full_command);

    /**
     * Print copyright
     */
    void printCopyright(std::ostream &out);
}

#endif
