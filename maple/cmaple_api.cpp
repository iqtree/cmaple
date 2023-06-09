#include "cmaple_api.h"
#include "../utils/tools.h"
#include "../maple/cmaple.h"

using namespace cmaple;

/**
    Initialize library instance with default settings
 */
void __attribute__ ((constructor)) initLibrary(void)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
}

/**
    clean up the library instance
 */
void __attribute__ ((destructor)) cleanUpLibrary(void)
{
    // do nothing
}

void cmaple::printCMapleCopyright(std::ostream &out)
{
    printCopyright(out);
}

void cmaple::setAln(char* aln)
{
    Params& params = Params::getInstance();
    
    // specify the input file
    params.aln_path = aln;
}

void cmaple::setModel(std::string& model)
{
    Params& params = Params::getInstance();
    
    // specify the model
    params.model_name = model;
}

void cmaple::runInference()
{
    Params& params = Params::getInstance();
    
    // run CMaple
    runCMaple(params);
}

void cmaple::inferTree(char* input)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
    
    // specify the input file
    params.aln_path = input;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::inferTree(char* input, std::string& model)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
    
    // specify the input file
    params.aln_path = input;
    
    // specify the model
    params.model_name = model;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::inferTree(char* input, Params& params)
{
    // specify the input file
    params.aln_path = input;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::extractDiff(char* input_aln)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    
    // specify the input file
    params.aln_path = input_aln;
    
    // set a flag to only extract the DIFF file from an input alignment
    params.only_extract_diff = true;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::extractDiff(char* input_aln, char* ref)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    
    // specify the input file
    params.aln_path = input_aln;
    
    // specify the ref
    params.ref_path = ref;
    
    // set a flag to only extract the DIFF file from an input alignment
    params.only_extract_diff = true;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::reconstructAln(char* input_diff)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    
    // specify the input file
    params.aln_path = input_diff;
    
    // Init the path to the output alignment
    std::string output_aln_str(input_diff);
    output_aln_str += ".fa";
    params.output_aln = new char[output_aln_str.length() + 1];
    strcpy(params.output_aln, output_aln_str.c_str());
    
    // run CMaple
    runCMaple(params);
}

void cmaple::extendTree(char* tree, char* input)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
    
    // specify the input file
    params.aln_path = input;
    
    // specify the input tree
    params.input_treefile = tree;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::extendTree(char* tree, char* input, std::string& model)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
    
    // specify the input file
    params.aln_path = input;
    
    // specify the input tree
    params.input_treefile = tree;
    
    // specify the model
    params.model_name = model;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::extendTree(char* tree, char* input, Params& params)
{
    // specify the input file
    params.aln_path = input;
    
    // specify the input tree
    params.input_treefile = tree;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::computeBranchSupports(char* input, uint32_t num_threads, int64_t random_seed)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
    
    // set random seed (if provided)
    if (random_seed != -1) params.ran_seed = random_seed;
    cmaple::init_random(params.ran_seed, true);
    
    // set number of threads
    params.num_threads = num_threads;
    
    // specify the input file
    params.aln_path = input;
    
    // turn on the flag to compute aLRT-SH
    params.compute_aLRT_SH = true;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::computeBranchSupports(char* tree, char* input, uint32_t num_threads, int64_t random_seed)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
    
    // set random seed (if provided)
    if (random_seed != -1) params.ran_seed = random_seed;
    cmaple::init_random(params.ran_seed, true);
    
    // set number of threads
    params.num_threads = num_threads;
    
    // specify the input file
    params.aln_path = input;
    
    // specify the input tree
    params.input_treefile = tree;
    
    // turn on the flag to compute aLRT-SH
    params.compute_aLRT_SH = true;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::computeBranchSupports(char* tree, char* input, std::string& model, uint32_t num_threads, int64_t random_seed)
{
    // init params ~ inference settings
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.overwrite_output = true;
    
    // set random seed (if provided)
    if (random_seed != -1) params.ran_seed = random_seed;
    cmaple::init_random(params.ran_seed, true);
    
    // set number of threads
    params.num_threads = num_threads;
    
    // specify the input file
    params.aln_path = input;
    
    // specify the input tree
    params.input_treefile = tree;
    
    // specify the model
    params.model_name = model;
    
    // turn on the flag to compute aLRT-SH
    params.compute_aLRT_SH = true;
    
    // run CMaple
    runCMaple(params);
}

void cmaple::computeBranchSupports(char* tree, char* input, Params& params, uint32_t num_threads, int64_t random_seed)
{
    // set random seed (if provided)
    if (random_seed != -1) params.ran_seed = random_seed;
    cmaple::init_random(params.ran_seed, true);
    
    // set number of threads
    params.num_threads = num_threads;
    
    // specify the input file
    params.aln_path = input;
    
    // specify the input tree
    params.input_treefile = tree;
    
    // turn on the flag to compute aLRT-SH
    params.compute_aLRT_SH = true;
    
    // run CMaple
    runCMaple(params);
}
