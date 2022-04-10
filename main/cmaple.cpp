/*
 *  cmaple.h
 *  Created on: Jan 25, 2022
 *      Author: Nhan Ly-Trong
 */

#include "cmaple.h"
#include "utils/timeutil.h"

CMaple::CMaple()
{
    params = NULL;
    aln = new Alignment();
}

CMaple::CMaple(Params *n_params)
{
    params = n_params;
    aln = new Alignment();
}

CMaple::~CMaple()
{
    if (aln)
    {
        delete aln;
        aln = NULL;
    }
}

void CMaple::extractDiffFile()
{
    // record the starting time
    auto start = getRealTime();
   
    // read input sequences
    ASSERT(params->aln_path);
    StrVector sequences;
    StrVector seq_names;
    aln->readSequences(params->aln_path, sequences, seq_names);
    
    // validate the input sequences
    if (sequences.size() == 0)
        outError("Empty input sequences. Please check and try again!");
    // make sure all sequences have the same length
    if (detectInputFile(params->aln_path) == IN_FASTA)
        for (int i = 0; i < sequences.size(); i++)
        {
            if (sequences[i].length() != sequences[0].length())
                outError("Sequence " + seq_names[i] + " has a different length compared to the first sequence.");
        }
    
    // generate reference sequence from the input sequences
    string ref_sequence;
    // read the reference sequence from file (if the user supplies it)
    if (params->ref_path)
        ref_sequence = aln->readRef(params->ref_path, params->only_extract_diff);
    else
        ref_sequence = aln->generateRef(sequences, seq_names, params->only_extract_diff);
    
    // prepare output (Diff) file
    // init diff_path if it's null
    if (!params->diff_path)
    {
        string diff_path_str(params->aln_path);
        diff_path_str += ".diff";
        
        params->diff_path = new char[diff_path_str.length() + 1];
        strcpy(params->diff_path, diff_path_str.c_str());
    }
    // ask for overwriting existing file (if necessary)
    if (fileExists(params->diff_path) && !overwriteFile(params->diff_path))
        exit(EXIT_SUCCESS);
    // open the output file
    ofstream out = ofstream(params->diff_path);
    
    // write reference sequence to the output file
    out << ">" << REF_NAME << endl;
    out << ref_sequence << endl;
    
    // extract and write mutations of each sequence to file
    aln->extractMutations(sequences, seq_names, ref_sequence, out, params->only_extract_diff);
    
    // close the output file
    out.close();
    
    // record the end time and show the runtime
    auto end = getRealTime();
    cout << "The input alignment is converted into DIFF format at " << params->diff_path << endl;
    cout << " - Converting time: " << end-start << endl;
}

void runCMaple(Params &params)
{
    
    CMaple cmaple = CMaple(&params);
    
    // extract sequences (in vectors of regions) from an input alignment (in PHYLIP or FASTA format)
    if (params.aln_path)
    {
        cmaple.extractDiffFile();
        
        // stop if the user only wants to extract diff file
        if (params.only_extract_diff)
            return;
        
        // convert Sequences (from vector of Mutations into vector Regions) for further inference
        cmaple.aln->convertSequences();
    }
    // or read sequences (in vectors of regions) from a DIFF file
    else
        cmaple.aln->readDiff(params.diff_path, params.ref_path);
    
    // sort sequences by their distances to the reference sequence
    cmaple.aln->sortSeqsByDistances(params.hamming_weight);
    
    cout << "dsfds" << endl;
    
}
