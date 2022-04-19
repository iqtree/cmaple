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
    model = new Model();
    
    cumulative_rate = NULL;
    pseu_mutation_count = NULL;
}

CMaple::CMaple(Params *n_params)
{
    params = n_params;
    aln = new Alignment();
    model = new Model();
    
    cumulative_rate = NULL;
    pseu_mutation_count = NULL;
}

CMaple::~CMaple()
{
    if (aln)
    {
        delete aln;
        aln = NULL;
    }
    
    if (cumulative_rate)
    {
        delete[] cumulative_rate;
        cumulative_rate = NULL;
    }
    
    if (model)
    {
        delete model;
        model = NULL;
    }
    
    if (pseu_mutation_count)
    {
        delete[] pseu_mutation_count;
        pseu_mutation_count = NULL;
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
        for (PositionType i = 0; i < sequences.size(); i++)
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

void CMaple::loadInput()
{
    // extract sequences (in vectors of regions) from an input alignment (in PHYLIP or FASTA format)
    if (params->aln_path)
    {
        extractDiffFile();
        
        // convert Sequences (from vector of Mutations into vector Regions) for further inference
        if (!params->only_extract_diff)
            aln->convertSequences();
    }
    // or read sequences (in vectors of regions) from a DIFF file
    else
    {
        if (params->only_extract_diff)
            outError("To export a Diff file, please supple an alignment via --aln <ALIGNMENT>");
        
        aln->readDiff(params->diff_path, params->ref_path);
    }
}

void CMaple::computeCumulativeRate()
{
    PositionType sequence_length = aln->ref_seq.size();
    ASSERT(sequence_length > 0);
    cumulative_rate = new double[sequence_length];
    
    cumulative_rate[0] = model->mutation_mat[aln->ref_seq[0] * (aln->num_states + 1)];
    
    for (PositionType i = 1; i < sequence_length; i++)
        cumulative_rate[i] = cumulative_rate[i - 1] + model->mutation_mat[aln->ref_seq[i] * (aln->num_states + 1)];
        
}

void CMaple::preInference()
{
    // sort sequences by their distances to the reference sequence
    aln->sortSeqsByDistances(params->hamming_weight);
    
    // extract related info (freqs, log_freqs) of the ref sequence
    model->extractRefInfo(aln->ref_seq, aln->num_states);
    
    // init the mutation matrix from a model name
    model->initMutationMat(params->model_name, aln->num_states, pseu_mutation_count);
    
    // compute cummulative rates of the ref sequence
    computeCumulativeRate();
}

void CMaple::tmpTestingMethod()
{
    // show a comparison matrix
    for (PositionType i = 0; i < aln->size(); i++)
        cout << "\t" << aln->data()[i]->seq_name;
    for (PositionType i = 0; i < aln->size(); i++)
    {
        Sequence* seq1 = aln->data()[i];
        cout << endl;
        cout << seq1->seq_name << "\t";
        
        for (PositionType j = 0; j < aln->size(); j++)
        {
            Sequence* seq2 = aln->data()[j];
            cout << aln->compareSequences(seq1->mutations, seq2->mutations) << "\t";
        }
    }
    cout << endl;
}

void runCMaple(Params &params)
{
    auto start = getRealTime();
    CMaple cmaple = CMaple(&params);
    
    // load input data
    cmaple.loadInput();
    
    // terminate if the user only wants to export a Diff file from an alignment
    if (params.only_extract_diff)
        return;
    
    // prepare for the inference
    cmaple.preInference();
    
    // just test new method
    cmaple.tmpTestingMethod();
    
    // show runtime
    auto end = getRealTime();
    cout << "Runtime: " << end - start << "s" << endl;
}
