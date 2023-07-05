//
//  alignment.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "alignmentbase.h"

using namespace std;
using namespace cmaple;

char cmaple::symbols_protein[] = "ARNDCQEGHILKMFPSTWYVX"; // X for unknown AA
char cmaple::symbols_dna[]     = "ACGT";
char cmaple::symbols_rna[]     = "ACGU";
char cmaple::symbols_morph[] = "0123456789ABCDEFGHIJKLMNOPQRSTUV";

cmaple::AlignmentBase::AlignmentBase() = default;

cmaple::AlignmentBase::~AlignmentBase() = default;

void cmaple::AlignmentBase::processSeq(string &sequence, string &line, PositionType line_num) {
    for (string::iterator it = line.begin(); it != line.end(); ++it) {
        if ((*it) <= ' ') continue;
        if (isalnum(*it) || (*it) == '-' || (*it) == '?'|| (*it) == '.' || (*it) == '*' || (*it) == '~')
            sequence.append(1, toupper(*it));
        else if (*it == '(' || *it == '{') {
            auto start_it = it;
            while (*it != ')' && *it != '}' && it != line.end())
                ++it;
            if (it == line.end())
                throw "Line " + convertIntToString(line_num) + ": No matching close-bracket ) or } found";
            sequence.append(1, '?');
            cout << "NOTE: Line " << line_num << ": " << line.substr(start_it-line.begin(), (it-start_it)+1) << " is treated as unknown character" << endl;
        } else {
            throw "Line " + convertIntToString(line_num) + ": Unrecognized character "  + *it;
        }
    }
}

void cmaple::AlignmentBase::readFasta(const char *aln_path, StrVector &sequences, StrVector &seq_names, bool check_min_seqs){
    ostringstream err_str;
    igzstream in;
    PositionType line_num = 1;
    string line;

    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(aln_path);
    // remove the failbit
    in.exceptions(ios::badbit);

    {
        cout << "Reading FASTA file" <<endl;
        
        for (; !in.eof(); ++line_num) {
            safeGetline(in, line);
            if (line == "") {
                continue;
            }

            if (line[0] == '>') { // next sequence
                string::size_type pos = line.find_first_of("\n\r");
                seq_names.push_back(line.substr(1, pos-1));
                trimString(seq_names.back());
                sequences.push_back("");
                continue;
            }
            
            // read sequence contents
            if (sequences.empty()) {
                throw "First line must begin with '>' to define sequence name";
            }
            
            processSeq(sequences.back(), line, line_num);
        }
    }

    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();
    
    if (sequences.size() < MIN_NUM_TAXA && check_min_seqs)
        throw "There must be at least " + convertIntToString(MIN_NUM_TAXA) + " sequences";

    // now try to cut down sequence name if possible
    PositionType i, step = 0;
    StrVector new_seq_names, remain_seq_names;
    new_seq_names.resize(seq_names.size());
    remain_seq_names = seq_names;

    for (step = 0; step < 4; ++step) {
        bool duplicated = false;
        unordered_set<string> namesSeenThisTime;
        //Set of shorted names seen so far, this iteration
        for (i = 0; i < (PositionType) seq_names.size(); ++i) {
            if (remain_seq_names[i].empty()) continue;
            size_t pos = remain_seq_names[i].find_first_of(" \t");
            if (pos == string::npos) {
                new_seq_names[i] += remain_seq_names[i];
                remain_seq_names[i] = "";
            } else {
                new_seq_names[i] += remain_seq_names[i].substr(0, pos);
                remain_seq_names[i] = "_" + remain_seq_names[i].substr(pos+1);
            }
            if (!duplicated) {
                //add the shortened name for sequence i to the
                //set of shortened names seen so far, and set
                //duplicated to true if it was already there.
                duplicated = !namesSeenThisTime.insert(new_seq_names[i]).second;
            }
        }
        if (!duplicated) break;
    }
    
    if (step > 0) {
        for (i = 0; i < (PositionType) seq_names.size(); ++i)
            if (seq_names[i] != new_seq_names[i]) {
                cout << "NOTE: Change sequence name '" << seq_names[i] << "' -> " << new_seq_names[i] << endl;
            }
    }

    seq_names = new_seq_names;
}

void cmaple::AlignmentBase::readPhylip(const char *aln_path, StrVector &sequences, StrVector &seq_names, bool check_min_seqs)
{
    ostringstream err_str;
    igzstream in;
    PositionType line_num = 1;
    
    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(aln_path);
    PositionType seq_id = 0;
    PositionType nseq = 0;
    PositionType nsite = 0;
    string line;
    
    // remove the failbit
    in.exceptions(ios::badbit);

    for (; !in.eof(); ++line_num) {
        safeGetline(in, line);
        line = line.substr(0, line.find_first_of("\n\r"));
        if (line == "") continue;

        if (nseq == 0) { // read number of sequences and sites
            istringstream line_in(line);
            if (!(line_in >> nseq >> nsite))
                throw "Invalid PHYLIP format. First line must contain number of sequences and sites";
           
            if (nseq < MIN_NUM_TAXA && check_min_seqs)
                throw "There must be at least " + convertIntToString(MIN_NUM_TAXA) + " sequences";
            if (nsite < 1)
                throw "No alignment columns";

            seq_names.resize(nseq, "");
            sequences.resize(nseq, "");

        } else { // read sequence contents
            if (seq_names[seq_id] == "") { // cut out the sequence name
                string::size_type pos = line.find_first_of(" \t");
                if (pos == string::npos) pos = 10; //  assume standard phylip
                seq_names[seq_id] = line.substr(0, pos);
                line.erase(0, pos);
            }
            
            PositionType old_len = sequences[seq_id].length();
            processSeq(sequences[seq_id], line, line_num);
            
            if (sequences[seq_id].length() != sequences[0].length()) {
                err_str << "Line " << line_num << ": Sequence " << seq_names[seq_id] << " has wrong sequence length " << sequences[seq_id].length() << endl;
                throw err_str.str();
            }
            if ((PositionType) sequences[seq_id].length() > old_len)
                ++seq_id;
            if (seq_id == nseq) {
                seq_id = 0;
            }
        }
    }
    
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();
}

void cmaple::AlignmentBase::readSequences(const char* aln_path, StrVector &sequences, StrVector &seq_names, bool check_min_seqs)
{
    // detect the input file format
    InputType intype = detectInputFile(aln_path);
    
    // read the input file
    cout << "Reading alignment file " << aln_path << " ... ";
    switch (intype) {
        case IN_FASTA:
            cout << "Fasta format detected" << endl;
            readFasta(aln_path, sequences, seq_names, check_min_seqs);
            break;
            
        case IN_PHYLIP:
            cout << "Phylip format detected" << endl;
            readPhylip(aln_path, sequences, seq_names, check_min_seqs);
            break;
            
        default:
            outError("Please input an alignment file in FASTA or PHYLIP format!");
            break;
    }
}

string cmaple::AlignmentBase::generateRef(StrVector &sequences)
{
    // validate the input sequences
    if (!sequences.size() || !sequences[0].length())
        outError("Empty input sequences. Please check & try again!");
    
    cout << "Generating a reference sequence from the input alignment..." << endl;
    
    // init dummy variables
    char NULL_CHAR = '\0';
    string ref_str (sequences[0].length(), NULL_CHAR);
    
    // determine a character for each site one by one
    PositionType threshold = sequences.size() * 0.5;
    for (PositionType i = 0; i < (PositionType) ref_str.length(); ++i)
    {
        // Init a map to count the number of times each character appears
        std::map<char, PositionType> num_appear;
        
        for (PositionType j = 0; j < (PositionType) sequences.size(); ++j)
        {
            // update num_appear for the current character
            PositionType count = num_appear[sequences[j][i]] + 1;
            num_appear[sequences[j][i]] = count;
            
            // stop counting if a character appear in more than 1/2 sequences at the current site
            if (count >= threshold)
            {
                ref_str[i] = sequences[j][i];
                break;
            }
        }
        
        // manually determine the most popular charater for the current site (if no character dominates all the others)
        if (ref_str[i] == NULL_CHAR)
        {
            for (const std::pair<const char, PositionType>& character : num_appear)
                if (ref_str[i] == NULL_CHAR || character.second > num_appear[ref_str[i]])
                    ref_str[i] = character.first;
        }
    }
    
    // parse ref_sequence into vector of states
    parseRefSeq(ref_str);
    
    // return the reference genome
    return ref_str;
}

string cmaple::AlignmentBase::readRef(const std::string& ref_path)
{
    ASSERT(ref_path.length());
    if (!fileExists(ref_path))
        outError("File not found " + ref_path);
    
    string ref_str = "";
    
    // read sequences from file
    cout << "Reading a reference sequence from file..." << endl;
    StrVector str_sequences;
    StrVector seq_names;
    readSequences(ref_path.c_str(), str_sequences, seq_names, false);
    
    // validate the input sequence(s)
    if (!str_sequences.size() || !str_sequences[0].length())
        outError("No sequence found for the reference!");
    
    // extract the ref_sequence
    ref_str = str_sequences[0];
    
    // parse ref_str into vector of states (if necessary)
    parseRefSeq(ref_str);
    
    return ref_str;
}

void cmaple::AlignmentBase::outputMutation(ofstream &out, Sequence* sequence, char state_char, PositionType pos, PositionType length)
{
    // output the mutation into a MAPLE file
    out << state_char << "\t" << (pos + 1);
    if (length != -1)
        out << "\t" << length;
    out << endl;
    
    // add the mutation into sequence
    if (sequence)
    {
        if (length == -1)
            sequence->emplace_back(convertChar2State(state_char), pos);
        else
            sequence->emplace_back(convertChar2State(state_char), pos, length);
    }
}

void cmaple::AlignmentBase::extractMutations(StrVector &str_sequences, StrVector &seq_names, string ref_sequence, ofstream &out, bool only_extract_maple)
{
    ASSERT(str_sequences.size() == seq_names.size() && str_sequences.size() > 0 && out);
    data.clear();
    Sequence* sequence = NULL;
    PositionType seq_length = ref_sequence.length();
    
    // extract mutations of sequences one by one
    for (PositionType i = 0; i < (PositionType) str_sequences.size(); ++i)
    {
        // validate the sequence length
        string str_sequence = str_sequences[i];
        if (seq_length != (PositionType) str_sequence.length())
            outError("The sequence length of " + seq_names[i] + " (" + convertIntToString(str_sequence.length()) + ") is different from that of the reference sequence (" + convertIntToString(ref_sequence.length()) + ")!");
        
        // write taxon name
        out << ">" << seq_names[i] << endl;
        
        // init new sequence instance for the inference process afterwards
        if (!only_extract_maple)
        {
            data.emplace_back(seq_names[i]);
            sequence = &data.back();
        }
        
        // init dummy variables
        int state = 0;
        PositionType length = 0;
        for (PositionType pos = 0; pos < seq_length; ++pos)
        {
            switch (state)
            {
                case 0: // previous character is neither 'N' nor '-'
                    if (str_sequence[pos] != ref_sequence[pos])
                    {
                        length = 1;
                        
                        // starting a sequence of 'N'
                        if (toupper(str_sequence[pos]) == 'N' && getSeqType() == SEQ_DNA)
                            state = 1;
                        // starting a sequence of '-'
                        else if (str_sequence[pos] == '-')
                            state = 2;
                        // output a mutation
                        else
                            outputMutation(out, sequence, str_sequence[pos], pos);
                    }
                    break;
                case 1: // previous character is 'N'
                    // inscrease the length if the current character is still 'N'
                    if (toupper(str_sequence[pos]) == 'N' && str_sequence[pos] != ref_sequence[pos])
                        ++length;
                    else
                    {
                        // output the previous sequence of 'N'
                        outputMutation(out, sequence, str_sequence[pos-1], pos - length, length);
                        
                        // reset state
                        state = 0;
                        
                        // handle new character different from the reference
                        if (str_sequence[pos] != ref_sequence[pos])
                        {
                            length = 1;
                            // starting a sequence of '-'
                            if (str_sequence[pos] == '-')
                                state = 2;
                            // output a mutation
                            else
                            {
                                outputMutation(out, sequence, str_sequence[pos], pos);
                                state = 0;
                            }
                        }
                    }
                    break;
                case 2: // previous character is '-'
                    // inscrease the length if the current character is still '-'
                    if (toupper(str_sequence[pos]) == '-' && str_sequence[pos] != ref_sequence[pos])
                        ++length;
                    else
                    {
                        // output the previous sequence of '-'
                        outputMutation(out, sequence, str_sequence[pos-1], pos - length, length);
                        
                        // reset state
                        state = 0;
                        
                        // handle new character different from the reference
                        if (str_sequence[pos] != ref_sequence[pos])
                        {
                            length = 1;
                            // starting a sequence of 'N'
                            if (toupper(str_sequence[pos]) == 'N' && getSeqType() == SEQ_DNA)
                                state = 1;
                            // output a mutation
                            else
                            {
                                outputMutation(out, sequence, str_sequence[pos], pos);
                                state = 0;
                            }
                        }
                    }
                    break;
            }
        }
        
        //  output the last sequence of 'N' or '-' (if any)
        if (state != 0)
            outputMutation(out, sequence, str_sequence[str_sequence.length() - 1], str_sequence.length() - length, length);
    }
}

void cmaple::AlignmentBase::parseRefSeq(string& ref_sequence)
{
    ref_seq.resize(ref_sequence.length());
    
    // transform ref_sequence to uppercase
    transform(ref_sequence.begin(), ref_sequence.end(), ref_sequence.begin(), ::toupper);
    
    for (PositionType i = 0; i < (PositionType) ref_sequence.length(); ++i)
    {
        ref_seq[i] = convertChar2State(ref_sequence[i]);
        
        // validate the ref state
        if (ref_seq[i] >= num_states)
        {
            ref_sequence[i] = convertState2Char(0);
            outWarning("Invalid reference state found at site " + convertPosTypeToString(i) + " was replaced by a default state " + ref_sequence[i]);
            ref_seq[i] = 0;
        }
    }
}

void cmaple::AlignmentBase::readMapleFile(const std::string& aln_filename, const std::string& ref_path)
{
    ASSERT(aln_filename.length());
    
    if (!fileExists(aln_filename))
        outError("File not found " + aln_filename);
    
    // read the reference sequence (if the user supplies it)
    string ref_sequence{};
    if (ref_path.length())
        ref_sequence = readRef(ref_path);
    
    // init dummy variables
    string seq_name = "";
    vector<Mutation> mutations;
    ifstream in = ifstream(aln_filename);
    PositionType line_num = 1;
    string line;

    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    // remove the failbit
    in.exceptions(ios::badbit);

    cout << "Reading a MAPLE file" << endl;
    
    // extract reference sequence first
    for (; !in.eof(); ++line_num)
    {
        safeGetline(in, line);
        if (line == "") continue;
        
        // read the first line (">REF")
        if (line[0] == '>')
        {
            string::size_type pos = line.find_first_of("\n\r");
            seq_name = line.substr(1, pos-1);
            
            // transform seq_name to upper case
            transform(seq_name.begin(), seq_name.end(), seq_name.begin(), ::toupper);
            
            if (seq_name != REF_NAME && seq_name != "REFERENCE")
                outError("MAPLE file must start by >REF. Please check and try again!");
        }
        // read the reference sequence
        else
        {
            // make sure the first line was found
            if (seq_name != REF_NAME && seq_name != "REFERENCE")
                outError("MAPLE file must start by >REF. Please check and try again!");
            
            // choose the ref sequence if the user already supplies a reference sequence
            if (ref_path.length())
            {
                outWarning("Skipping the reference sequence in the MAPLE file since the reference sequence is already specified via '-ref' option.");
                line = ref_sequence;
            }
            
            // transform ref_sequence to uppercase
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            
            // detect the seq_type from the ref_sequences
            if (getSeqType() == SEQ_UNKNOWN)
            {
                StrVector tmp_str_vec;
                tmp_str_vec.push_back(line);
                setSeqType(detectSequenceType(tmp_str_vec));
            }
            
            // parse the reference sequence into vector of state
            parseRefSeq(line);
            
            // reset the seq_name
            seq_name = "";
            
            // break to read sequences of other taxa
            break;
        }
    }
    
    // extract sequences of other taxa one by one
    for (; !in.eof(); ++line_num)
    {
        safeGetline(in, line);
        if (line == "") continue;
        
        // Read sequence name
        if (line[0] == '>')
        {
            // record the sequence of the previous taxon
            if (seq_name.length())
            {
                data.emplace_back(seq_name, mutations);
                
                // reset dummy variables
                seq_name = "";
                mutations.clear();
            }
            
            // Read new sequence name
            string::size_type pos = line.find_first_of("\n\r");
            seq_name = line.substr(1, pos-1);
            if (!seq_name.length())
                outError("Empty sequence name found at line " + convertIntToString(line_num) + ". Please check and try again!");
        }
        // Read a Mutation
        else
        {
            // validate the input
            char separator = '\t';
            size_t num_items = std::count(line.begin(), line.end(), separator) + 1;
            if (num_items < 2 || num_items > 3)
                outError("Invalid input. Each difference must be presented be <Type>    <Position>  [<Length>]. Please check and try again!");
            
            // extract mutation info
            stringstream ssin(line);
            string tmp;
            
            // extract <Type>
            ssin >> tmp;
            StateType state = convertChar2State(toupper(tmp[0]));
            
            // extract <Position>
            ssin >> tmp;
            PositionType pos = convert_positiontype(tmp.c_str());
            if (pos <= 0 || pos > (PositionType) ref_seq.size())
                outError("<Position> must be greater than 0 and less than the reference sequence length (" + convertPosTypeToString(ref_seq.size()) + ")!");
            
            // extract <Length>
            PositionType length = 1;
            if (ssin.good())
            {
                ssin >> tmp;
                if (state == TYPE_N || state == TYPE_DEL)
                {
                    length = convert_positiontype(tmp.c_str());
                    if (length <= 0)
                        outError("<Length> must be greater than 0!");
                    if (length + pos - 1 > (PositionType) ref_seq.size())
                        outError("<Length> + <Position> must be less than the reference sequence length (" + convertPosTypeToString(ref_seq.size()) + ")!");
                }
                else
                    outWarning("Ignoring <Length> of " + tmp + ". <Length> is only appliable for 'N' or '-'.");
            }
            
            // add a new mutation into mutations
            if (state == TYPE_N || state == TYPE_DEL)
                mutations.emplace_back(state, pos - 1, length);
            else
                mutations.emplace_back(state, pos - 1);
        }
    }
    
    // Record the sequence of  the last taxon
    if (seq_name.length())
        data.emplace_back(seq_name, mutations);
    
    // validate the input
    if (ref_seq.size() == 0)
        outError("Reference sequence is not found!");
    if (data.size() < MIN_NUM_TAXA)
        outError("The number of taxa must be at least " + convertIntToString(MIN_NUM_TAXA));

    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();
}

void cmaple::AlignmentBase::reconstructAln(const std::string& aln_filename, const std::string& output_file, const cmaple::Params& params)
{
    ASSERT(aln_filename.length() && output_file.length());
    
    if (!fileExists(aln_filename))
        outError("File not found " + aln_filename);
    
    // check whether the output file already exists
    if (!params.overwrite_output && fileExists(output_file))
        outError("File " + output_file + " already exists. Use `--overwrite` option if you want overwrite it.\n");
    
    // init dummy variables
    string seq_name = "";
    PositionType current_pos = 1;
    PositionType seq_length = 0;
    ifstream in = ifstream(aln_filename);
    ofstream out = ofstream(output_file);
    PositionType line_num = 1;
    string line;
    string ref_str = "";

    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    // remove the failbit
    in.exceptions(ios::badbit);

    cout << "Reading a MAPLE file" << endl;
    
    // extract reference sequence first
    for (; !in.eof(); ++line_num)
    {
        safeGetline(in, line);
        if (line == "") continue;
        
        // read the first line (">REF")
        if (line[0] == '>')
        {
            string::size_type pos = line.find_first_of("\n\r");
            seq_name = line.substr(1, pos-1);
            
            // transform seq_name to upper case
            transform(seq_name.begin(), seq_name.end(), seq_name.begin(), ::toupper);
            
            if (seq_name != REF_NAME && seq_name != "REFERENCE")
                outError("MAPLE file must start by >REF. Please check and try again!");
        }
        // read the reference sequence
        else
        {
            // make sure the first line was found
            if (seq_name != REF_NAME && seq_name != "REFERENCE")
                outError("MAPLE file must start by >REF. Please check and try again!");
            
            // get ref_str
            ref_str = line;
            
            // get seq_length
            seq_length = ref_str.length();
            
            // reset the seq_name
            seq_name = "";
            
            // break to read sequences of other taxa
            break;
        }
    }
    
    // extract sequences of other taxa one by one
    for (; !in.eof(); ++line_num)
    {
        safeGetline(in, line);
        if (line == "") continue;
        
        // Read sequence name
        if (line[0] == '>')
        {
            // record the sequence of the previous taxon
            if (seq_name.length())
            {
                string tmp_str = "";
                if (current_pos <= seq_length)
                {
                    tmp_str.resize(seq_length - current_pos + 1, '-');
                    for (int i = 0; i < seq_length - current_pos + 1; ++i)
                        tmp_str[i] = ref_str[current_pos - 1 + i];
                }
                out << tmp_str << endl;
                
                // reset dummy variables
                seq_name = "";
                current_pos = 1;
            }
            
            // Read new sequence name
            string::size_type pos = line.find_first_of("\n\r");
            seq_name = line.substr(1, pos-1);
            if (!seq_name.length())
                outError("Empty sequence name found at line " + convertIntToString(line_num) + ". Please check and try again!");
            
            // write out the sequence name
            out << line << endl;
        }
        // Read a Mutation
        else
        {
            // validate the input
            char separator = '\t';
            size_t num_items = std::count(line.begin(), line.end(), separator) + 1;
            if (num_items < 2 || num_items > 3)
                outError("Invalid input. Each difference must be presented be <Type>    <Position>  [<Length>]. Please check and try again!");
            
            // extract mutation info
            stringstream ssin(line);
            string tmp;
            
            // extract <Type>
            ssin >> tmp;
            char state = toupper(tmp[0]);
            
            // extract <Position>
            ssin >> tmp;
            PositionType pos = convert_positiontype(tmp.c_str());
            if (pos <= 0 || pos > seq_length)
                outError("<Position> must be greater than 0 and less than the reference sequence length (" + convertPosTypeToString(ref_seq.size()) + ")!");
            
            // extract <Length>
            PositionType length = 1;
            if (ssin.good())
            {
                ssin >> tmp;
                if (state == '-' || state == 'N')
                {
                    length = convert_positiontype(tmp.c_str());
                    if (length <= 0)
                        outError("<Length> must be greater than 0!");
                    if (length + pos - 1 > seq_length)
                        outError("<Length> + <Position> must be less than the reference sequence length (" + convertPosTypeToString(seq_length) + ")!");
                }
                else
                    outWarning("Ignoring <Length> of " + tmp + ". <Length> is only appliable for 'N' or '-'.");
            }
            
            // add ref str if any
            if (current_pos < pos)
            {
                string tmp_str = "";
                tmp_str.resize(pos - current_pos, '-');
                for (int i = 0; i < pos - current_pos; ++i)
                    tmp_str[i] = ref_str[current_pos - 1 + i];
                out << tmp_str;
            }
            
            // add a new mutation
            string tmp_str = "";
            tmp_str.resize(length, '-');
            for (int i = 0; i < length; ++i)
                tmp_str[i] = state;
            out << tmp_str;
            
            // update current_pos
            current_pos = pos + length;
        }
    }
    
    // Record the sequence of  the last taxon
    string tmp_str = "";
    if (current_pos <= seq_length)
    {
        tmp_str.resize(seq_length - current_pos + 1, '-');
        for (int i = 0; i < seq_length - current_pos + 1; ++i)
            tmp_str[i] = ref_str[current_pos - 1 + i];
    }
    out << tmp_str << endl;
    
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();
    out.close();
}

char cmaple::AlignmentBase::convertState2Char(StateType state) {
    if (state == TYPE_N || state == TYPE_DEL) return '-';
    if (state > TYPE_INVALID) return '?';

    switch (getSeqType()) {
    case SEQ_BINARY:
        switch (state) {
        case 0:
            return '0';
        case 1:
            return '1';
        default:
            return '?'; // unrecognize character
        }
    case SEQ_DNA: // DNA
        switch (state) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 1+4+3:
            return 'R'; // A or G, Purine
        case 2+8+3:
            return 'Y'; // C or T, Pyrimidine
        case 1+8+3:
            return 'W'; // A or T, Weak
        case 2+4+3:
            return 'S'; // G or C, Strong
        case 1+2+3:
            return 'M'; // A or C, Amino
        case 4+8+3:
            return 'K'; // G or T, Keto
        case 2+4+8+3:
            return 'B'; // C or G or T
        case 1+2+8+3:
            return 'H'; // A or C or T
        case 1+4+8+3:
            return 'D'; // A or G or T
        case 1+2+4+3:
            return 'V'; // A or G or C
        default:
            return '?'; // unrecognize character
        }
        return state;
    case SEQ_PROTEIN: // Protein
        if (state < 20)
            return symbols_protein[(StateType)state];
        else if (state == 20) return 'B';
        else if (state == 21) return 'Z';
        else if (state == 22) return 'J';
//        else if (state == 4+8+19) return 'B';
//        else if (state == 32+64+19) return 'Z';
        else
            return '-';
    case SEQ_MORPH:
        // morphological state
        if (state < strlen(symbols_morph))
            return symbols_morph[state];
        else
            return '-';
    default:
        // unknown
        return '*';
    }
}

StateType cmaple::AlignmentBase::convertChar2State(char state) {
    if (state == '-')
        return TYPE_DEL;
    if (state == '?' || state == '.' || state == '~')
        return TYPE_N;

    char *loc;

    switch (getSeqType()) {
    case SEQ_BINARY:
        switch (state) {
        case '0':
            return 0;
        case '1':
            return 1;
        default:
            {
                string invalid_state_msg = "Invalid state ";
                invalid_state_msg += state;
                invalid_state_msg += ". Please check and try again!";
                outError(invalid_state_msg);
                return TYPE_INVALID;
            }
        }
        break;
    case SEQ_DNA: // DNA
        switch (state) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'U':
            return 3;
        case 'R':
            return 1+4+3; // A or G, Purine
        case 'Y':
            return 2+8+3; // C or T, Pyrimidine
        case 'O':
        case 'N':
        case 'X':
            return TYPE_N;
        case 'W':
            return 1+8+3; // A or T, Weak
        case 'S':
            return 2+4+3; // G or C, Strong
        case 'M':
            return 1+2+3; // A or C, Amino
        case 'K':
            return 4+8+3; // G or T, Keto
        case 'B':
            return 2+4+8+3; // C or G or T
        case 'H':
            return 1+2+8+3; // A or C or T
        case 'D':
            return 1+4+8+3; // A or G or T
        case 'V':
            return 1+2+4+3; // A or G or C
        default:
            {
                string invalid_state_msg = "Invalid state ";
                invalid_state_msg += state;
                invalid_state_msg += ". Please check and try again!";
                outError(invalid_state_msg);
                return TYPE_INVALID; // unrecognize character
            }
        }
        return state;
    case SEQ_PROTEIN: // Protein
//        if (state == 'B') return 4+8+19;
//        if (state == 'Z') return 32+64+19;
        if (state == 'B') return 20;
        if (state == 'Z') return 21;
        if (state == 'J') return 22;
        if (state == '*') return TYPE_N; // stop codon
        if (state == 'U') return TYPE_N; // 21st amino-acid
        if (state == 'O') return TYPE_N; // 22nd amino-acid
        loc = strchr(symbols_protein, state);

        if (!loc)
        {
            string invalid_state_msg = "Invalid state ";
            invalid_state_msg += state;
            invalid_state_msg += ". Please check and try again!";
            outError(invalid_state_msg);
            return TYPE_INVALID; // unrecognize character
        }
        state = loc - symbols_protein;
        if (state < 20)
            return state;
        else
            return TYPE_N;
    case SEQ_MORPH: // Standard morphological character
        loc = strchr(symbols_morph, state);

        if (!loc)
        {
            string invalid_state_msg = "Invalid state ";
            invalid_state_msg += state;
            invalid_state_msg += ". Please check and try again!";
            outError(invalid_state_msg);
            return TYPE_INVALID; // unrecognize character
        }
        state = loc - symbols_morph;
        return state;
    default:
        {
            string invalid_state_msg = "Invalid state ";
            invalid_state_msg += state;
            invalid_state_msg += ". Please check and try again!";
            outError(invalid_state_msg);
            return TYPE_INVALID;
        }
    }
}

PositionType cmaple::AlignmentBase::computeSeqDistance(Sequence& sequence, RealNumType hamming_weight)
{
    // dummy variables
    PositionType num_ambiguities = 0;
    PositionType num_diffs = 0;
    
    // browse mutations one by one
    for (const auto& mutation: sequence)
    {
        switch (mutation.type)
        {
            case TYPE_R: // Type R does not exist in Mutations (only in Regions)
                outError("Sorry! Something went wrong. Invalid mutation type (type R).");
                break;
            case TYPE_N: // handle unsequenced sites ('-' or 'N')
            case TYPE_DEL:
            case TYPE_O: // handle ambiguity
                num_ambiguities += mutation.getLength();
                ++num_diffs;
                break;
            default:
                // handle normal character
                if (mutation.type < num_states)
                    ++num_diffs;
                else
                {
                    ++num_diffs;
                    ++num_ambiguities;
                }
                break;
        }
    }
    
    // calculate and record the distance of the current sequence
    return num_diffs * hamming_weight + num_ambiguities;
}

void cmaple::AlignmentBase::sortSeqsByDistances(RealNumType hamming_weight)
{
   // init dummy variables
    PositionType num_seqs = data.size();
    PositionType *distances = new PositionType[num_seqs];
    PositionType *sequence_indexes = new PositionType[num_seqs];
    
    // calculate the distances of each sequence
    for (PositionType i = 0; i < num_seqs; ++i)
    {
        // dummy variables
        sequence_indexes[i] = i;
        
        // calculate and record the distance of the current sequence
        distances[i] = computeSeqDistance(data[i], hamming_weight);
        
        // NHANLT: debug
        //distances[i] *= 1000;
    }
    
    // NHANLT: debug

    
    // sort distances
    quicksort(distances, 0, num_seqs - 1, sequence_indexes);

    // re-order sequences by distances
    vector<Sequence> tmp_sequences(move(data));
    data.reserve(num_seqs);
    for (PositionType i = 0; i < num_seqs; ++i)
        data.push_back(move(tmp_sequences[sequence_indexes[i]]));
    
    // delete distances, sequence_indexes
    delete[] distances;
    delete[] sequence_indexes;
}

void cmaple::AlignmentBase::extractMapleFile(const std::string& aln_filename, const std::string& output_filename, const cmaple::Params& params, const bool only_extract_maple)
{   
    // read input sequences
    ASSERT(aln_filename.length());
    StrVector sequences;
    StrVector seq_names;
    readSequences(aln_filename.c_str(), sequences, seq_names);
    
    // validate the input sequences
    if (sequences.size() == 0)
        outError("Empty input sequences. Please check and try again!");
    // make sure all sequences have the same length
    if (params.aln_format == IN_FASTA || (params.aln_format == IN_UNKNOWN && detectInputFile(aln_filename.c_str()) == IN_FASTA))
        for (PositionType i = 0; i < (PositionType) sequences.size(); ++i)
        {
            if (sequences[i].length() != sequences[0].length())
                outError("Sequence " + seq_names[i] + " has a different length compared to the first sequence.");
        }
    
    // detect the type of the input sequences
    setSeqType(params.seq_type);
    if (getSeqType() == SEQ_UNKNOWN)
        setSeqType(detectSequenceType(sequences));
    
    // generate reference sequence from the input sequences
    string ref_sequence;
    // read the reference sequence from file (if the user supplies it)
    if (params.ref_path.length())
        ref_sequence = readRef(params.ref_path);
    else
        ref_sequence = generateRef(sequences);

    // check whether the output file already exists
    if (!params.overwrite_output && fileExists(output_filename))
        outError("File " + output_filename + " already exists. Use `--overwrite` option if you want overwrite it.\n");
    
    // open the output file
    ofstream out = ofstream(output_filename);
    
    // write reference sequence to the output file
    out << ">" << REF_NAME << endl;
    out << ref_sequence << endl;
    
    // extract and write mutations of each sequence to file
    extractMutations(sequences, seq_names, ref_sequence, out, only_extract_maple);
    
    // close the output file
    out.close();
}

SeqType cmaple::AlignmentBase::detectSequenceType(StrVector& sequences)
{
    size_t num_nuc   = 0;
    size_t num_ungap = 0;
    size_t num_bin   = 0;
    size_t num_alpha = 0;
    size_t num_digit = 0;
    double detectStart = getRealTime();
    size_t sequenceCount = sequences.size();

#pragma omp parallel for reduction(+:num_nuc,num_ungap,num_bin,num_alpha,num_digit)
    for (size_t seqNum = 0; seqNum < sequenceCount; ++seqNum) {
        auto start = sequences.at(seqNum).data();
        auto stop  = start + sequences.at(seqNum).size();
        for (auto i = start; i!=stop; ++i) {
            if ((*i) == 'A' || (*i) == 'C' || (*i) == 'G' || (*i) == 'T' || (*i) == 'U') {
                ++num_nuc;
                ++num_ungap;
                continue;
            }
            if ((*i)=='?' || (*i)=='-' || (*i) == '.' ) {
                continue;
            }
            if (*i != 'N' && *i != 'X' &&  (*i) != '~') {
                num_ungap++;
                if (isdigit(*i)) {
                    num_digit++;
                    if ((*i) == '0' || (*i) == '1') {
                        num_bin++;
                    }
                }
            }
            if (isalpha(*i)) {
                num_alpha++;
            }
        }
    } // end OMP parallel for

    if (verbose_mode >= VB_MED) {
        cout << "Sequence Type detection took " << (getRealTime()-detectStart) << " seconds." << endl;
    }
    if (((double)num_nuc) / num_ungap > 0.9)
    {
        std::cout << "DNA data detected." << std::endl;
        return SEQ_DNA;
    }
    if (((double)num_bin) / num_ungap > 0.9)
    {
        std::cout << "Binary data detected." << std::endl;
        return SEQ_BINARY;
    }
    if (((double)num_alpha + num_nuc) / num_ungap > 0.9)
    {
        std::cout << "Protein data detected." << std::endl;
        return SEQ_PROTEIN;
    }
    if (((double)(num_alpha + num_digit + num_nuc)) / num_ungap > 0.9)
    {
        std::cout << "Morphological data detected." << std::endl;
        return SEQ_MORPH;
    }
    return SEQ_UNKNOWN;
}

void cmaple::AlignmentBase::updateNumStates()
{
    switch (getSeqType()) {
        case SEQ_PROTEIN:
            num_states = 20;
            break;
            
        default: // dna
            num_states = 4;
            break;
    }
}

InputType cmaple::AlignmentBase::getAlignmentFormat(const std::string& n_format)
{
    string format(n_format);
    transform(format.begin(), format.end(), format.begin(), ::toupper);
    if (format == "MAPLE")
        return IN_MAPLE;
    if (format == "PHYLIP")
        return IN_PHYLIP;
    if (format == "FASTA")
        return IN_FASTA;
    
    // default
    return IN_UNKNOWN;
}

SeqType cmaple::AlignmentBase::getSeqType(const std::string& n_seqtype_str)
{
    string seqtype_str(n_seqtype_str);
    transform(seqtype_str.begin(), seqtype_str.end(), seqtype_str.begin(), ::toupper);
    if (seqtype_str == "DNA")
        return SEQ_DNA;
    if (seqtype_str == "AA")
        return SEQ_PROTEIN;
    
    // default
    return SEQ_UNKNOWN;
}