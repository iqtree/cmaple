//
//  alignment.cpp
//  alignment
//
//  Created by NhanLT on 31/3/2022.
//

#include "alignment.h"

using namespace std;
using namespace cmaple;

char cmaple::symbols_protein[] = "ARNDCQEGHILKMFPSTWYVX";  // X for unknown AA
char cmaple::symbols_dna[] = "ACGT";
char cmaple::symbols_rna[] = "ACGU";
char cmaple::symbols_morph[] = "0123456789ABCDEFGHIJKLMNOPQRSTUV";

cmaple::Alignment::Alignment()
    : data(std::vector<Sequence>()),
      ref_seq(std::vector<cmaple::StateType>()),
      num_states(4),
      attached_trees(std::unordered_set<void*>()){};

cmaple::Alignment::Alignment(std::istream& aln_stream,
                             const std::string& ref_seq,
                             const InputType format,
                             const cmaple::SeqRegion::SeqType seqtype)
    : Alignment() {
  // Initialize an alignment instance from the input stream
  read(aln_stream, ref_seq, format, seqtype);
}

cmaple::Alignment::Alignment(const std::string& aln_filename,
                             const std::string& ref_seq,
                             const InputType format,
                             const cmaple::SeqRegion::SeqType seqtype)
    : Alignment() {
  // Initialize an alignment instance from the input file
  read(aln_filename, ref_seq, format, seqtype);
}

cmaple::Alignment::~Alignment() = default;

void cmaple::Alignment::read(std::istream& aln_stream,
                             const std::string& ref_seq,
                             const InputType format,
                             const cmaple::SeqRegion::SeqType seqtype) {
  if (cmaple::verbose_mode >= cmaple::VB_MED) {
    std::cout << "Reading an alignment" << std::endl;
  }

  // Reset aln_base
  reset();

  // Set format (if specified)
  if (format != IN_AUTO && format != IN_UNKNOWN) {
    aln_format = format;
    // Otherwise, detect the format from the alignment
  } else {
    aln_format = detectInputFile(aln_stream);

    // validate the format
    if (aln_format == IN_UNKNOWN) {
      throw std::invalid_argument(
          "Failed to detect the format from the alignment!");
    }
  }

  // Set seqtype. If it's auto (not specified), we'll dectect it later when
  // reading the alignment
  if (seqtype != SeqRegion::SEQ_UNKNOWN) {
    setSeqType(seqtype);
  } else {
    setSeqType(SeqRegion::SEQ_AUTO);
  }

  // Read the input alignment
  try {
    // in FASTA or PHYLIP format
    if (aln_format != IN_MAPLE) {
      readFastaOrPhylip(aln_stream, ref_seq);
      // in MAPLE format
    } else {
      if (ref_seq.length() && cmaple::verbose_mode > cmaple::VB_QUIET) {
        outWarning(
            "Ignore the input reference as it must be already "
            "specified in the MAPLE format");
      }
      readMaple(aln_stream);
    }

    // sort sequences by their distances to the reference sequence
    sortSeqsByDistances();

    // avoid using DNA build for protein data
    if (NUM_STATES < num_states) {
      throw std::invalid_argument(
          "Look like you're using the wrong (DNA) compilation version for "
          "Protein data. Please use the version for Protein data (*-aa) "
          "instead.");
    }
  } catch (std::logic_error e) {
    throw std::invalid_argument(e.what());
  }
}

void cmaple::Alignment::read(const std::string& aln_filename,
                             const std::string& ref_seq,
                             const InputType format,
                             const cmaple::SeqRegion::SeqType seqtype) {
  if (!aln_filename.length()) {
    throw std::invalid_argument("Please specify an alignnment file");
  }
  assert(aln_filename.length() > 0);

  // Create a stream from the input alignment
  ifstream aln_stream;
  try {
    aln_stream.exceptions(ios::failbit | ios::badbit);
    aln_stream.open(aln_filename);
  } catch (ios::failure) {
    std::string err_msg(ERR_READ_INPUT);
    throw ios::failure(err_msg + aln_filename);
  }

  // Initialize an alignment instance from the input stream
  read(aln_stream, ref_seq, format, seqtype);

  // close aln_stream
  aln_stream.close();
}

void cmaple::Alignment::write(std::ostream& aln_stream,
                              const InputType& format) {
  assert(data.size() > 0);
  assert(format != IN_AUTO && format != IN_UNKNOWN);
    
  // Handle empty alignment
  if (!data.size()) {
    throw std::logic_error("Alignment is empty. Please call read(...) first!");
  }

  switch (format) {
    case IN_MAPLE:
      writeMAPLE(aln_stream);
      break;
    case IN_FASTA:
      writeFASTA(aln_stream);
      break;
    case IN_PHYLIP:
      writePHYLIP(aln_stream);
      break;
    default:  // IN_AUTO or IN_UNKNOWN
      throw std::invalid_argument(
          "Unsupported format for outputting the alignment!");
      break;
  }
}

void cmaple::Alignment::write(const std::string& aln_filename,
                              const InputType& format,
                              const bool overwrite) {
  // Validate the input
  if (!aln_filename.length()) {
    throw std::invalid_argument(
        "Please specify a filename to output the alignment");
  }

  assert(aln_filename.length() > 0);
  assert(data.size() > 0);
  assert(format != IN_AUTO && format != IN_UNKNOWN);
    
  // Check whether the output file already exists
  if (!overwrite && fileExists(aln_filename)) {
    throw ios::failure(
        "File " + aln_filename +
        " already exists. Please set overwrite = true to overwrite it.");
  }

  // Open a stream to write the output
  std::ofstream aln_stream = ofstream(aln_filename);

  // Write alignment to the stream
  write(aln_stream, format);

  // Close the stream
  aln_stream.close();
}

std::ostream& cmaple::operator<<(std::ostream& out_stream,
                                 cmaple::Alignment& aln) {
  assert(aln.data.size() > 0);
    
  // write the alignment to the stream in MAPLE format
  aln.write(out_stream);

  // return the stream
  return out_stream;
}

std::istream& cmaple::operator>>(std::istream& in_stream,
                                 cmaple::Alignment& aln) {
  // go back to the beginning og the stream
  in_stream.clear();
  in_stream.seekg(0);

  // read the alignment from a stream
  aln.read(in_stream);

  // return the stream
  return in_stream;
}

void cmaple::Alignment::reset() {
  setSeqType(cmaple::SeqRegion::SEQ_AUTO);
  data.clear();
  ref_seq.clear();
  aln_format = IN_AUTO;
  attached_trees.clear();
}

void cmaple::Alignment::processSeq(string& sequence,
                                   string& line,
                                   PositionType line_num) {
  for (string::iterator it = line.begin(); it != line.end(); ++it) {
    if ((*it) <= ' ') {
      continue;
    }
    if (isalnum(*it) || (*it) == '-' || (*it) == '?' || (*it) == '.' ||
        (*it) == '*' || (*it) == '~') {
      sequence.append(1, toupper(*it));
    } else if (*it == '(' || *it == '{') {
      auto start_it = it;
      while (*it != ')' && *it != '}' && it != line.end()) {
        ++it;
      }
      if (it == line.end()) {
        throw std::logic_error("Line " + convertIntToString(line_num) +
                               ": No matching close-bracket ) or } found");
      }
      sequence.append(1, '?');
      if (cmaple::verbose_mode > cmaple::VB_QUIET) {
        cout << "NOTE: Line " << line_num << ": "
             << line.substr((std::basic_string<char>::size_type) (start_it - line.begin()), (std::basic_string<char>::size_type) (it - start_it) + 1)
             << " is treated as unknown character" << endl;
      }
    } else {
      throw std::logic_error("Line " + convertIntToString(line_num) +
                             ": Unrecognized character " + *it);
    }
  }
}

void cmaple::Alignment::readFasta(std::istream& in,
                                  StrVector& sequences,
                                  StrVector& seq_names,
                                  bool check_min_seqs) {
  ostringstream err_str;
  PositionType line_num = 1;
  string line;

  // set the failbit and badbit
  in.exceptions(ios::failbit | ios::badbit);
  // remove the failbit
  in.exceptions(ios::badbit);

  {
    for (; !in.eof(); ++line_num) {
      safeGetline(in, line);
      if (line == "") {
        continue;
      }

      if (line[0] == '>') {  // next sequence
        string::size_type pos = line.find_first_of("\n\r");
        seq_names.push_back(line.substr(1, pos - 1));
        trimString(seq_names.back());
        sequences.push_back("");
        continue;
      }

      // read sequence contents
      if (sequences.empty()) {
        throw std::logic_error(
            "First line must begin with '>' to define sequence name");
      }

      processSeq(sequences.back(), line, line_num);
    }
  }

  // set the failbit again
  in.exceptions(ios::failbit | ios::badbit);
  // reset the stream
  resetStream(in);

  if (sequences.size() < MIN_NUM_TAXA && check_min_seqs) {
    throw std::logic_error("There must be at least " +
                           convertIntToString(MIN_NUM_TAXA) + " sequences");
  }

  // now try to cut down sequence name if possible
  std::vector<std::string>::size_type i;
  PositionType step = 0;
  StrVector new_seq_names, remain_seq_names;
  new_seq_names.resize(seq_names.size());
  remain_seq_names = seq_names;

  for (step = 0; step < 4; ++step) {
    bool duplicated = false;
    unordered_set<string> namesSeenThisTime;
    // Set of shorted names seen so far, this iteration
    for (i = 0; i < (PositionType)seq_names.size(); ++i) {
      if (remain_seq_names[i].empty()) {
        continue;
      }
      size_t pos = remain_seq_names[i].find_first_of(" \t");
      if (pos == string::npos) {
        new_seq_names[i] +=
          remain_seq_names[i];
        remain_seq_names[i] = "";
      } else {
        new_seq_names[i] +=
          remain_seq_names[i].substr(0, pos);
        remain_seq_names[i] = "_" +
          remain_seq_names[i].substr(pos + 1);
      }
      if (!duplicated) {
        // add the shortened name for sequence i to the
        // set of shortened names seen so far, and set
        // duplicated to true if it was already there.
        duplicated = !namesSeenThisTime.insert(new_seq_names[i]).second;
      }
    }
    if (!duplicated) {
      break;
    }
  }

  if (step > 0) {
    for (i = 0; i < seq_names.size(); ++i) {
      if (seq_names[i] != new_seq_names[i] &&
          cmaple::verbose_mode > cmaple::VB_QUIET) {
        cout << "NOTE: Change sequence name '" << seq_names[i] << "' -> "
             << new_seq_names[i] << endl;
      }
    }
  }

  seq_names = new_seq_names;
}

void cmaple::Alignment::readPhylip(std::istream& in,
                                   StrVector& sequences,
                                   StrVector& seq_names,
                                   bool check_min_seqs) {
  ostringstream err_str;
  PositionType line_num = 1;

  // set the failbit and badbit
  in.exceptions(ios::failbit | ios::badbit);
  std::vector<std::string>::size_type seq_id = 0;
  std::vector<std::string>::size_type nseq = 0;
  std::vector<std::string>::size_type nsite = 0;
  string line;

  // remove the failbit
  in.exceptions(ios::badbit);

  for (; !in.eof(); ++line_num) {
    safeGetline(in, line);
    line = line.substr(0, line.find_first_of("\n\r"));
    if (line == "") {
      continue;
    }

    if (nseq == 0) {  // read number of sequences and sites
      istringstream line_in(line);
      if (!(line_in >> nseq >> nsite)) {
        throw std::logic_error(
            "Invalid PHYLIP format. First line must contain number of "
            "sequences and sites");
      }

      if (nseq < MIN_NUM_TAXA && check_min_seqs) {
        throw std::logic_error("There must be at least " +
                               convertIntToString(MIN_NUM_TAXA) + " sequences");
      }
      if (nsite < 1) {
        throw std::logic_error("No alignment columns");
      }

      seq_names.resize(nseq, "");
      sequences.resize(nseq, "");

    } else {                          // read sequence contents
      if (seq_names[seq_id] == "") {  // cut out the sequence name
        string::size_type pos = line.find_first_of(" \t");
        if (pos == string::npos) {
          pos = 10;  //  assume standard phylip
        }
        seq_names[seq_id] = line.substr(0, pos);
        line.erase(0, pos);
      }

      std::basic_string<char>::size_type old_len = sequences[seq_id].length();
      processSeq(sequences[seq_id], line, line_num);

      if (sequences[seq_id].length() != sequences[0].length()) {
        err_str << "Line " << line_num << ": Sequence " << seq_names[seq_id]
                << " has wrong sequence length " << sequences[seq_id].length()
                << endl;
        throw std::logic_error(err_str.str());
      }
      if (sequences[seq_id].length() > old_len) {
        ++seq_id;
      }
      if (seq_id == nseq) {
        seq_id = 0;
      }
    }
  }

  // validate the number of seqs
  if (seq_id) {
    throw logic_error(
        "Number of sequences " + convertIntToString((int) nseq) +
        " (in the header) is different from the actual number of sequences " +
        convertIntToString((int) seq_id));
  }

  // set the failbit again
  in.exceptions(ios::failbit | ios::badbit);

  // reset the stream
  resetStream(in);
}

void cmaple::Alignment::readSequences(std::istream& aln_stream,
                                      StrVector& sequences,
                                      StrVector& seq_names,
                                      InputType n_aln_format,
                                      bool check_min_seqs) {
  // detect the input file format
  if (n_aln_format == IN_AUTO || n_aln_format == IN_UNKNOWN) {
    n_aln_format = detectInputFile(aln_stream);
  }
  assert(n_aln_format != IN_AUTO && n_aln_format != IN_UNKNOWN);

  // read the input file
  if (cmaple::verbose_mode >= cmaple::VB_MAX) {
    std::cout << "Reading alignment from a stream..." << std::endl;
  }
  switch (n_aln_format) {
    case IN_FASTA:
      readFasta(aln_stream, sequences, seq_names, check_min_seqs);
      break;

    case IN_PHYLIP:
      readPhylip(aln_stream, sequences, seq_names, check_min_seqs);
      break;

    default:
      throw std::logic_error(
          "Please input an alignment in FASTA or PHYLIP format!");
      break;
  }
}

auto cmaple::Alignment::generateRef(StrVector& sequences) -> string {
  assert(sequences.size() > 0);
  assert(sequences[0].length() > 0);
    
  // validate the input sequences
  if (!sequences.size() || !sequences[0].length()) {
    throw std::logic_error("Empty input sequences. Please check & try again!");
  }

  if (cmaple::verbose_mode >= cmaple::VB_MAX) {
    cout << "Generating a reference sequence from the input alignment..."
         << endl;
  }

  // init dummy variables
  const char NULL_CHAR = '\0';
  const char GAP = '-';
  string ref_str(sequences[0].length(), NULL_CHAR);
  const char DEFAULT_CHAR = cmaple::Alignment::convertState2Char(0, seq_type_);

  // determine a character for each site one by one
  PositionType threshold = (PositionType) (sequences.size() * 0.5);
  for (std::basic_string<char>::size_type i = 0; i < ref_str.length(); ++i) {
    // Init a map to count the number of times each character appears
    std::map<char, PositionType> num_appear;

    for (std::vector<std::string>::size_type j = 0; j < sequences.size(); ++j) {
      // update num_appear for the current character
      PositionType count = num_appear[sequences[j][i]] + 1;
      num_appear[sequences[j][i]] = count;

      // stop counting if a non-gap character appear in more than 1/2 sequences
      // at the current site
      if (count >= threshold && sequences[j][i] != GAP) {
        ref_str[i] = sequences[j][i];
        break;
      }
    }

    // manually determine the most popular charater for the current site (if no
    // character dominates all the others)
    if (ref_str[i] == NULL_CHAR) {
      for (const std::pair<const char, PositionType>& character : num_appear) {
        if (character.first != GAP &&
            (ref_str[i] == NULL_CHAR ||
             character.second > num_appear[ref_str[i]])) {
          ref_str[i] = character.first;
        }
      }
    }

    // if not found -> all characters in this site are gaps -> choose the
    // default state
    if (ref_str[i] == NULL_CHAR) {
      ref_str[i] = DEFAULT_CHAR;
    }
  }
    
  assert(ref_str.length() == sequences[0].length());

  // return the reference genome
  return ref_str;
}

auto cmaple::Alignment::readRefSeq(const std::string& ref_filename,
                                   const std::string& ref_name) -> string {
  if (!fileExists(ref_filename)) {
    throw ios::failure("File not found " + ref_filename);
  }
  if (!ref_name.length()) {
    throw std::invalid_argument(
        "Please specify the name of the reference sequence!");
  }
  assert(ref_filename.length() > 0);
  assert(ref_name.length() > 0);

  // convert ref_name to uppercase
  std::string ref_name_upcase(ref_name);
  transform(ref_name_upcase.begin(), ref_name_upcase.end(),
            ref_name_upcase.begin(), ::toupper);

  // read sequences from file
  if (cmaple::verbose_mode >= cmaple::VB_MAX) {
    cout << "Reading a reference sequence from an alignment file..." << endl;
  }
  StrVector str_sequences;
  StrVector seq_names;
  // Create a stream from the input alignment
  ifstream ref_stream;
  try {
    ref_stream.exceptions(ios::failbit | ios::badbit);
    ref_stream.open(ref_filename);
  } catch (ios::failure) {
    std::string err_msg(ERR_READ_INPUT);
    throw std::logic_error(err_msg + ref_filename);
  }

  // Read sequences from the alignment
  readSequences(ref_stream, str_sequences, seq_names, IN_AUTO, false);

  // close ref_stream
  ref_stream.close();

  // validate the input sequence(s)
  if (!str_sequences.size() || !str_sequences[0].length()) {
    throw std::logic_error("No sequence found for the reference!");
  }

  // extract the ref_sequence
  string ref_str = "";
  for (std::vector<std::string>::size_type i = 0; i < seq_names.size(); ++i) {
    std::string seq_name = seq_names[(std::vector<std::string>::size_type) i];
    transform(seq_name.begin(), seq_name.end(), seq_name.begin(), ::toupper);
    if (seq_name == ref_name_upcase) {
      ref_str = str_sequences[i];
      break;
    }
  }

    assert(ref_str.length() > 0);

  // Validate the output
  if (!ref_str.length()) {
    throw std::logic_error("The reference sequence named " + ref_name +
                           " is not found or empty");
  }

  return ref_str;
}

void cmaple::Alignment::addMutation(Sequence* sequence,
                                    char state_char,
                                    PositionType pos,
                                    PositionType length) {
  assert(sequence);
  assert(pos >= 0);
    
  // add the mutation into sequence
  if (length == -1) {
    sequence->emplace_back(convertChar2State(state_char), pos);
  } else {
    sequence->emplace_back(convertChar2State(state_char), pos, length);
  }
}

void cmaple::Alignment::extractMutations(const StrVector& str_sequences,
                                         const StrVector& seq_names,
                                         const string& ref_sequence) {
  // Validate the inputs
  assert(str_sequences.size() > 0);
  assert(str_sequences.size() == seq_names.size());
  assert(ref_sequence.length() > 0);
    
  if (!str_sequences.size()) {
    throw std::invalid_argument("The vector of sequences is empty");
  }

  data.clear();
  Sequence* sequence = nullptr;
  std::vector<char>::size_type seq_length = ref_sequence.length();

  // extract mutations of sequences one by one
  for (std::vector<std::string>::size_type i = 0; i < str_sequences.size(); ++i) {
    // validate the sequence length
    string str_sequence = str_sequences[i];
    if (seq_length != str_sequence.length()) {
      throw std::logic_error(
          "The sequence length of " + seq_names[(std::vector<std::string>::size_type) i] + " (" +
          convertIntToString((int) str_sequence.length()) +
          ") is different from that of the reference sequence (" +
          convertIntToString((int) ref_sequence.length()) + ")!");
    }

    // init new sequence instance for the inference process afterwards
    data.push_back(string(seq_names[(std::vector<std::string>::size_type) i]));
    sequence = &data.back();

    // init dummy variables
    int state = 0;
    PositionType length = 0;
    for (std::basic_string<char>::size_type pos = 0; pos < seq_length; ++pos) {
      switch (state) {
        case 0:  // previous character is neither 'N' nor '-'
          if (str_sequence[pos] != ref_sequence[pos]) {
            length = 1;

            // starting a sequence of 'N'
            if (toupper(str_sequence[pos]) == 'N' &&
                getSeqType() == cmaple::SeqRegion::SEQ_DNA) {
              state = 1;
              // starting a sequence of '-'
            } else if (str_sequence[pos] == '-') {
              state = 2;
              // output a mutation
            } else {
              addMutation(sequence, str_sequence[pos], (PositionType) pos);
            }
          }
          break;
        case 1:  // previous character is 'N'
          // inscrease the length if the current character is still 'N'
          if (toupper(str_sequence[pos]) == 'N' &&
              str_sequence[pos] != ref_sequence[pos]) {
            ++length;
          } else {
            // output the previous sequence of 'N'
            addMutation(sequence, str_sequence[pos - 1], ((PositionType) pos) - length, length);

            // reset state
            state = 0;

            // handle new character different from the reference
            if (str_sequence[pos] != ref_sequence[pos]) {
              length = 1;
              // starting a sequence of '-'
              if (str_sequence[pos] == '-') {
                state = 2;
                // output a mutation
              } else {
                addMutation(sequence, str_sequence[pos], (PositionType) pos);
                state = 0;
              }
            }
          }
          break;
        case 2:  // previous character is '-'
          // inscrease the length if the current character is still '-'
          if (toupper(str_sequence[pos]) == '-' &&
              str_sequence[pos] != ref_sequence[pos]) {
            ++length;
          } else {
            // output the previous sequence of '-'
            addMutation(sequence, str_sequence[pos - 1], ((PositionType) pos) - length, length);

            // reset state
            state = 0;

            // handle new character different from the reference
            if (str_sequence[pos] != ref_sequence[pos]) {
              length = 1;
              // starting a sequence of 'N'
              if (toupper(str_sequence[pos]) == 'N' &&
                  getSeqType() == cmaple::SeqRegion::SEQ_DNA) {
                state = 1;
                // output a mutation
              } else {
                addMutation(sequence, str_sequence[pos], (PositionType) pos);
                state = 0;
              }
            }
          }
          break;
      }
    }

    //  output the last sequence of 'N' or '-' (if any)
    if (state != 0) {
      addMutation(sequence, str_sequence[str_sequence.length() - 1],
                  ((PositionType) str_sequence.length()) - length, length);
    }
  }
}

void cmaple::Alignment::parseRefSeq(string& ref_sequence) {
  assert(ref_sequence.length() > 0);
    
  ref_seq.resize(ref_sequence.length());
  const char DEFAULT_CHAR = cmaple::Alignment::convertState2Char(0, seq_type_);

  // transform ref_sequence to uppercase
  transform(ref_sequence.begin(), ref_sequence.end(), ref_sequence.begin(),
            ::toupper);

  for (std::vector<StateType>::size_type i = 0; i < ref_sequence.length(); ++i) {
    ref_seq[i] = convertChar2State(ref_sequence[i]);

    // validate the ref state
    if (ref_seq[i] >= num_states) {
      ref_sequence[i] = DEFAULT_CHAR;
      if (cmaple::verbose_mode > cmaple::VB_QUIET) {
        outWarning("Invalid reference state found at site " +
                   convertPosTypeToString((PositionType) i) +
                   " was replaced by a default state " + DEFAULT_CHAR);
      }
      ref_seq[i] = 0;
    }
  }
}

void cmaple::Alignment::readMaple(std::istream& in) {
  // init dummy variables
  string seq_name;
  vector<Mutation> mutations;
  PositionType line_num = 1;
  string line;

  // set the failbit and badbit
  in.exceptions(ios::failbit | ios::badbit);
  // remove the failbit
  in.exceptions(ios::badbit);

  if (cmaple::verbose_mode >= cmaple::VB_MAX) {
    cout << "Reading an alignment in MAPLE format from a stream" << endl;
  }

  // extract reference sequence first
  for (; !in.eof(); ++line_num) {
    safeGetline(in, line);
    if (line == "") {
      continue;
    }

    // read the first line (">REF")
    if (line[0] == '>') {
      string::size_type pos = line.find_first_of("\n\r");
      seq_name = line.substr(1, pos - 1);

      // transform seq_name to upper case
      transform(seq_name.begin(), seq_name.end(), seq_name.begin(), ::toupper);

      if (seq_name != REF_NAME && seq_name != "REFERENCE") {
        throw std::logic_error(
            "MAPLE file must start by >REF. Please check and try again!");
      }
    }
    // read the reference sequence
    else {
      // make sure the first line was found
      if (seq_name != REF_NAME && seq_name != "REFERENCE") {
        throw std::logic_error(
            "MAPLE file must start by >REF. Please check and try again!");
      }

      // transform ref_sequence to uppercase
      transform(line.begin(), line.end(), line.begin(), ::toupper);

      // detect the seq_type from the ref_sequences
      const cmaple::SeqRegion::SeqType current_seq_type = getSeqType();
      if (current_seq_type == cmaple::SeqRegion::SEQ_AUTO ||
          current_seq_type == cmaple::SeqRegion::SEQ_UNKNOWN) {
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
  for (; !in.eof(); ++line_num) {
    safeGetline(in, line);
    if (line == "") {
      continue;
    }

    // Read sequence name
    if (line[0] == '>') {
      // record the sequence of the previous taxon
      if (seq_name.length()) {
        data.emplace_back(std::move(seq_name), std::move(mutations));

        // reset dummy variables
        seq_name.clear();
        mutations.clear();
      }

      // Read new sequence name
      string::size_type pos = line.find_first_of("\n\r");
      seq_name = line.substr(1, pos - 1);
      if (!seq_name.length()) {
        throw std::logic_error("Empty sequence name found at line " +
                               convertIntToString(line_num) +
                               ". Please check and try again!");
      }
    }
    // Read a Mutation
    else {
      // validate the input
      char separator = '\t';
      long num_items = std::count(line.begin(), line.end(), separator) + 1;
      if (num_items < 2 || num_items > 3) {
        throw std::logic_error(
            "Invalid input. Each difference must be presented be <Type>    "
            "<Position>  [<Length>]. Please check and try again!");
      }

      // extract mutation info
      stringstream ssin(line);
      string tmp;

      // extract <Type>
      ssin >> tmp;
      StateType state = convertChar2State(toupper(tmp[0]));

      // extract <Position>
      ssin >> tmp;
      PositionType pos = convert_positiontype(tmp.c_str());
      if (pos <= 0 || pos > (PositionType)ref_seq.size()) {
        throw std::logic_error(
            "<Position> must be greater than 0 and less than the reference "
            "sequence length (" +
            convertPosTypeToString((PositionType) ref_seq.size()) + ")!");
      }

      // extract <Length>
      PositionType length = 1;
      if (ssin.good()) {
        ssin >> tmp;
        if (state == TYPE_N || state == TYPE_DEL) {
          length = convert_positiontype(tmp.c_str());
          if (length <= 0) {
            throw std::logic_error("<Length> must be greater than 0!");
          }
          if (length + pos - 1 > (PositionType)ref_seq.size()) {
            throw std::logic_error(
                "<Length> + <Position> must be less than the reference "
                "sequence length (" +
                convertPosTypeToString((PositionType) ref_seq.size()) + ")!");
          }
        } else if (cmaple::verbose_mode >= cmaple::VB_MED) {
          outWarning("Ignoring <Length> of " + tmp +
                     ". <Length> is only appliable for 'N' or '-'.");
        }
      }

      // add a new mutation into mutations
      if (state == TYPE_N || state == TYPE_DEL) {
        mutations.emplace_back(state, pos - 1, length);
      } else {
        mutations.emplace_back(state, pos - 1);
      }
    }
  }

  // Record the sequence of  the last taxon
  if (seq_name.length()) {
    data.emplace_back(std::move(seq_name), std::move(mutations));
  }

  // validate the input
  assert(ref_seq.size() > 0);
  if (ref_seq.size() == 0) {
    throw std::logic_error("Reference sequence is not found!");
  }
  if (data.size() < MIN_NUM_TAXA) {
    throw std::logic_error("The number of taxa must be at least " +
                           convertIntToString(MIN_NUM_TAXA));
  }

  resetStream(in);
}

auto cmaple::Alignment::convertState2Char(
    const cmaple::StateType& state,
    const cmaple::SeqRegion::SeqType& seqtype) -> char {
  if (state == TYPE_N || state == TYPE_DEL) {
    return '-';
  }
  if (state > TYPE_INVALID) {
    return '?';
  }

  switch (seqtype) {
    /*case cmaple::SeqRegion::SEQ_BINARY:
        switch (state) {
        case 0:
            return '0';
        case 1:
            return '1';
        default:
            return '?'; // unrecognize character
        }*/
    case cmaple::SeqRegion::SEQ_DNA:  // DNA
      switch (state) {
        case 0:
          return 'A';
        case 1:
          return 'C';
        case 2:
          return 'G';
        case 3:
          return 'T';
        case 1 + 4 + 3:
          return 'R';  // A or G, Purine
        case 2 + 8 + 3:
          return 'Y';  // C or T, Pyrimidine
        case 1 + 8 + 3:
          return 'W';  // A or T, Weak
        case 2 + 4 + 3:
          return 'S';  // G or C, Strong
        case 1 + 2 + 3:
          return 'M';  // A or C, Amino
        case 4 + 8 + 3:
          return 'K';  // G or T, Keto
        case 2 + 4 + 8 + 3:
          return 'B';  // C or G or T
        case 1 + 2 + 8 + 3:
          return 'H';  // A or C or T
        case 1 + 4 + 8 + 3:
          return 'D';  // A or G or T
        case 1 + 2 + 4 + 3:
          return 'V';  // A or G or C
        default:
          return '?';  // unrecognize character
      }
      return state;
    case cmaple::SeqRegion::SEQ_PROTEIN:  // Protein
      if (state < 20) {
        return symbols_protein[(StateType)state];
      } else if (state == 20) {
        return 'B';
      } else if (state == 21) {
        return 'Z';
      } else if (state == 22) {
        return 'J';
        //        else if (state == 4+8+19) return 'B';
        //        else if (state == 32+64+19) return 'Z';
      } else {
        return '-';
      }
    /*case SEQ_MORPH:
        // morphological state
        if (state < strlen(symbols_morph))
            return symbols_morph[state];
        else
            return '-';*/
    default:
      // unknown
      return '*';
  }
}

auto cmaple::Alignment::convertChar2State(char state) -> StateType {
  if (state == '-') {
    return TYPE_DEL;
  }
  if (state == '?' || state == '.' || state == '~') {
    return TYPE_N;
  }

  char* loc;

  switch (getSeqType()) {
    /*case SEQ_BINARY:
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
                throw std::invalid_argument(invalid_state_msg);
                return TYPE_INVALID;
            }
        }
        break;*/
    case cmaple::SeqRegion::SEQ_DNA:  // DNA
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
          return 1 + 4 + 3;  // A or G, Purine
        case 'Y':
          return 2 + 8 + 3;  // C or T, Pyrimidine
        case 'O':
        case 'N':
        case 'X':
          return TYPE_N;
        case 'W':
          return 1 + 8 + 3;  // A or T, Weak
        case 'S':
          return 2 + 4 + 3;  // G or C, Strong
        case 'M':
          return 1 + 2 + 3;  // A or C, Amino
        case 'K':
          return 4 + 8 + 3;  // G or T, Keto
        case 'B':
          return 2 + 4 + 8 + 3;  // C or G or T
        case 'H':
          return 1 + 2 + 8 + 3;  // A or C or T
        case 'D':
          return 1 + 4 + 8 + 3;  // A or G or T
        case 'V':
          return 1 + 2 + 4 + 3;  // A or G or C
        default: {
          string invalid_state_msg = "Invalid state ";
          invalid_state_msg += state;
          invalid_state_msg += ". Please check and try again!";
          throw std::invalid_argument(invalid_state_msg);
          return TYPE_INVALID;  // unrecognize character
        }
      }
      return (StateType) state;
    case cmaple::SeqRegion::SEQ_PROTEIN:  // Protein
      //        if (state == 'B') return 4+8+19;
      //        if (state == 'Z') return 32+64+19;
      if (state == 'B') {
        return 20;
      }
      if (state == 'Z') {
        return 21;
      }
      if (state == 'J') {
        return 22;
      }
      if (state == '*') {
        return TYPE_N;  // stop codon
      }
      if (state == 'U') {
        return TYPE_N;  // 21st amino-acid
      }
      if (state == 'O') {
        return TYPE_N;  // 22nd amino-acid
      }
      loc = strchr(symbols_protein, state);

      if (!loc) {
        string invalid_state_msg = "Invalid state ";
        invalid_state_msg += state;
        invalid_state_msg += ". Please check and try again!";
        throw std::invalid_argument(invalid_state_msg);
        return TYPE_INVALID;  // unrecognize character
      }
      state = loc - symbols_protein;
      if (state < 20) {
        return (StateType) state;
      } else {
        return TYPE_N;
      }
    /*case SEQ_MORPH: // Standard morphological character
        loc = strchr(symbols_morph, state);

        if (!loc)
        {
            string invalid_state_msg = "Invalid state ";
            invalid_state_msg += state;
            invalid_state_msg += ". Please check and try again!";
            throw std::invalid_argument(invalid_state_msg);
            return TYPE_INVALID; // unrecognize character
        }
        state = loc - symbols_morph;
        return state;*/
    default: {
      string invalid_state_msg = "Invalid state ";
      invalid_state_msg += state;
      invalid_state_msg += ". Please check and try again!";
      throw std::invalid_argument(invalid_state_msg);
      return TYPE_INVALID;
    }
  }
}

auto cmaple::Alignment::computeSeqDistance(Sequence& sequence,
                                           RealNumType hamming_weight)
    -> PositionType {        
  // dummy variables
  PositionType num_ambiguities = 0;
  PositionType num_diffs = 0;

  // browse mutations one by one
  for (const auto& mutation : sequence) {
    switch (mutation.type) {
      case TYPE_R:  // Type R does not exist in Mutations (only in Regions)
        throw std::logic_error(
            "Sorry! Something went wrong. Invalid mutation type (type R).");
        break;
      case TYPE_N:  // handle unsequenced sites ('-' or 'N')
      case TYPE_DEL:
      case TYPE_O:  // handle ambiguity
        num_ambiguities += mutation.getLength();
        ++num_diffs;
        break;
      default:
        // handle normal character
        if (mutation.type < num_states) {
          ++num_diffs;
        } else {
          ++num_diffs;
          ++num_ambiguities;
        }
        break;
    }
  }

  // calculate and record the distance of the current sequence
  return ((PositionType) num_diffs * hamming_weight) + num_ambiguities;
}

void cmaple::Alignment::sortSeqsByDistances() {
  // init dummy variables
  const RealNumType hamming_weight = 1000;
  std::vector<cmaple::Sequence>::size_type num_seqs = data.size();
  auto* distances = new PositionType[num_seqs];
  auto* sequence_indexes = new PositionType[num_seqs];

  // calculate the distances of each sequence
  for (PositionType i = 0; i < num_seqs; ++i) {
    // dummy variables
    sequence_indexes[i] = i;

    // calculate and record the distance of the current sequence
    distances[i] = computeSeqDistance(data[(std::vector<cmaple::Sequence>::size_type) i], hamming_weight);

    // NHANLT: debug
    // distances[i] *= 1000;
  }

  // NHANLT: debug

  // sort distances
  quicksort(distances, 0, num_seqs - 1, sequence_indexes);

  // re-order sequences by distances
  vector<Sequence> tmp_sequences(std::move(data));
  data.reserve(num_seqs);
  for (PositionType i = 0; i < num_seqs; ++i) {
    data.push_back(std::move(tmp_sequences[(std::vector<cmaple::Sequence>::size_type) sequence_indexes[i]]));
  }

  // delete distances, sequence_indexes
  delete[] distances;
  delete[] sequence_indexes;
}

std::string cmaple::Alignment::getRefSeqStr() {
  const std::basic_string<char>::size_type seq_length = ref_seq.size();
  std::string ref_sequence(seq_length, ' ');
  for (std::vector<cmaple::StateType>::size_type i = 0; i < seq_length; ++i) {
    ref_sequence[i] =
        cmaple::Alignment::convertState2Char(ref_seq[i], seq_type_);
  }
  
  assert(ref_sequence.length() == seq_length);
  return ref_sequence;
}

std::string cmaple::Alignment::getSeqString(const std::string& ref_seq_str,
                                            Sequence* sequence) {
  // clone the sequence
  std::string sequence_str = ref_seq_str;
  // apply mutations in sequence_str
  Mutation* mutation = &sequence->front();
  for (std::vector<cmaple::Mutation>::size_type j = 0; j < sequence->size(); ++j, ++mutation) {
    char state =
        cmaple::Alignment::convertState2Char(mutation->type, seq_type_);

    // replace characters in sequence_str
    for (std::basic_string<char>::size_type pos = mutation->position;
         pos < mutation->position + mutation->getLength(); ++pos) {
      sequence_str[pos] = state;
    }
  }

  return sequence_str;
}

void cmaple::Alignment::writeMAPLE(std::ostream& out) {
  assert(data.size() > 0);
    
  // write reference sequence to the output file
  out << ">" << REF_NAME << endl;
  out << getRefSeqStr() << endl;

  // Write sequences one by one
  const std::vector<cmaple::Sequence>::size_type num_seqs = data.size();
  Sequence* sequence = &data.front();
  for (std::vector<cmaple::Sequence>::size_type i = 0; i < num_seqs; ++i, ++sequence) {
    // write the sequence name
    out << ">" << sequence->seq_name << endl;

    // write mutations
    Mutation* mutation = &sequence->front();
    for (std::vector<cmaple::Mutation>::size_type j = 0; j < sequence->size(); ++j, ++mutation) {
      const StateType type = mutation->type;
      out << cmaple::Alignment::convertState2Char(type, seq_type_) << "\t"
          << (mutation->position + 1);
      if (type == TYPE_N || type == TYPE_DEL) {
        out << "\t" << mutation->getLength();
      }
      out << endl;
    }
  }
}

void cmaple::Alignment::writeFASTA(std::ostream& out) {
  assert(data.size() > 0);
    
  // Get reference sequence
  const std::string ref_sequence = getRefSeqStr();

  // Write sequences one by one
  const std::vector<cmaple::Sequence>::size_type num_seqs = data.size();
  Sequence* sequence = &data.front();
  for (std::vector<cmaple::Sequence>::size_type i = 0; i < num_seqs; ++i, ++sequence) {
    // write the sequence name
    out << ">" << sequence->seq_name << endl;

    // write the sequence
    out << getSeqString(ref_sequence, sequence) << std::endl;
  }
}

void cmaple::Alignment::writePHYLIP(std::ostream& out) {
  assert(data.size() > 0);
    
  // Write the header <num_seqs> <seq_length>
  const std::vector<cmaple::StateType>::size_type seq_length = ref_seq.size();
  const std::vector<cmaple::Sequence>::size_type num_seqs = data.size();
  out << num_seqs << "\t" << seq_length << std::endl;

  // Get reference sequence
  const std::string ref_sequence = getRefSeqStr();

  // Get max length of sequence names
  std::basic_string<char>::size_type max_name_length = 9;
  Sequence* sequence = &data.front();
  for (std::vector<cmaple::Sequence>::size_type i = 0; i < num_seqs; ++i, ++sequence) {
    if (sequence->seq_name.length() > max_name_length) {
      max_name_length = sequence->seq_name.length();
    }
  }

  // Add one extra space
  ++max_name_length;

  // Write sequences one by one
  sequence = &data.front();
  for (std::vector<cmaple::Sequence>::size_type i = 0; i < num_seqs; ++i, ++sequence) {
    // write the sequence name
    std::string seq_name = sequence->seq_name;
    seq_name.resize(max_name_length, ' ');
    out << seq_name;

    // write the sequence
    out << getSeqString(ref_sequence, sequence) << std::endl;
  }
}

void cmaple::Alignment::readFastaOrPhylip(std::istream& aln_stream,
                                          const std::string& ref_seq) {
  assert(aln_format != IN_UNKNOWN);
    
  // read input sequences
  if (aln_format == IN_UNKNOWN) {
    throw std::logic_error("Unknown alignment format");
  }
  StrVector sequences;
  StrVector seq_names;
  readSequences(aln_stream, sequences, seq_names, aln_format);

  // validate the input sequences
  assert(sequences.size() > 0);
  if (sequences.size() == 0) {
    throw std::logic_error(
        "Empty input sequences. Please check and try again!");
  }
  // make sure all sequences have the same length
  if (aln_format == IN_FASTA) {
    for (std::vector<cmaple::Sequence>::size_type i = 0; i < sequences.size(); ++i) {
      if (sequences[i].length() != sequences[0].length()) {
        throw std::logic_error(
            "Sequence " + seq_names[(std::vector<std::string>::size_type) i] +
            " has a different length compared to the first sequence.");
      }
    }
  }

  // detect the type of the input sequences
  const cmaple::SeqRegion::SeqType current_seq_type = getSeqType();
  if (current_seq_type == cmaple::SeqRegion::SEQ_AUTO ||
      current_seq_type == cmaple::SeqRegion::SEQ_UNKNOWN) {
    setSeqType(detectSequenceType(sequences));
  }

  // generate reference sequence from the input sequences
  string ref_sequence;
  // read the reference sequence from file (if the user supplies it)
  if (ref_seq.length()) {
    ref_sequence = ref_seq;
  } else {
    ref_sequence = generateRef(sequences);
  }

  // parse ref_sequence into vector of states
  parseRefSeq(ref_sequence);
    
  assert(ref_sequence.length() > 0);

  // extract and write mutations of each sequence to file
  extractMutations(sequences, seq_names, ref_sequence);
}

auto cmaple::Alignment::detectSequenceType(StrVector& sequences)
    -> cmaple::SeqRegion::SeqType {
  size_t num_nuc = 0;
  size_t num_ungap = 0;
  size_t num_bin = 0;
  size_t num_alpha = 0;
  size_t num_digit = 0;
  double detectStart = getRealTime();
  size_t sequenceCount = sequences.size();
  assert(sequenceCount > 0);

  // #pragma omp parallel for
  // reduction(+:num_nuc,num_ungap,num_bin,num_alpha,num_digit)
  for (size_t seqNum = 0; seqNum < sequenceCount; ++seqNum) {
    auto start = sequences.at(seqNum).data();
    auto stop = start + sequences.at(seqNum).size();
    for (auto i = start; i != stop; ++i) {
      if ((*i) == 'A' || (*i) == 'C' || (*i) == 'G' || (*i) == 'T' ||
          (*i) == 'U') {
        ++num_nuc;
        ++num_ungap;
        continue;
      }
      if ((*i) == '?' || (*i) == '-' || (*i) == '.') {
        continue;
      }
      if (*i != 'N' && *i != 'X' && (*i) != '~') {
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
  }  // end OMP parallel for

  if (verbose_mode >= VB_DEBUG) {
    cout << "Sequence Type detection took " << (getRealTime() - detectStart)
         << " seconds." << endl;
  }
  if (((double)num_nuc) / num_ungap > 0.9) {
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG) {
      std::cout << "DNA data detected." << std::endl;
    }
    return cmaple::SeqRegion::SEQ_DNA;
  }
  /*if (((double)num_bin) / num_ungap > 0.9)
  {
      if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
          std::cout << "Binary data detected." << std::endl;
      return SEQ_BINARY;
  }*/
  if (((double)num_alpha + num_nuc) / num_ungap > 0.9) {
    if (cmaple::verbose_mode >= cmaple::VB_DEBUG) {
      std::cout << "Protein data detected." << std::endl;
    }
    return cmaple::SeqRegion::SEQ_PROTEIN;
  }
  /*if (((double)(num_alpha + num_digit + num_nuc)) / num_ungap > 0.9)
  {
      if (cmaple::verbose_mode >= cmaple::VB_DEBUG)
          std::cout << "Morphological data detected." << std::endl;
      return SEQ_MORPH;
  }*/
  return cmaple::SeqRegion::SEQ_UNKNOWN;
}

void cmaple::Alignment::updateNumStates() {
  switch (getSeqType()) {
    case cmaple::SeqRegion::SEQ_PROTEIN:
      num_states = 20;
      break;

    default:  // dna
      num_states = 4;
      break;
  }
}

auto cmaple::Alignment::detectMAPLEorFASTA(std::istream& in)
    -> cmaple::Alignment::InputType {
  PositionType num_taxa = 0;
  string line;

  // set the failbit and badbit
  in.exceptions(ios::failbit | ios::badbit);
  // remove the failbit
  in.exceptions(ios::badbit);

  // extract reference sequence first
  for (; !in.eof();) {
    safeGetline(in, line);
    if (line == "") {
      continue;
    }

    // count lines that begin with >
    if (line[0] == '>') {
      ++num_taxa;
      // handle other lines that presents a sequence (in FASTA) or a mutation
      // (in MAPLE)
    } else {
      // only check the line if we pass the first two lines that begins with >
      if (num_taxa >= 2) {
        // mutation must contain either a tab or a space in the middle
        auto pos_tab = line.find("\t");
        auto pos_space = line.find(" ");
        auto last_index = line.length() - 1;

        if ((pos_tab != string::npos && pos_tab < last_index) ||
            (pos_space != string::npos && pos_space < last_index)) {
          // reset aln_stream
          resetStream(in);
          return cmaple::Alignment::IN_MAPLE;
        } else {
          // reset aln_stream
          resetStream(in);

          return cmaple::Alignment::IN_FASTA;
        }
      }
    }
  }

  // reset aln_stream
  resetStream(in);

  return cmaple::Alignment::IN_FASTA;
}

auto cmaple::Alignment::detectInputFile(std::istream& in)
    -> cmaple::Alignment::InputType {
  unsigned char ch = ' ';
  unsigned char ch2 = ' ';
  int count = 0;
  do {
    in >> ch;
  } while (ch <= 32 && !in.eof() && count++ < 20);
  in >> ch2;
  // reset aln_stream
  resetStream(in);
  switch (ch) {
    // case '#': return IN_NEXUS;
    // case '(': return IN_NEWICK;
    // case '[': return IN_NEWICK;
    case '>':
      return detectMAPLEorFASTA(in);
    /*case 'C': if (ch2 == 'L') return IN_CLUSTAL;
              else if (ch2 == 'O') return IN_COUNTS;
              else return IN_OTHER;
    case '!': if (ch2 == '!') return IN_MSF; else return IN_OTHER;*/
    default:
      if (isdigit(ch)) {
        return cmaple::Alignment::IN_PHYLIP;
      }
      return cmaple::Alignment::IN_UNKNOWN;
  }

  return cmaple::Alignment::IN_UNKNOWN;
}

auto cmaple::Alignment::parseAlnFormat(const std::string& n_format)
    -> cmaple::Alignment::InputType {
  // transform to uppercase
  string format(n_format);
  transform(format.begin(), format.end(), format.begin(), ::toupper);
  if (format == "MAPLE") {
    return cmaple::Alignment::IN_MAPLE;
  }
  if (format == "PHYLIP") {
    return cmaple::Alignment::IN_PHYLIP;
  }
  if (format == "FASTA") {
    return cmaple::Alignment::IN_FASTA;
  }
  if (format == "AUTO") {
    return cmaple::Alignment::IN_AUTO;
  }

  // default
  return cmaple::Alignment::IN_UNKNOWN;
}
