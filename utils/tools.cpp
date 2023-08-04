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



#include "tools.h"
#include "timeutil.h"

//#include <filesystem>

using namespace std;
using namespace cmaple;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

VerboseMode cmaple::verbose_mode = VB_MED;

void cmaple::printCopyright(ostream &out) {
    out << "CMAPLE version ";
    out << cmaple_VERSION_MAJOR << "." << cmaple_VERSION_MINOR << cmaple_VERSION_PATCH;
    out << " for " << getOSName();
    out << " built " << __DATE__ << std::endl;
    
    out << "Developed by Nhan Ly-Trong, Chris Bielow, Nicola De Maio, Bui Quang Minh." << std::endl << std::endl;
}

/**
        Output an error to screen, then exit program
        @param error error message
 */
void cmaple::outError(const char *error, bool quit) {
	if (error == ERR_NO_MEMORY) {
        //print_stacktrace(cerr);
	}
	cerr << error << endl;
    if (quit)
    	exit(2);
}

/**
        Output an error to screen, then exit program
        @param error error message
 */
void cmaple::outError(const string &error, bool quit) {
    outError(error.c_str(), quit);
}

void cmaple::outError(const char *error, const char *msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

void cmaple::outError(const char *error, const string &msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

/**
        Output a warning message to screen
        @param error warning message
 */
void cmaple::outWarning(const char *warn) {
    cout << "WARNING: " << warn << endl;
}

void cmaple::outWarning(const string &warn) {
    outWarning(warn.c_str());
}

std::istream& cmaple::safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

//From Tung
string cmaple::convertPosTypeToString(PositionType number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string cmaple::convertIntToString(int number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string cmaple::convertInt64ToString(int64_t number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string cmaple::convertDoubleToString(RealNumType number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

std::string cmaple::convertDoubleToString(RealNumType number, uint8_t precision)
{
    stringstream ss; //create a stringstream
    ss << std::setprecision(precision) << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

bool cmaple::iEquals(const string &a, const string &b)
{
    unsigned int sz = a.size();
    if (b.size() != sz)
        return false;
    for (unsigned int i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}

//From Tung

bool cmaple::copyFile(const char SRC[], const char DEST[]) {
    std::ifstream src; // the source file
    std::ofstream dest; // the destination file

    src.open(SRC, std::ios::binary); // open in binary to prevent jargon at the end of the buffer
    dest.open(DEST, std::ios::binary); // same again, binary
    if (!src.is_open() || !dest.is_open())
        return false; // could not be copied

    dest << src.rdbuf(); // copy the content
    dest.close(); // close destination file
    src.close(); // close source file

    return true; // file copied successfully
}

bool cmaple::fileExists(const string &strFilename) {
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0) {
        // We were able to get the file attributes
        // so the file obviously exists.
        blnReturn = true;
    } else {
        // We were not able to get the file attributes.
        // This may mean that we don't have permission to
        // access the folder which contains this file. If you
        // need to do that level of checking, lookup the
        // return values of stat which will give you
        // more details on why stat failed.
        blnReturn = false;
    }
    return (blnReturn);
}

/*int isDirectory(const char *path) {
  return std::filesystem::is_directory(path);
}

int isFile(const char *path) {
  return std::filesystem::is_regular_file(path);
}*/

int cmaple::convert_int(const char *str, int &end_pos) {
	char *endptr;
	int i = strtol(str, &endptr, 10);

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw std::invalid_argument(err);
	}
	end_pos = endptr - str;
	return i;
}

int cmaple::convert_int(const char *str) {
    char *endptr;
    int i = strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }

    return i;
}

PositionType cmaple::convert_positiontype(const char *str) {
    char *endptr;
    PositionType i = (PositionType) strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || i == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }

    return i;
}

void cmaple::convert_int_vec(const char *str, IntVector &vec) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		int i = strtol(beginptr, &endptr, 10);

		if ((i == 0 && endptr == beginptr) || abs(i) == HUGE_VALL) {
			string err = "Expecting integer, but found \"";
			err += beginptr;
			err += "\" instead";
            throw std::invalid_argument(err);
		}
		vec.push_back(i);
		if (*endptr == ',') endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}


int64_t cmaple::convert_int64(const char *str) {
    char *endptr;
    int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting large integer , but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }

    return i;
}

int64_t cmaple::convert_int64(const char *str, int &end_pos) {
	char *endptr;
	int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting large integer, but found \"";
		err += str;
		err += "\" instead";
        throw std::invalid_argument(err);
	}
	end_pos = endptr - str;
	return i;
}


RealNumType cmaple::convert_real_number(const char *str) {
    char *endptr;
    RealNumType d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }
    return d;
}

RealNumType cmaple::convert_real_number(const char *str, int &end_pos) {
	char *endptr;
	RealNumType d = strtod(str, &endptr);
	if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
        throw std::invalid_argument(err);
	}
	end_pos = endptr - str;
	return d;
}

void cmaple::convert_real_numbers(RealNumType* &arr, string input_str)
{
    // count the number of input real_numbers
    int number_count = count(input_str.begin(), input_str.end(), ' ') + 1;
    
    // init array
    arr = new RealNumType[number_count];
    
    // parse rates
    stringstream ss(input_str);
    int index = 0;
    while (ss.good())
    {
        ss >> arr[index];
        ++index;
    }
}

void cmaple::convert_real_number_vec(const char *str, RealNumberVector &vec, char separator) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		RealNumType d = strtod(beginptr, &endptr);

		if ((d == 0.0 && endptr == beginptr) || fabs(d) == HUGE_VALF) {
			string err = "Expecting floating-point number, but found \"";
			err += beginptr;
			err += "\" instead";
            throw std::invalid_argument(err);
		}
		vec.push_back(d);
		if (*endptr == separator) endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}

string cmaple::convert_time(const RealNumType sec) {
    int sec_int = (int) floor(sec);
    int secs = sec_int % 60;
    int mins = (sec_int % 3600) / 60;
    int hours = sec_int / 3600;
    stringstream ss;
    ss << hours << "h:" << mins << "m:" << secs << "s";
    return ss.str();
}

void cmaple::convert_range(const char *str, int &lower, int &upper, int &step_size) {
    char *endptr;

    // parse the lower bound of the range
    int d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }
    //lower = d;
    int d_save = d;
    upper = d;
    if (*endptr == 0) return;


    // parse the upper bound of the range
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }
    step_size = d;
}

void cmaple::convert_range(const char *str, RealNumType &lower, RealNumType &upper, RealNumType &step_size) {
    char *endptr;

    // parse the lower bound of the range
    RealNumType d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }
    //lower = d;
    RealNumType d_save = d;
    upper = d;
    if (*endptr == 0) return;


    // parse the upper bound of the range
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw std::invalid_argument(err);
    }
    step_size = d;
}

void cmaple::reinitDoubleArr(RealNumType* &arr, StateType size, bool delete_first, bool set_zero)
{
    // delete the current array
    if (delete_first && arr)
        delete [] arr;
    
    // request memory allocation for the new array
    arr = new RealNumType[size];
    if (set_zero)
        for (StateType i = 0; i < size; ++i)
            arr[i] = 0;
}

void cmaple::convert_string_vec(const char *str, StrVector &vec, char separator) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    string elem;
    do {
    	endptr = strchr(beginptr, separator);
    	if (!endptr) {
    		elem.assign(beginptr);
    		vec.push_back(elem);
    		return;
    	}
    	elem.assign(beginptr, endptr-beginptr);
    	vec.push_back(elem);
		beginptr = endptr+1;
    } while (*endptr != 0);

}

void cmaple::normalize_frequencies_from_index(RealNumType* freqs, int num_states, int starting_index)
{
    ASSERT(num_states > 0);
    // calculate the total_freqs
    RealNumType total_freqs = 0;
    for (int i = starting_index; i < starting_index+num_states; ++i)
        total_freqs += freqs[i];
    
    // normalize the freqs
    if (fabs(total_freqs) < 1e-5)
        throw std::logic_error("Sum of state frequencies must be greater than zero!");
    total_freqs = 1.0 / total_freqs;
    for (int i = starting_index; i < starting_index+num_states; ++i)
        freqs[i] *= total_freqs;
}

bool cmaple::is_number(const std::string& s)
{
    char* end = nullptr;
    double val = strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && val != HUGE_VAL;
}

cmaple::Params::Params()
{
    initDefaultValue();
}

void cmaple::Params::initDefaultValue()
{
    aln_path = "";
    aln_format_str = "AUTO";
    ref_path = "";
    ref_seqname = "";
    sub_model_str = "GTR";
    fixed_blengths = false;
    overwrite_output = false;
    threshold_prob = 1e-8;
    mutation_update_period = 25;
    failure_limit_sample = 5;
    failure_limit_subtree = 4;
    failure_limit_subtree_short_search = 1;
    strict_stop_seeking_placement_sample = false;
    strict_stop_seeking_placement_subtree = false;
    strict_stop_seeking_placement_subtree_short_search = true;
    thresh_log_lh_sample = 200;
    thresh_log_lh_subtree = 160;
    thresh_log_lh_subtree_short_search = 40;
    thresh_log_lh_failure = 0.01;
    min_blength_factor = 0.2;
    min_blength_mid_factor = 4.1;
    max_blength_factor = 40;
    thresh_diff_update = 1e-7;
    thresh_diff_fold_update = 1.001;
    output_aln = "";
    output_aln_format_str = "MAPLE";
    num_tree_improvement = 1;
    thresh_entire_tree_improvement = 1;
    thresh_placement_cost = -1e-5;
    thresh_placement_cost_short_search = -1;
    tree_format_str = "BIN";
    shallow_tree_search = false;
    output_testing = NULL;
    compute_aLRT_SH = false;
    aLRT_SH_replicates = 1000;
    aLRT_SH_half_epsilon = 0.05;
    num_threads = 1;
    input_treefile = "";
    output_prefix = "";
    allow_replace_input_tree = false;
    fixed_min_blength = -1;
    seq_type_str = "AUTO";
    tree_search_type_str = "NORMAL";
    make_consistent = false;
    
    // initialize random seed based on current time
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    //ran_seed = (unsigned) (tv.tv_sec+tv.tv_usec);
    ran_seed = (tv.tv_usec);
}

void cmaple::parseArg(int argc, char *argv[], Params &params) {
    for (int cnt = 1; cnt < argc; ++cnt) {
        try {
            if (strcmp(argv[cnt], "--alignment") == 0 || strcmp(argv[cnt], "-aln") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -aln <ALN_FILENAME>");
                
                params.aln_path = argv[cnt];

                continue;
            }
            // "--diff" should be removed, I temporarily keep it for compatibility purpose
            if (strcmp(argv[cnt], "--diff") == 0 || strcmp(argv[cnt], "-diff") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -diff <ALN_FILENAME>");
                
                params.aln_path = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--prefix") == 0 || strcmp(argv[cnt], "-prf") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -prf <OUTPUT_PREFIX>");
                
                params.output_prefix = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--output-aln") == 0 || strcmp(argv[cnt], "-out-aln") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -out-aln <ALN_FILENAME>,<ALN_FORMAT>. Note <ALN_FORMAT> could be MAPLE, PHYLIP, FASTA, or AUTO");
                
                // parse inputs
                std::string inputs = argv[cnt];
                std::string delimiter = ",";
                size_t pos = inputs.find(delimiter);
                if (pos != std::string::npos)
                {
                    params.output_aln = inputs.substr(0, pos);
                    // validate output_aln
                    if (!params.output_aln.length())
                        outError("<ALN_FILENAME> is empty!");
                    inputs.erase(0, pos + delimiter.length());
                    params.output_aln_format_str = inputs;
                }
                else
                    outError("Use -out-aln <ALN_FILENAME>,<ALN_FORMAT>. Note <ALN_FORMAT> could be MAPLE, PHYLIP, FASTA, or AUTO");

                continue;
            }
            if (strcmp(argv[cnt], "-st") == 0 || strcmp(argv[cnt], "--seqtype") == 0) {
                cnt++;
                if (cnt >= argc)
                    outError("Use -st DNA or -st AA or -st AUTO");
                params.seq_type_str = argv[cnt];
                
                continue;
            }
            if (strcmp(argv[cnt], "-verbose") == 0 || strcmp(argv[cnt], "--verbose-mode") == 0) {
                cnt++;
                if (cnt >= argc)
                    outError("Use -verbose QUIET/MIN/MED/MAX/DEBUG");
                // parse verbose
                std::string verbose = argv[cnt];
                transform(verbose.begin(), verbose.end(), verbose.begin(), ::toupper);
                if (verbose == "QUIET")
                    verbose_mode = VB_QUIET;
                else if (verbose == "MIN")
                    verbose_mode = VB_MIN;
                else if (verbose == "MED")
                    verbose_mode = VB_MED;
                else if (verbose == "MAX")
                    verbose_mode = VB_MAX;
                else if (verbose == "DEBUG")
                    verbose_mode = VB_DEBUG;
                else
                    outError("Use -verbose QUIET/MIN/MED/MAX/DEBUG");
                continue;
            }
            if (strcmp(argv[cnt], "-format") == 0 || strcmp(argv[cnt], "--aln-format") == 0) {
                cnt++;
                if (cnt >= argc)
                    outError("Use -format MAPLE, PHYLIP, FASTA, or AUTO");
                params.aln_format_str = argv[cnt];
                
                continue;
            }
            if (strcmp(argv[cnt], "--tree") == 0 || strcmp(argv[cnt], "-t") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -t <INPUT_TREEFILE>");
                
                params.input_treefile = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--reference") == 0 || strcmp(argv[cnt], "-ref") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -ref <REF_FILENAME>,<REF_SEQNAME>");
                
                // parse inputs
                std::string inputs = argv[cnt];
                std::string delimiter = ",";
                size_t pos = inputs.find(delimiter);
                if (pos != std::string::npos)
                {
                    params.ref_path = inputs.substr(0, pos);
                    // validate ref_path
                    if (!params.ref_path.length())
                        outError("<REF_FILENAME> is empty!");
                    inputs.erase(0, pos + delimiter.length());
                    params.ref_seqname = inputs;
                    // validate ref_seqname
                    if (!params.ref_seqname.length())
                        outError("<REF_SEQNAME> is empty!");
                }
                else
                    outError("Use -ref <REF_FILENAME>,<REF_SEQNAME>");
            
                continue;
            }
            if (strcmp(argv[cnt], "--min-blength") == 0 || strcmp(argv[cnt], "-min-bl") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -min-bl <NUMBER>");
                try
                {
                    params.fixed_min_blength = convert_real_number(argv[cnt]);
                } catch (std::exception e)
                {
                    outError(e.what());
                }
                
                if (params.fixed_min_blength <= 0)
                    outError("<NUMBER> following -min-bl must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--model") == 0 || strcmp(argv[cnt], "-m") == 0) {
                ++cnt;
                if (cnt >= argc)
                    outError("Use -m <model_name>");
                
                params.sub_model_str = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "--tree-search") == 0 || strcmp(argv[cnt], "-tree-search") == 0) {
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -tree-search <FAST|NORMAL|MORE_ACCURATE>");
                
                params.tree_search_type_str = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "-blfix") == 0 || strcmp(argv[cnt], "-fixbr") == 0 || strcmp(argv[cnt], "--fixed-blength") == 0) {
                params.fixed_blengths = true;
                continue;
            }
            if (strcmp(argv[cnt], "--overwrite") == 0 || strcmp(argv[cnt], "-overwrite") == 0) {
                params.overwrite_output = true;
                continue;
            }
            if (strcmp(argv[cnt], "--make-consistent") == 0 || strcmp(argv[cnt], "-consistent") == 0) {
                params.make_consistent = true;
                continue;
            }
            if (strcmp(argv[cnt], "--replace-input-tree") == 0 || strcmp(argv[cnt], "-rep-tree") == 0) {
                params.allow_replace_input_tree = true;
                continue;
            }
            if (strcmp(argv[cnt], "--threshold-prob") == 0 || strcmp(argv[cnt], "-thresh-prob") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -thresh-prob <PROB_THRESH>");
                
                try
                {
                    params.threshold_prob = convert_real_number(argv[cnt]);
                } catch (std::exception e)
                {
                    outError(e.what());
                }
                
                if (params.threshold_prob <= 0)
                    outError("<PROB_THRESH> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--mutation-update") == 0 || strcmp(argv[cnt], "-mut-update") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -mut-update <NUMBER>");
                
                try
                {
                    params.mutation_update_period = convert_int(argv[cnt]);
                }
                catch(std::exception e)
                {
                    outError(e.what());
                }
                
                if (params.mutation_update_period <= 0)
                    outError("<NUMBER> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--failure-limit") == 0 || strcmp(argv[cnt], "-fail-limit") == 0) {
                
                ++cnt;
                
                try
                {
                    params.failure_limit_sample = convert_int(argv[cnt]);
                }
                catch(std::exception e)
                {
                    outError(e.what());
                }
                
                if (params.failure_limit_sample <= 0)
                    outError("<NUMBER> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--failure-limit-subtree") == 0 || strcmp(argv[cnt], "-fail-limit-stree") == 0) {
                
                ++cnt;
                
                try
                {
                    params.failure_limit_subtree = convert_int(argv[cnt]);
                }
                catch(std::exception e)
                {
                    outError(e.what());
                }
                
                if (params.failure_limit_subtree <= 0)
                    outError("<NUMBER> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--strict-stop-init") == 0 || strcmp(argv[cnt], "-strict-stop-init") == 0) {
                
                params.strict_stop_seeking_placement_sample = true;

                continue;
            }
            if (strcmp(argv[cnt], "--unstrict-stop-subtree") == 0 || strcmp(argv[cnt], "-unstrict-stop-stree") == 0) {
                
                params.strict_stop_seeking_placement_subtree = false;

                continue;
            }
            if (strcmp(argv[cnt], "--output-multifurcating-tree") == 0 || strcmp(argv[cnt], "-out-mul-tree") == 0) {
                
                params.tree_format_str = "MUL";

                continue;
            }
            if (strcmp(argv[cnt], "--shallow-tree-search") == 0 || strcmp(argv[cnt], "-shallow-search") == 0) {
                
                params.shallow_tree_search = true;

                continue;
            }
            if (strcmp(argv[cnt], "--output-testing") == 0 || strcmp(argv[cnt], "-out-test") == 0) {
                
                ++cnt;
                
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -out-test <FILE_PATH>");
                
                params.output_testing = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--branch-support") == 0 || strcmp(argv[cnt], "-br-supp") == 0) {
                
                params.compute_aLRT_SH = true;

                continue;
            }
            if (strcmp(argv[cnt], "--replicates") == 0 || strcmp(argv[cnt], "-rep") == 0) {
                ++cnt;
                if (cnt >= argc)
                    outError("Use --replicates <NUM_REPLICATES>");
                
                try
                {
                    params.aLRT_SH_replicates = convert_int(argv[cnt]);
                }
                catch(std::exception e)
                {
                    outError(e.what());
                }
                
                if (params.aLRT_SH_replicates <= 0)
                    outError("<NUM_REPLICATES> must be positive!");
                continue;
            }
            if (strcmp(argv[cnt], "--epsilon") == 0 || strcmp(argv[cnt], "-eps") == 0) {
                ++cnt;
                if (cnt >= argc)
                    outError("Use --epsilon <FLOATING_NUM>");
                try
                {
                    params.aLRT_SH_half_epsilon = convert_real_number(argv[cnt]) * 0.5;
                } catch (std::exception e)
                {
                    outError(e.what());
                }
                
                continue;
            }
            if (strcmp(argv[cnt], "-seed") == 0 || strcmp(argv[cnt], "--seed") == 0) {
                cnt++;
                if (cnt >= argc)
                    outError("Use -seed <random_seed>");
                
                try
                {
                    params.ran_seed = abs(convert_int(argv[cnt]));
                }
                catch(std::exception e)
                {
                    outError(e.what());
                }
                continue;
            }
            if (strcmp(argv[cnt], "-nt") == 0 || strcmp(argv[cnt], "-c") == 0 ||
                strcmp(argv[cnt], "-T") == 0  || strcmp(argv[cnt], "--threads") == 0) {
                cnt++;
                if (cnt >= argc)
                    outError("Use -nt <num_threads|AUTO>");
                if (iEquals(argv[cnt], "AUTO"))
                    params.num_threads = 0;
                else {
                    try
                    {
                        params.num_threads = convert_int(argv[cnt]);
                    }
                    catch(std::exception e)
                    {
                        outError(e.what());
                    }
                    if (params.num_threads < 1)
                      outError("At least 1 thread please");
                }
                continue;
            }
            
            // return invalid option
            string err = "Invalid \"";
            err += argv[cnt];
            err += "\" option.";
            outError(err);
        }
        // try
        catch (const char *str) {
                exit(EXIT_SUCCESS);
        } catch (string str) {
                exit(EXIT_SUCCESS);
        } catch (...) {
            string err = "Unknown argument \"";
            err += argv[cnt];
            err += "\"";
            exit(EXIT_SUCCESS);
        }

    }
    
    // validate options
    if (!params.aln_path.length())
        outError("Please supply an alignment file via -aln <ALN_FILENAME>");
    
    if (argc <= 1) {
        quickStartGuide();
    }
}

void cmaple::quickStartGuide() {
    cmaple::printCopyright(cout);
    cout << "Quick Start Guide" << endl;
    exit(0);
}

void cmaple::trimString(string &str) {
    str.erase(0, str.find_first_not_of(" \n\r\t"));
    str.erase(str.find_last_not_of(" \n\r\t")+1);
}

bool cmaple::renameString(string& name) {
    bool renamed = false;
    for (string::iterator i = name.begin(); i != name.end(); i++) {
        if (!isalnum(*i) && (*i) != '_' && (*i) != '-' && (*i) != '.' && (*i) != '|' && (*i) != '/') {
            (*i) = '_';
            renamed = true;
        }
    }
    return renamed;
}

int cmaple::countPhysicalCPUCores() {
    #ifdef _OPENMP
    return omp_get_num_procs();
    #else
    return std::thread::hardware_concurrency();
    #endif
}

Params& Params::getInstance() {
    static Params instance;
    return instance;
}

void cmaple::setNumThreads(const int num_threads)
{
    std::string msg = "";
    // setup the number of threads for openmp
#ifdef _OPENMP
    int max_procs = countPhysicalCPUCores();
    msg += "OpenMP: ";
    if (num_threads >= 1) {
        omp_set_num_threads(num_threads);
        msg += convertIntToString(num_threads) + " threads";
    }
    else // num_threads == 0
    {   // not calling 'omp_set_num_threads' uses all cores automatically
        msg += "auto-detect threads";
    }
    msg += " (" + convertIntToString(max_procs) + " CPU cores detected)";
    if (num_threads  > max_procs) {
        throw std::invalid_argument("\nYou have specified more threads than CPU cores available!");
    }
    #ifndef WIN32  // not supported on Windows (only <=OpenMP2.0)
    omp_set_max_active_levels(1);
    #endif
#else
    if (num_threads != 1) {
        throw std::invalid_argument("\n\nNumber of threads must be 1 for sequential version.");
    }
#endif
    if (cmaple::verbose_mode > cmaple::VB_QUIET)
        std::cout << msg << std::endl;
}

void cmaple::resetStream(std::istream& instream)
{
    instream.clear();
    instream.seekg(0, ios::beg);
}
