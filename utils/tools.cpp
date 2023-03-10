/***************************************************************************
 *   Copyright (C) 2022 by                                            *
 *   BUI Quang Minh <m.bui@anu.edu.au>                                *
 *   Nhan Ly-Trong <trongnhan.uit@gmail.com>                                    *
 *                                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/



#include "tools.h"
#include "timeutil.h"

//#include <filesystem>

using namespace std;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

#define RAN_STANDARD 1
#define RAN_SPRNG    2
#define RAN_RAND4    3

#define RAN_TYPE 2

#if RAN_TYPE == RAN_STANDARD

int init_random(int seed) {
    srand(seed);
    cout << "(Using rand() - Standard Random Number Generator)" << endl;
    // init random generator for AliSim
    Params::getInstance().generator.seed(seed);
    return seed;
}

int finish_random() {
    return 0;
}


#elif RAN_TYPE == RAN_RAND4
/******************************************************************************/
/* random numbers generator  (Numerical recipes)                              */
/******************************************************************************/

/* variable */
long _idum;

/* definitions */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double randomunitintervall()
/* Long period (> 2e18) random number generator. Returns a uniform random
   deviate between 0.0 and 1.0 (exclusive of endpoint values).

   Source:
   Press et al., "Numerical recipes in C", Cambridge University Press, 1992
   (chapter 7 "Random numbers", ran2 random number generator) */ {
    int j;
    long k;
    static long _idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (_idum <= 0) {
        if (-(_idum) < 1)
            _idum = 1;
        else
            _idum = -(_idum);
        _idum2 = (_idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (_idum) / IQ1;
            _idum = IA1 * (_idum - k * IQ1) - k*IR1;
            if (_idum < 0)
                _idum += IM1;
            if (j < NTAB)
                iv[j] = _idum;
        }
        iy = iv[0];
    }
    k = (_idum) / IQ1;
    _idum = IA1 * (_idum - k * IQ1) - k*IR1;
    if (_idum < 0)
        _idum += IM1;
    k = _idum2 / IQ2;
    _idum2 = IA2 * (_idum2 - k * IQ2) - k*IR2;
    if (_idum2 < 0)
        _idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - _idum2;
    iv[j] = _idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
} /* randomunitintervall */

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

int init_random(int seed) /* RAND4 */ {
    //    srand((unsigned) time(NULL));
    //    if (seed < 0)
    //     seed = rand();
    _idum = -(long) seed;
#ifndef PARALLEL
    cout << "(Using RAND4 Random Number Generator)" << endl;
#else /* PARALLEL */
    {
        int n;
        if (PP_IamMaster) {
            cout << "(Using RAND4 Random Number Generator with leapfrog method)" << endl;
        }
        for (n = 0; n < PP_Myid; n++)
            (void) randomunitintervall();
        if (verbose_mode >= VB_MED) {
            cout << "(" << PP_Myid << ") !!! random seed set to " << seed << ", " << n << " drawn !!!" << endl;
        }
    }
#endif
    // init random generator for AliSim
    Params::getInstance().generator.seed(seed);
    return (seed);
} /* initrandom */

int finish_random() {
    return 0;
}
/******************/

#else /* SPRNG */

/******************/

int *randstream;

int init_random(int seed, bool write_info, int** rstream) {
    //    srand((unsigned) time(NULL));
    if (seed < 0)
        seed = make_sprng_seed();
#ifndef PARALLEL
    if (write_info)
        cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    if (rstream) {
        *rstream = init_sprng(0, 1, seed, SPRNG_DEFAULT); /*init stream*/
    } else {
        randstream = init_sprng(0, 1, seed, SPRNG_DEFAULT); /*init stream*/
        if (verbose_mode >= VB_MED) {
            print_sprng(randstream);
        }
    }
#else /* PARALLEL */
    if (PP_IamMaster && write_info) {
        cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    }
    /* MPI_Bcast(&seed, 1, MPI_UNSIGNED, PP_MyMaster, MPI_COMM_WORLD); */
    if (rstream) {
        *rstream = init_sprng(PP_Myid, PP_NumProcs, seed, SPRNG_DEFAULT); /*initialize stream*/
    } else {
        randstream = init_sprng(PP_Myid, PP_NumProcs, seed, SPRNG_DEFAULT); /*initialize stream*/
        if (verbose_mode >= VB_MED) {
            cout << "(" << PP_Myid << ") !!! random seed set to " << seed << " !!!" << endl;
            print_sprng(randstream);
        }
    }
#endif /* PARALLEL */
    return (seed);
} /* initrandom */

int finish_random(int *rstream) {
    if (rstream)
        return free_sprng(rstream);
    else
        return free_sprng(randstream);
}

#endif /* USE_SPRNG */

/******************/

/* returns a random integer in the range [0; n - 1] */
int random_int(int n, int *rstream) {
    return (int) floor(random_double(rstream) * n);
} /* randominteger */

/* returns a random integer in the range [a; b] */
int random_int(int a, int b) {
    ASSERT(b > a);
    //return a + (RAND_MAX * rand() + rand()) % (b + 1 - a);
    return a + random_int(b - a);
}

double random_double(int *rstream) {
#ifndef FIXEDINTRAND
#ifndef PARALLEL
#if RAN_TYPE == RAN_STANDARD
    return ((double) rand()) / ((double) RAND_MAX + 1);
#elif RAN_TYPE == RAN_SPRNG
    if (rstream)
        return sprng(rstream);
    else
        return sprng(randstream);
#else /* NO_SPRNG */
    return randomunitintervall();
#endif /* NO_SPRNG */
#else /* NOT PARALLEL */
#if RAN_TYPE == RAN_SPRNG
    if (rstream)
        return sprng(rstream);
    else
        return sprng(randstream);
#else /* NO_SPRNG */
    int m;
    for (m = 1; m < PP_NumProcs; m++)
        (void) randomunitintervall();
    PP_randn += (m - 1);
    PP_rand++;
    return randomunitintervall();
#endif /* NO_SPRNG */
#endif /* NOT PARALLEL */
#else /* FIXEDINTRAND */
    cerr << "!!! fixed \"random\" integers for testing purposes !!!" << endl;
    return 0.0;
#endif /* FIXEDINTRAND */

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

VerboseMode verbose_mode;

void printCopyright(ostream &out) {
     out << "CMAPLE";
}

/**
        Output an error to screen, then exit program
        @param error error message
 */
void outError(const char *error, bool quit) {
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
void outError(const string &error, bool quit) {
    outError(error.c_str(), quit);
}

void outError(const char *error, const char *msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

void outError(const char *error, const string &msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

/**
        Output a warning message to screen
        @param error warning message
 */
void outWarning(const char *warn) {
    cout << "WARNING: " << warn << endl;
}

void outWarning(const string &warn) {
    outWarning(warn.c_str());
}

std::istream& safeGetline(std::istream& is, std::string& t)
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
string convertPosTypeToString(PositionType number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertIntToString(int number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertInt64ToString(int64_t number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertDoubleToString(RealNumType number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

bool iEquals(const string &a, const string &b)
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

bool copyFile(const char SRC[], const char DEST[]) {
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

bool fileExists(const string &strFilename) {
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

int convert_int(const char *str, int &end_pos) {
	char *endptr;
	int i = strtol(str, &endptr, 10);

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		outError(err);
	}
	end_pos = endptr - str;
	return i;
}

int convert_int(const char *str) {
    char *endptr;
    int i = strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        outError(err);
    }

    return i;
}

PositionType convert_positiontype(const char *str) {
    char *endptr;
    PositionType i = (PositionType) strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || i == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        outError(err);
    }

    return i;
}

void convert_int_vec(const char *str, IntVector &vec) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		int i = strtol(beginptr, &endptr, 10);

		if ((i == 0 && endptr == beginptr) || abs(i) == HUGE_VALL) {
			string err = "Expecting integer, but found \"";
			err += beginptr;
			err += "\" instead";
            outError(err);
		}
		vec.push_back(i);
		if (*endptr == ',') endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}


int64_t convert_int64(const char *str) {
    char *endptr;
    int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting large integer , but found \"";
        err += str;
        err += "\" instead";
        outError(err);
    }

    return i;
}

int64_t convert_int64(const char *str, int &end_pos) {
	char *endptr;
	int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting large integer, but found \"";
		err += str;
		err += "\" instead";
        outError(err);
	}
	end_pos = endptr - str;
	return i;
}


RealNumType convert_real_number(const char *str) {
    char *endptr;
    RealNumType d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        outError(err);
    }
    return d;
}

RealNumType convert_real_number(const char *str, int &end_pos) {
	char *endptr;
	RealNumType d = strtod(str, &endptr);
	if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
        outError(err);
	}
	end_pos = endptr - str;
	return d;
}

void convert_real_numbers(RealNumType* &arr, string input_str)
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

void convert_real_number_vec(const char *str, RealNumberVector &vec, char separator) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		RealNumType d = strtod(beginptr, &endptr);

		if ((d == 0.0 && endptr == beginptr) || fabs(d) == HUGE_VALF) {
			string err = "Expecting floating-point number, but found \"";
			err += beginptr;
			err += "\" instead";
            outError(err);
		}
		vec.push_back(d);
		if (*endptr == separator) endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}

string convert_time(const RealNumType sec) {
    int sec_int = (int) floor(sec);
    int secs = sec_int % 60;
    int mins = (sec_int % 3600) / 60;
    int hours = sec_int / 3600;
    stringstream ss;
    ss << hours << "h:" << mins << "m:" << secs << "s";
    return ss.str();
}

void convert_range(const char *str, int &lower, int &upper, int &step_size) {
    char *endptr;

    // parse the lower bound of the range
    int d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        outError(err);
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
        outError(err);
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
        outError(err);
    }
    step_size = d;
}

void convert_range(const char *str, RealNumType &lower, RealNumType &upper, RealNumType &step_size) {
    char *endptr;

    // parse the lower bound of the range
    RealNumType d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        outError(err);
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
        outError(err);
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
        outError(err);
    }
    step_size = d;
}

void reinitDoubleArr(RealNumType* &arr, StateType size, bool delete_first, bool set_zero)
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

void convert_string_vec(const char *str, StrVector &vec, char separator) {
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

void normalize_frequencies_from_index(RealNumType* freqs, int num_states, int starting_index)
{
    ASSERT(num_states > 0);
    // calculate the total_freqs
    RealNumType total_freqs = 0;
    for (int i = starting_index; i < starting_index+num_states; ++i)
        total_freqs += freqs[i];
    
    // normalize the freqs
    if (fabs(total_freqs) < 1e-5)
        outError("Sum of state frequencies must be greater than zero!");
    total_freqs = 1.0 / total_freqs;
    for (int i = starting_index; i < starting_index+num_states; ++i)
        freqs[i] *= total_freqs;
}

bool is_number(const std::string& s)
{
    char* end = nullptr;
    double val = strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && val != HUGE_VAL;
}

void quickStartGuide();

void initDefaultValue(Params &params)
{
    params.aln_path = NULL;
    params.diff_path = NULL;
    params.ref_path = NULL;
    params.only_extract_diff = false;
    params.hamming_weight = 1000;
    params.model_name = "GTR";
    params.redo_inference = false;
    params.threshold_prob = 1e-8;
    params.mutation_update_period = 25;
    params.failure_limit_sample = 5;
    params.failure_limit_subtree = 4;
    params.failure_limit_subtree_short_search = 1;
    params.strict_stop_seeking_placement_sample = false;
    params.strict_stop_seeking_placement_subtree = false;
    params.strict_stop_seeking_placement_subtree_short_search = true;
    params.thresh_log_lh_sample = 200;
    params.thresh_log_lh_subtree = 160;
    params.thresh_log_lh_subtree_short_search = 40;
    params.thresh_log_lh_failure = 0.01;
    params.min_blength_factor = 0.2;
    params.min_blength_mid_factor = 4.1;
    params.max_blength_factor = 40;
    params.thresh_diff_update = 1e-7;
    params.thresh_diff_fold_update = 1.001;
    params.output_aln = NULL;
    params.num_tree_improvement = 1;
    params.thresh_entire_tree_improvement = 1;
    params.thresh_placement_cost = -1e-5;
    params.thresh_placement_cost_short_search = -1;
    params.export_binary_tree = true;
    params.optimize_branch_length = true;
    params.short_range_topo_search = false;
    params.output_testing = NULL;
    params.compute_aLRT_SH = false;
    params.aLRT_SH_replicates = 10000;
    params.aLRT_SH_epsilon = 0.1;
    params.num_threads = 1;
    params.input_treefile = NULL;
    
    // initialize random seed based on current time
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    //params.ran_seed = (unsigned) (tv.tv_sec+tv.tv_usec);
    params.ran_seed = (tv.tv_usec);
}

void parseArg(int argc, char *argv[], Params &params) {
    // init parameters
    initDefaultValue(params);
    
    for (int cnt = 1; cnt < argc; ++cnt) {
        try {
            if (strcmp(argv[cnt], "--aln") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --aln <ALIGNMENT_PATH>");
                
                params.aln_path = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--diff") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --diff <DIFF_PATH>");
                
                params.diff_path = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--output-aln") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --output-aln <ALIGNMENT_PATH>");
                
                params.output_aln = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--tree") == 0 || strcmp(argv[cnt], "-t") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use -t <INPUT_TREEFILE>");
                
                params.input_treefile = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--ref") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --ref <REF_PATH>");
                
                params.ref_path = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--extract-diff") == 0) {
                
                params.only_extract_diff = true;

                continue;
            }
            if (strcmp(argv[cnt], "--hamming-weight") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --hamming-weight <WEIGHT>");
                
                params.hamming_weight = convert_real_number(argv[cnt]);
                
                if (params.hamming_weight < 0)
                    outError("<WEIGHT> must not be negative!");

                continue;
            }
            if (strcmp(argv[cnt], "--model") == 0 || strcmp(argv[cnt], "-m") == 0) {
                ++cnt;
                if (cnt >= argc)
                    outError("Use --model <model_name>");
                
                params.model_name = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "-redo") == 0 || strcmp(argv[cnt], "--redo") == 0) {
                params.redo_inference = true;
                continue;
            }
            if (strcmp(argv[cnt], "--thresh-prob") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --thresh-prob <PROB_THRESH>");
                
                params.threshold_prob = convert_real_number(argv[cnt]);
                
                if (params.threshold_prob <= 0)
                    outError("<PROB_THRESH> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--mutation-update") == 0) {
                
                ++cnt;
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --mutation-update <NUMBER>");
                
                params.mutation_update_period = convert_int(argv[cnt]);
                
                if (params.mutation_update_period <= 0)
                    outError("<NUMBER> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--failure-limit") == 0) {
                
                ++cnt;
                
                params.failure_limit_sample = convert_int(argv[cnt]);
                
                if (params.failure_limit_sample <= 0)
                    outError("<NUMBER> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--failure-limit-subtree") == 0) {
                
                ++cnt;
                
                params.failure_limit_subtree = convert_int(argv[cnt]);
                
                if (params.failure_limit_subtree <= 0)
                    outError("<NUMBER> must be positive!");

                continue;
            }
            if (strcmp(argv[cnt], "--strict-stop-init") == 0) {
                
                params.strict_stop_seeking_placement_sample = true;

                continue;
            }
            if (strcmp(argv[cnt], "--unstrict-stop-subtree") == 0) {
                
                params.strict_stop_seeking_placement_subtree = false;

                continue;
            }
            if (strcmp(argv[cnt], "--multifurcating-tree") == 0) {
                
                params.export_binary_tree = false;

                continue;
            }
            if (strcmp(argv[cnt], "--no-optimize-blength") == 0) {
                
                params.optimize_branch_length = false;

                continue;
            }
            if (strcmp(argv[cnt], "--short-topo-search") == 0) {
                
                params.short_range_topo_search = true;

                continue;
            }
            if (strcmp(argv[cnt], "--output-testing") == 0) {
                
                ++cnt;
                
                if (cnt >= argc || argv[cnt][0] == '-')
                    outError("Use --output-testing <FILE_PATH>");
                
                params.output_testing = argv[cnt];

                continue;
            }
            if (strcmp(argv[cnt], "--branch-support") == 0) {
                
                params.compute_aLRT_SH = true;

                continue;
            }
            if (strcmp(argv[cnt], "--replicates") == 0) {
                ++cnt;
                if (cnt >= argc)
                    outError("Use --replicates <NUM_REPLICATES>");
                
                params.aLRT_SH_replicates = convert_int(argv[cnt]);
                
                if (params.aLRT_SH_replicates <= 0)
                    outError("<NUM_REPLICATES> must be positive!");
                continue;
            }
            if (strcmp(argv[cnt], "--epsilon") == 0) {
                ++cnt;
                if (cnt >= argc)
                    outError("Use --epsilon <FLOATING_NUM>");
                
                params.aLRT_SH_epsilon = convert_real_number(argv[cnt]);
                
                continue;
            }
            if (strcmp(argv[cnt], "-seed") == 0 || strcmp(argv[cnt], "--seed") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -seed <random_seed>";
                params.ran_seed = abs(convert_int(argv[cnt]));
                continue;
            }
            if (strcmp(argv[cnt], "-nt") == 0 || strcmp(argv[cnt], "-c") == 0 ||
                strcmp(argv[cnt], "-T") == 0  || strcmp(argv[cnt], "--threads") == 0) {
                cnt++;
                if (cnt >= argc)
                throw "Use -nt <num_threads|AUTO>";
                if (iEquals(argv[cnt], "AUTO"))
                    params.num_threads = 0;
                else {
                    params.num_threads = convert_int(argv[cnt]);
                    if (params.num_threads < 1)
                        throw "At least 1 thread please";
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
    if (!params.diff_path && !params.aln_path)
        outError("Please supply either an alignment or a Diff file to start!");
        
    if (params.only_extract_diff && !params.aln_path)
        outError("Please supply an input alignment via --aln <ALIGNMENT_PATH>");
    
    if (argc <= 1) {
        quickStartGuide();
    }
}

void quickStartGuide() {
    printCopyright(cout);
    cout << "Quick Start Guide" << endl;
    exit(0);
}

InputType detectInputFile(const char *input_file) {

    if (!fileExists(input_file))
        outError("File not found ", input_file);

    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(input_file);

        unsigned char ch = ' ';
        unsigned char ch2 = ' ';
        int count = 0;
        do {
            in >> ch;
        } while (ch <= 32 && !in.eof() && count++ < 20);
        in >> ch2;
        in.close();
        switch (ch) {
            case '#': return IN_NEXUS;
            case '(': return IN_NEWICK;
            case '[': return IN_NEWICK;
            case '>': return IN_FASTA;
            case 'C': if (ch2 == 'L') return IN_CLUSTAL;
                      else if (ch2 == 'O') return IN_COUNTS;
                      else return IN_OTHER;
            case '!': if (ch2 == '!') return IN_MSF; else return IN_OTHER;
            default:
                if (isdigit(ch)) return IN_PHYLIP;
                return IN_OTHER;
        }
    } catch (ios::failure const&) {
        outError("Cannot read file ", input_file);
    } catch (...) {
        outError("Cannot read file ", input_file);
    }
    return IN_OTHER;
}

bool overwriteFile(char *filename) {
    ifstream infile(filename);
    if (infile.is_open()) {
        cout << "Overwrite " << filename << " (y/n)? ";
        char ch;
        cin >> ch;
        if (ch != 'Y' && ch != 'y') {
            infile.close();
            return false;
        }
    }
    infile.close();
    return true;
}

void trimString(string &str) {
    str.erase(0, str.find_first_not_of(" \n\r\t"));
    str.erase(str.find_last_not_of(" \n\r\t")+1);
}

bool renameString(string& name) {
    bool renamed = false;
    for (string::iterator i = name.begin(); i != name.end(); i++) {
        if (!isalnum(*i) && (*i) != '_' && (*i) != '-' && (*i) != '.' && (*i) != '|' && (*i) != '/') {
            (*i) = '_';
            renamed = true;
        }
    }
    return renamed;
}

int countPhysicalCPUCores() {
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
