#include "logstream.h"
using namespace std;

outstreambuf* outstreambuf::open( const char* name, ios::openmode mode) {
    // if (!(Params::getInstance().suppress_output_flags & OUT_LOG) && MPIHelper::getInstance().isMaster()) {
        fout.open(name, mode);
        if (!fout.is_open()) {
            cerr << "ERROR: Could not open " << name << " for logging" << endl;
            exit(EXIT_FAILURE);
            return NULL;
        }
        fout_buf = fout.rdbuf();
    // }
    
    cout_buf = cout.rdbuf();
    cout.rdbuf(this);
    return this;
}

bool outstreambuf::is_open() {
    return fout.is_open();
}

outstreambuf* outstreambuf::close() {
    cout.rdbuf(cout_buf);
    if ( fout.is_open()) {
        sync();
        fout.close();
        return this;
    }
    return NULL;
}

int outstreambuf::overflow( int c) { // used for output buffer only
    /*if ((verbose_mode >= VB_MIN && MPIHelper::getInstance().isMaster()) || verbose_mode >= VB_MED)
        if (cout_buf->sputc(c) == EOF) return EOF;
    if (Params::getInstance().suppress_output_flags & OUT_LOG)
        return c;
    if (!MPIHelper::getInstance().isMaster())
        return c;
    if (fout_buf->sputc(c) == EOF) return EOF;
    return c;*/
    if (cout_buf->sputc(c) == EOF) return EOF;
    if (fout_buf->sputc(c) == EOF) return EOF;
    return c;
}

int outstreambuf::sync() { // used for output buffer only
    /*if ((verbose_mode >= VB_MIN && MPIHelper::getInstance().isMaster()) || verbose_mode >= VB_MED)
        cout_buf->pubsync();
    if ((Params::getInstance().suppress_output_flags & OUT_LOG) || !MPIHelper::getInstance().isMaster())
        return 0;
    return fout_buf->pubsync();*/
    
    cout_buf->pubsync();
    return fout_buf->pubsync();
}

/**##################################################**/

void LogStream::startLogFile(cmaple::Params& params) {
    // Initialize the name of the log file -> use aln_ as the output prefix if users didn't specify it
    const std::string prefix = (params.output_prefix.length() ? params.output_prefix :  params.aln_path);
    log_file_ = prefix + ".log";
    
    _out_buf.open(log_file_.c_str());
    _err_buf.init(_out_buf.get_fout_buf());
    _must_buf.init(_out_buf.get_cout_buf(), _out_buf.get_fout_buf());
}

void LogStream::endLogFile() {
    if (_out_buf.is_open())
        _out_buf.close();
}

void LogStream::funcExit(void) {
    /*if(_exit_wait_optn) {
        printf("\npress [return] to finish: ");
        fflush(stdout);
        while (getchar() != '\n');
    }*/
    
    endLogFile();
    // MPIHelper::getInstance().finalize();
}

string LogStream::getLogFilePath()
{
    return log_file_;
}
