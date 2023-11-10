#include "logstream.h"
using namespace std;

auto outstreambuf::open(const char *name, ios::openmode mode)
    -> outstreambuf * {
  // if (!(Params::getInstance().suppress_output_flags & OUT_LOG) &&
  // MPIHelper::getInstance().isMaster()) {
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

auto outstreambuf::is_open() -> bool { return fout.is_open(); }

auto outstreambuf::close() -> outstreambuf * {
  cout.rdbuf(cout_buf);
  if (fout.is_open()) {
    sync();
    fout.close();
    return this;
  }
  return NULL;
}

auto outstreambuf::overflow(int c) -> int { // used for output buffer only
  /*if ((verbose_mode >= VB_MIN && MPIHelper::getInstance().isMaster()) ||
  verbose_mode >= VB_MED) if (cout_buf->sputc(c) == EOF) return EOF; if
  (Params::getInstance().suppress_output_flags & OUT_LOG) return c; if
  (!MPIHelper::getInstance().isMaster()) return c; if (fout_buf->sputc(c) ==
  EOF) return EOF; return c;*/
  if (cout_buf->sputc(c) == EOF) {
    return EOF;
  }
  if (fout_buf->sputc(c) == EOF) {
    return EOF;
  }
  return c;
}

auto outstreambuf::sync() -> int { // used for output buffer only
  /*if ((verbose_mode >= VB_MIN && MPIHelper::getInstance().isMaster()) ||
  verbose_mode >= VB_MED) cout_buf->pubsync(); if
  ((Params::getInstance().suppress_output_flags & OUT_LOG) ||
  !MPIHelper::getInstance().isMaster()) return 0; return fout_buf->pubsync();*/

  cout_buf->pubsync();
  return fout_buf->pubsync();
}

/**##################################################**/

void LogStream::startLogFile(cmaple::Params& params) {
  // Initialize the name of the log file -> use aln_ as the output prefix if
  // users didn't specify it
  const std::string prefix =
      (params.output_prefix.length() ? params.output_prefix : params.aln_path);
  log_file_ = prefix + ".log";

  _out_buf.open(log_file_.c_str());
  _err_buf.init(_out_buf.get_fout_buf());
  _must_buf.init(_out_buf.get_cout_buf(), _out_buf.get_fout_buf());
}

void LogStream::endLogFile() {
  if (_out_buf.is_open()) {
    _out_buf.close();
  }
}

void LogStream::funcExit() {
  /*if(_exit_wait_optn) {
      printf("\npress [return] to finish: ");
      fflush(stdout);
      while (getchar() != '\n');
  }*/

  endLogFile();
  // MPIHelper::getInstance().finalize();
}

auto LogStream::getLogFilePath() -> string { return log_file_; }
