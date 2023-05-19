#pragma once
#include <streambuf>
#include <fstream>
#include <iostream>
#include "tools.h"

/** Output stream buffer */
class outstreambuf : public std::streambuf {
public:
    /**
        Open the log file
     */
    outstreambuf* open(const char* name, std::ios::openmode mode = std::ios::out);
    
    /**
        TRUE if the log file is open
     */
    bool is_open();
    
    /**
        Close the outstreambuf
     */
    outstreambuf* close();
    
    /**
        Destructor
     */
    ~outstreambuf() { close(); }
    
    /**
        Get file output buffer?
     */
    std::streambuf* get_fout_buf() {
        return fout_buf;
    }
    
    /**
        Get cout buffer
     */
    std::streambuf* get_cout_buf() {
        return cout_buf;
    }
    
    /**
        Get the ofstream
     */
    std::ofstream* get_fout() {
        return &fout;
    }
    
protected:
    std::ofstream fout;
    std::streambuf* cout_buf;
    std::streambuf* fout_buf;
    virtual int overflow( int c = EOF);
    virtual int sync();
};

/** Error stream buffer */
class errstreambuf : public std::streambuf {
public:
    /**
        Initialize the errstreambuf
     */
    void init(std::streambuf* fout_buf) {
        this->fout_buf = fout_buf;
        cerr_buf = std::cerr.rdbuf();
        std::cerr.rdbuf(this);
        new_line = true;
    }
    
    /**
        Destructor
     */
    ~errstreambuf() {
        std::cerr.rdbuf(cerr_buf);
    }
    
protected:
    std::streambuf* cerr_buf;
    std::streambuf* fout_buf;
    bool new_line;
    
    virtual int overflow( int c = EOF) {
        if (new_line)
            cerr_buf->sputn("ERROR: ", 7);
        if (cerr_buf->sputc(c) == EOF) {
            new_line = false;
            if (c == '\n') new_line = true;
            return EOF;
        }
        /*if ((Params::getInstance().suppress_output_flags & OUT_LOG)) {
            new_line = false;
            if (c == '\n') new_line = true;
            return c;
        }*/
        if (new_line)
            fout_buf->sputn("ERROR: ", 7);
        new_line = false;
        if (c == '\n') new_line = true;
        if (fout_buf->sputc(c) == EOF) return EOF;
        return c;
    }
    
    virtual int sync() {
        cerr_buf->pubsync();
        /*if (Params::getInstance().suppress_output_flags & OUT_LOG)
            return 0;*/
        return fout_buf->pubsync();
    }
};

/** Must stream buffer ? */
class muststreambuf : public std::streambuf {
public:
    void init(std::streambuf* cout_buf, std::streambuf* fout_buf) {
        this->fout_buf = fout_buf;
        this->cout_buf = cout_buf;
    }
    
protected:
    std::streambuf* cout_buf;
    std::streambuf* fout_buf;
    
    virtual int overflow( int c = EOF) {
        if (cout_buf->sputc(c) == EOF) {
            return EOF;
        }
        if (fout_buf->sputc(c) == EOF) return EOF;
        return c;
    }
    
    virtual int sync() {
        cout_buf->pubsync();
        return fout_buf->pubsync();
    }
};


/* ################# Redirect output to log file ###################### */
class LogStream
{
private:
    outstreambuf _out_buf;
    errstreambuf _err_buf;
    muststreambuf _must_buf;
    std::string log_file_;
    
public:
    void startLogFile(cmaple::Params& params);
    void endLogFile();
    void funcExit(void);
    std::string getLogFilePath();
};
