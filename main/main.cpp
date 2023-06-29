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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmaple_config.h>

#include "../utils/timeutil.h"
#include "../utils/tools.h"
#include "../utils/operatingsystem.h" //for getOSName()
#include "../utils/logstream.h"
#include "../maple/cmaple.h"

#include <stdlib.h>
#include <stdio.h>

#include <csignal>

using namespace std;
using namespace cmaple;

/** ############## Redirect output to a log file ############## **/
LogStream logstream;
void funcAbort(int signal_number)
{
    /*Your code goes here. You can output debugging info.
      If you return from this function, and it was called
      because abort() was called, your program will exit or crash anyway
      (with a dialog box on Windows).
     */
#if (defined(__GNUC__) || defined(__clang__)) && !defined(WIN32) && !defined(WIN64) && !defined(__CYGWIN__)
    // print_stacktrace(cerr);
#endif

    cerr << endl << "*** CMAPLE CRASHES WITH SIGNAL ";
    switch (signal_number) {
        case SIGABRT: cerr << "ABORTED"; break;
        case SIGFPE:  cerr << "ERRONEOUS NUMERIC"; break;
        case SIGILL:  cerr << "ILLEGAL INSTRUCTION"; break;
        case SIGSEGV: cerr << "SEGMENTATION FAULT"; break;
#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__ && !defined WIN64
        case SIGBUS: cerr << "BUS ERROR"; break;
#endif
    }
    cerr << endl;
    cerr << "*** For bug report please send to developers:" << endl << "***    Log file: " << logstream.getLogFilePath();
    cerr << endl << "***    Alignment files (if possible)" << endl;
    logstream.funcExit();
    signal(signal_number, SIG_DFL);
}
/** ############## Redirect output to a log file ############## **/

int main(int argc, char *argv[]) {
    parseArg(argc, argv, Params::getInstance());
    
    // Config log file
    Params& params = Params::getInstance();
    // atexit(logstream.funcExit);
    logstream.startLogFile(params);
    signal(SIGABRT, &funcAbort);
    signal(SIGFPE, &funcAbort);
    signal(SIGILL, &funcAbort);
    signal(SIGSEGV, &funcAbort);
#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__ && !defined WIN64
    signal(SIGBUS, &funcAbort);
#endif
    
    // print copyright
    printCopyright(cout);
    
    // show general information
    cout << "Command:";
    for (int i = 0; i < argc; i++)
        cout << " " << argv[i];
    cout << endl;
    
    // Show info
    cout << "Seed:    " << Params::getInstance().ran_seed <<  " " << std::endl;
    
    // setup the number of threads for openmp
#ifdef _OPENMP
    int max_procs = countPhysicalCPUCores();
    cout << "OpenMP: ";
    if (Params::getInstance().num_threads >= 1) {
        omp_set_num_threads(Params::getInstance().num_threads);
        cout << Params::getInstance().num_threads << " threads";
    }
    else // num_threads == 0
    {   // not calling 'omp_set_num_threads' uses all cores automatically
        cout << "auto-detect threads";
    }
    cout << " (" << max_procs << " CPU cores detected) \n";
    if (Params::getInstance().num_threads  > max_procs) {
        cout << endl;
        outError("You have specified more threads than CPU cores available");
    }
    #ifndef WIN32  // not supported on Windows (only <=OpenMP2.0)
    omp_set_max_active_levels(1);
    #endif
#else
    if (Params::getInstance().num_threads != 1) {
        cout << endl << endl;
        outError("Number of threads must be 1 for sequential version.");
    }
#endif
    std::cout << std::endl;
    
    // Measure runtime
    time_t start_time;
    
    // call the main function
    runCMaple(Params::getInstance());
    
    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    
    return EXIT_SUCCESS;
}
