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
#include "../tree/tree.h"

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
    // NHANLT TEST
    /*Model model("JC");
    std::map<std::string,std::string> model_params = model.getParams();
    std::cout << model_params[MODEL_NAME] << std::endl;
    std::cout << model_params[MODEL_FREQS] << std::endl;
    std::cout << model_params[MODEL_RATES] << std::endl;*/
    
    cmaple::Params& params = Params::getInstance();
    parseArg(argc, argv, params);
    
    // Config log file
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
    
    // Show the random seed number
    /*cout << "Seed:    " << params.ran_seed <<  " " << std::endl;
    
    // Show the number of threads
    setNumThreads(params.num_threads);*/
    
    // Measure runtime
    time_t start_time;
    
    // call the main function
    /*CMaple cmaple;
    cmaple::Params& cmaple_params = cmaple.getSettings();
    cmaple_params = params;
    cmaple.inferTree();*/
    Model model(params.model_name);
    // with an alignment file
    Alignment aln1(params.aln_path);
    // with an alignment file
    // Alignment aln2(""); // tested PASS
    // Alignment aln3("notfound"); // tested PASS
    // with a stream
    const std::string aln_filename = params.aln_path;
    std::ifstream aln_stream;
    try {
        aln_stream.exceptions(ios::failbit | ios::badbit);
        aln_stream.open(aln_filename);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, aln_filename);
    }
    Alignment aln(aln_stream);
    aln_stream.close();
    // Write alignment to file in MAPLE format
    aln.write("output.maple", "MAPLE", true);
    // write aln to a file in FASTA format
    aln.write("output.fa", "FASTA", true);
    aln.write("output.phy", "PHYLIP", true);
    // Read ref_seq from an alignment file (not yet exposed to APIs)
    ASSERT(params.ref_path.length() && params.ref_seqname.length());
    AlignmentBase aln_base;
    std::string ref_seq = aln_base.readRefSeq(params.ref_path, params.ref_seqname);
    Alignment aln2("output.fa", ref_seq);
    aln2.write("output1.maple", "MAPLE", true);
    Alignment aln3("output.fa");
    aln3.write("output2.maple", "MAPLE", true);
    Alignment aln4("output.phy", ref_seq);
    aln4.write("output_phy1.maple", "MAPLE", true);
    Alignment aln5("output.phy");
    aln5.write("output_phy2.maple", "MAPLE", true);
    // aln5.write("output_phy3.maple", "INVALID", true);
    // aln.write("output.maple");
    // without tree file
    Tree tree1(aln, model);
    // with a tree file
    const std::string tree_filename = params.input_treefile;
    Tree tree2(aln, model, tree_filename);
    // with a tree stream
    std::ifstream tree_stream;
    try {
        tree_stream.exceptions(ios::failbit | ios::badbit);
        tree_stream.open(tree_filename);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, tree_filename);
    }
    Tree tree(aln, model, tree_stream);
    tree_stream.close();
    tree.infer();
    
    std::cout << tree.infer();

    // std::cout << tree.exportString("BIN", true) << std::endl;
    //Tree tree(aln, model, "");
    std::cout << "Tree likelihood: " << tree.computeLh() << std::endl;
    tree.computeBranchSupports(8, 100);
    std::cout << tree.exportString("BIN", true) << std::endl;
    tree.infer();
    std::cout << tree.exportString() << std::endl;
    tree.computeBranchSupports(8, 100);
    std::cout << tree.exportString("BIN", true) << std::endl;
    std::cout << "Tree likelihood: " << tree.computeLh() << std::endl;
    
    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    
    return EXIT_SUCCESS;
}
