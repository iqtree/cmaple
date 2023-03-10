/***************************************************************************
 *   Copyright (C) 2022 by BUI Quang Minh, Nhan Ly-Trong  *
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmaple_config.h>

#include "utils/timeutil.h"
#include "utils/tools.h"
#include "utils/operatingsystem.h" //for getOSName()
#include "maple/cmaple.h"

#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main(int argc, char *argv[]) {
    parseArg(argc, argv, Params::getInstance());
    
    cout << "Command:";
    for (int i = 0; i < argc; i++)
        cout << " " << argv[i];
    cout << endl;
    
    // Show info
    cout << "Seed:    " << Params::getInstance().ran_seed <<  " ";
    init_random(Params::getInstance().ran_seed, true);
    
    // setup the number of threads for openmp
#ifdef _OPENMP
    if (Params::getInstance().num_threads >= 1) {
        omp_set_num_threads(Params::getInstance().num_threads);
        Params::getInstance().num_threads = omp_get_max_threads();
    }
//    int max_threads = omp_get_max_threads();
    int max_procs = countPhysicalCPUCores();
    cout << "OpenMP: ";
    if (Params::getInstance().num_threads > 0)
        cout << Params::getInstance().num_threads  << " threads";
    else
        cout << "auto-detect threads";
    cout << " (" << max_procs << " CPU cores detected) \n";
    if (Params::getInstance().num_threads  > max_procs) {
        cout << endl;
        outError("You have specified more threads than CPU cores available");
    }
    omp_set_max_active_levels(1);
#else
    if (Params::getInstance().num_threads != 1) {
        cout << endl << endl;
        outError("Number of threads must be 1 for sequential version.");
    }
#endif
    
    // Measure runtime
    time_t start_time;
    
    // call the main function
    runCMaple(Params::getInstance());
    
    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    
    return EXIT_SUCCESS;
}
