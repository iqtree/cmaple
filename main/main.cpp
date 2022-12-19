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

#include <stdio.h>
#include "utils/timeutil.h"
#include "utils/tools.h"
#include "utils/operatingsystem.h" //for getOSName()
#include <stdlib.h>
#include "cmaple.h"
#include "unitest/maintest.cpp"

using namespace std;

void printCopyright(ostream &out) {
     out << "CMAPLE";
}


int main(int argc, char *argv[]) {
    parseArg(argc, argv, Params::getInstance());
    
    // GTest
    testing::InitGoogleTest(&argc, (char**)argv);
    RUN_ALL_TESTS();
    
    // Measure runtime
    time_t start_time;
    
    // call the main function
    runCMaple(Params::getInstance());
    
    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    
    return EXIT_SUCCESS;
}
