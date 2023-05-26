#pragma once
#include <ostream>
#include "cmaple.h"

/** The namespace of CMaple library */
namespace cmaple {

    /** Print the copyright of CMaple */
    void printCMapleCopyright(std::ostream &out);

    /** execute CMaple from an input file (alignment of Maple file)*/
    void executeCMaple(char* input_file);
}
