#include "cmaple_api.h"
#include "../utils/tools.h"

using namespace cmaple;

void cmaple::printCMapleCopyright(std::ostream &out)
{
    printCopyright(out);
}

void cmaple::testCMaple()
{
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.redo_inference = true;
    params.overwrite_output = true;
    params.redo_blength = true;
    string diff_file = "test_5K.diff";
    params.diff_path = new char[diff_file.length() + 1];
    strcpy(params.diff_path, diff_file.c_str());
    runCMaple(params);
}


