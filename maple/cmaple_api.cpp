#include "cmaple_api.h"
#include "../utils/tools.h"
#include "../maple/cmaple.h"

using namespace cmaple;

void cmaple::printCMapleCopyright(std::ostream &out)
{
    printCopyright(out);
}

void cmaple::executeCMaple(char* input_file)
{
    Params& params = Params::getInstance();
    initDefaultValue(params);
    params.input_path = input_file;
    runCMaple(params);
}


