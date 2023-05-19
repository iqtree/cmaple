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
    std::cout << params.aLRT_SH_replicates << std::endl;
    //runCMaple(params);
}


