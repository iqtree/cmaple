#include "../tree/treewrapper.h"

#pragma once
namespace cmaple
{
    /*! \brief Get the version of CMaple
     * @return CMAPLE version in a string
     */
    std::string getVersion();
    
    /*! \brief Get the citation strings for CMaple/Maple manuscripts
     * @return A citation string
     */
    std::string getCitations();

    /** \brief Run CMaple with user-specified parameters (via command-line)
     * @param[in] params user-specified parameters (via command-line)
     */
    void runCMaple(cmaple::Params& params);

    /** \brief Function for testing only
     * @param[in] params user-specified (parameters via command-line)
     */
    void testing(cmaple::Params& params);
}
