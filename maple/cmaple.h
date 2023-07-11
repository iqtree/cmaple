#include "../tree/tree.h"

#pragma once
namespace cmaple
{
    /*! \brief Get the version of CMaple
     *
     * Get the version of CMaple
     * @return CMAPLE version
     */
    std::string getVersion();
    
    /*! \brief Get the citation string
     *
     * Get the citation strings for CMaple/Maple manuscripts
     * @return citation strings
     */
    std::string getCitations();

    /** \brief Run CMaple with user-specified (via command-line) parameters
     * @param[in] params user-specified (via command-line) parameters
     */
    void runCMaple(cmaple::Params& params);

    /** \brief Function for testing only
     * @param[in] params user-specified (via command-line) parameters
     */
    void testing(cmaple::Params& params);
}
