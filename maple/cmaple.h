#include "../tree/tree.h"

#pragma once
namespace cmaple
{
    /*!
     *  \addtogroup cmaple
     *  @{
     */

    /*! \brief Get the version of CMaple
     * @return CMAPLE version in a string
     */
    std::string getVersion();
    
    /*! \brief Get the citation strings for CMaple/Maple manuscripts
     * @return A citation string
     */
    std::string getCitations();

    /*! \brief Check if an alignment is suitable for using [C]Maple method
     * @param[in] aln an alignment
     * @return TRUE if the alignment is suitable for using [C]Maple method; otherwise, classical methods (e.g., IQ-TREE, RAXML) are recommended.
     * @throw std::invalid\_argument if the input alignment is empty or invalid
     */
    bool checkMapleSuitability(const Alignment& aln);

    /*! @} End of Doxygen Groups*/

    /** \brief Run CMaple with user-specified parameters (via command-line)
     * @private
     * @param[in] params user-specified parameters (via command-line)
     */
    void runCMaple(cmaple::Params& params);

    /** \brief Function for testing only
     * @private
     * @param[in] params user-specified (parameters via command-line)
     */
    void testing(cmaple::Params& params);
}
