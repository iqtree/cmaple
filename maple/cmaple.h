#include "../tree/tree.h"

#pragma once
namespace cmaple
{
    /*!
     *  \addtogroup cmaple
     *  @{
     */

    /*! \brief Get the version of CMAPLE
     * @return CMAPLE version in a string
     */
    std::string getVersion();
    
    /*! \brief Get a citation string of [C]Maple manuscript(s)
     * @return A citation string
     */
    std::string getCitations();

    /*! \brief Check if the [C]Maple algorithm is effective to analyse an input alignment
     * @param[in] aln an alignment
     * @return TRUE if the [C]Maple algorithm is effective to analyse the alignment;
     * otherwise, classical methods (e.g., IQ-TREE, RAXML) are recommended.
     * @throw std::invalid\_argument if the input alignment is empty or invalid
     */
    bool isEffective(const Alignment& aln);

    /*! @} End of Doxygen Groups*/

    /** \brief Run CMAPLE with user-specified parameters (via command-line)
     * @param[in] params user-specified parameters (via command-line)
     */
    void runCMAPLE(cmaple::Params& params);

    /** \brief Function for testing only
     * @param[in] params user-specified (parameters via command-line)
     */
    void testing(cmaple::Params& params);
}
