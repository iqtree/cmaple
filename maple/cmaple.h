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
     * @param[in] max_subs_per_site the maximum number of substitution per sites
     * that CMAPLE is effective (optional)
     * @param[in] mean_subs_per_site the mean number of substitution per sites
     * that CMAPLE is effective (optional)
     * @return TRUE if the [C]Maple algorithm is effective to analyse the alignment;
     * otherwise, classical methods (e.g., IQ-TREE, RAXML) are recommended.
     * @throw std::invalid\_argument if the input alignment is empty or invalid
     */
    bool isEffective(const Alignment& aln,
            const double max_subs_per_site = MAX_SUBS_PER_SITE,
            const double mean_subs_per_site = MEAN_SUBS_PER_SITE);

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
