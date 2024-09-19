#include "../utils/tools.h"

#pragma once

namespace cmaple {
    /** An alternative branch (when compute SPRTA)
     */
    class AltBranch {
    public:
        
        // Constructor
        AltBranch(RealNumType lh_val, cmaple::Index branch_id_val)
            : lh(lh_val), branch_id(branch_id_val) {}
    
        // Likelihood contribution
        RealNumType lh;
    
        // Branch id
        cmaple::Index branch_id;
    };
}
