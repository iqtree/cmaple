#include "model.h"
#include "model_aa.h"
#include "model_dna.h"

using namespace std;
using namespace cmaple;

cmaple::Model::Model(const cmaple::ModelBase::SubModel sub_model, const cmaple::SeqRegion::SeqType seqtype):model_base(nullptr)
{
    cmaple::ModelBase::SubModel n_sub_model = sub_model;
    cmaple::SeqRegion::SeqType n_seqtype = seqtype;
    
    // Make sure either sub_model or n_seqtype is non-auto
    if (n_sub_model == ModelBase::DEFAULT && n_seqtype == SeqRegion::SEQ_AUTO)
        throw std::invalid_argument("Either sub_model or seqtype must be non-AUTO");
    
    // If sub_model is DEFAULT, select the default model according to the seqtype
    if (n_sub_model == ModelBase::DEFAULT || n_sub_model == ModelBase::UNKNOWN)
    {
        switch (n_seqtype) {
            case SeqRegion::SEQ_DNA:
            {
                if (cmaple::verbose_mode >= cmaple::VB_MED)
                    std::cout << "Auto select GTR model for DNA data" << std::endl;
                n_sub_model = ModelBase::GTR;
                break;
            }
            case SeqRegion::SEQ_PROTEIN:
            {
                if (cmaple::verbose_mode >= cmaple::VB_MED)
                    std::cout << "Auto select LG model for AA data" << std::endl;
                n_sub_model = ModelBase::LG;
                break;
            }
            default:
            {
                n_sub_model = ModelBase::UNKNOWN;
                throw std::invalid_argument("Unable to select a substitution model from an unknown/auto-detect seqtype");
                break;
            }
        }
    }
        
    // If sequence type is not specified -> detect it from sub_model
    if (n_seqtype == cmaple::SeqRegion::SEQ_AUTO || n_seqtype == cmaple::SeqRegion::SEQ_UNKNOWN)
        n_seqtype = ModelBase::detectSeqType(n_sub_model);
    
    // Init model from the corresponding seqtype and sub_model
    switch (n_seqtype) {
        case cmaple::SeqRegion::SEQ_PROTEIN:
        {
            model_base = new ModelAA(n_sub_model);
            break;
        }
        case cmaple::SeqRegion::SEQ_DNA:
        {
            model_base = new ModelDNA(n_sub_model);
            break;
        }
        default:
        {
            throw std::invalid_argument("Unknown/Unsupported model");
            break;
        }
    }
}

cmaple::Model::~Model()
{
    if (model_base)
    {
        delete model_base;
        model_base = nullptr;
    }
}

void cmaple::Model::fixParameters(const bool& n_fixed_model_params)
{
    if (model_base)
        model_base->fixed_params = n_fixed_model_params;
}

cmaple::Model::ModelParams cmaple::Model::getParams()
{
    // Init an empty ModelParams
    ModelParams model_params = {"", "", ""};
    
    // Handle cases when model is not yet initialized
    if (model_base)
    {
        // model_name
        model_params.model_name = model_base->getModelName();
        
        // root frequencies
        model_params.state_freqs = model_base->exportRootFrequenciesStr();
        
        // Q matrix
        model_params.mut_rates = model_base->exportQMatrixStr();
    }
    
    return model_params;;
}
