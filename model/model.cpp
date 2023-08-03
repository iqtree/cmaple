#include "model.h"
#include "model_aa.h"
#include "model_dna.h"

using namespace std;
using namespace cmaple;

cmaple::Model::Model(const cmaple::ModelBase::SubModel sub_model, const cmaple::SeqRegion::SeqType n_seqtype):model_base(nullptr)
{
    cmaple::SeqRegion::SeqType seqtype = n_seqtype;
    // If sequence type is not specified -> detect it from sub_model
    if (seqtype == cmaple::SeqRegion::SEQ_UNKNOWN)
        seqtype = ModelBase::detectSeqType(sub_model);
    
    // Init model from the corresponding seqtype and sub_model
    switch (seqtype) {
        case cmaple::SeqRegion::SEQ_PROTEIN:
        {
            model_base = new ModelAA(sub_model);
            break;
        }
        case cmaple::SeqRegion::SEQ_DNA:
        {
            model_base = new ModelDNA(sub_model);
            break;
        }
        default:
        {
            throw std::invalid_argument("Unsupported model");
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
