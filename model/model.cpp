#include "model.h"
#include "model_aa.h"
#include "model_dna.h"

using namespace std;
using namespace cmaple;

cmaple::Model::Model(const SubModel sub_model, const SeqType n_seqtype):model_base(nullptr)
{
    SeqType seqtype = n_seqtype;
    // If sequence type is not specified -> detect it from sub_model
    if (seqtype == SEQ_UNKNOWN)
        seqtype = ModelBase::detectSeqType(sub_model);
    
    // Init model from the corresponding seqtype and sub_model
    switch (seqtype) {
        case SEQ_PROTEIN:
        {
            model_base = new ModelAA(sub_model);
            break;
        }
        case SEQ_DNA:
        {
            model_base = new ModelDNA(sub_model);
            break;
        }
        default:
        {
            outError("Unsupported model");
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

std::map<std::string, std::string> cmaple::Model::getParams()
{
    // Handle cases when model is not yet initialized
    if (!model_base)
    {
        std::map<std::string, std::string> empty_map;
        return empty_map;
    }
    
    return model_base->exportModelParams();
}
