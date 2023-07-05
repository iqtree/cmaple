#include "model.h"
#include "model_aa.h"
#include "model_dna.h"

using namespace std;
using namespace cmaple;

cmaple::Model::Model(const std::string& model_name):model_base(nullptr)
{
    // Detect seqtype from model_name
    SeqType seqtype = ModelBase::detectSeqType(model_name);
    
    // Init model from the corresponding seqtype and model_name
    switch (seqtype) {
        case SEQ_PROTEIN:
        {
            model_base = new ModelAA(model_name);
            break;
        }
        case SEQ_DNA:
        {
            model_base = new ModelDNA(model_name);
            break;
        }
        default: // DNA
        {
            outError("Unsupported model " + model_name);
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
