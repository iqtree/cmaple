#include "gtest/gtest.h"
#include "../alignment/seqregions.h"
#include "../model/model_dna.h"
#include "../model/model_dna_rate_variation.h"
#include "../tree/tree.h"

using namespace cmaple;

cmaple::Alignment loadAln5K();

TEST(RateVariation, initMutationMat)
{
    Alignment aln = loadAln5K();
    Model model(aln.ref_seq.size(), true, false, cmaple::ModelBase::JC);
    Tree tree(&aln, &model);

    ModelDNARateVariation* model_JC = (ModelDNARateVariation*) tree.model;
    model_JC->setAllMatricesToDefault();
    
    EXPECT_EQ(model_JC->sub_model, cmaple::ModelBase::JC);
    EXPECT_EQ(model_JC->getRootFreq(1), 0.25);
    EXPECT_EQ(model_JC->getRootFreq(3), 0.25);
    EXPECT_EQ(model_JC->getInverseRootFreq(0), 4);
    EXPECT_EQ(model_JC->getInverseRootFreq(2), 4);
    EXPECT_EQ(model_JC->getRootLogFreq(0), log(0.25));
    EXPECT_EQ(model_JC->getRootLogFreq(3), model_JC->getRootLogFreq(0));

    for(int i = 0; i < aln.ref_seq.size(); ++i) {
        EXPECT_EQ(model_JC->getDiagonalMutationMatrixEntry(1,i), -1);
        EXPECT_EQ(model_JC->getDiagonalMutationMatrixEntry(2,i), -1);
        EXPECT_EQ(model_JC->getMutationMatrix(i)[5], -1);
        EXPECT_EQ(model_JC->getMutationMatrixEntry(1,1,i), -1);
        EXPECT_EQ(model_JC->getMutationMatrixRow(1,i)[1], -1);
        EXPECT_EQ(model_JC->getMutationMatrix(i)[10], -1);
        EXPECT_EQ(model_JC->getMutationMatrixEntry(2,2,i), -1);
        EXPECT_EQ(model_JC->getMutationMatrixRow(2,i)[2], -1);
        EXPECT_EQ(model_JC->getTransposedMutationMatrix(i)[0], -1);
        EXPECT_EQ(model_JC->getTransposedMutationMatrixEntry(0,0,i), -1);
        EXPECT_EQ(model_JC->getTransposedMutationMatrixRow(0,i)[0], -1);
        EXPECT_EQ(model_JC->getTransposedMutationMatrix(i)[15], -1);
        EXPECT_EQ(model_JC->getTransposedMutationMatrixEntry(3,3,i), -1);
        EXPECT_EQ(model_JC->getTransposedMutationMatrixRow(3,i)[3], -1);
        EXPECT_EQ(model_JC->getFreqiFreqjQij(0, 5, i), -1);
        EXPECT_EQ(model_JC->getFreqiFreqjQij(0, 10, i), -1);
        EXPECT_EQ(model_JC->getFreqjTransposedijRow(0, i)[0], -0.25);
        EXPECT_EQ(model_JC->getFreqjTransposedijRow(3, i)[3], -0.25);
    }
}

TEST(RateVariation, testMatrices) 
{
    Alignment aln = loadAln5K();
    Model model(aln.ref_seq.size(), true, false, cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;

    ModelDNARateVariation* rv_model = (ModelDNARateVariation*) tree.model;
    rv_model->setAllMatricesToDefault();

    for(int i = 0; i < aln.ref_seq.size(); ++i) {
        for(StateType a = 0; a < 4; ++a) {
            for(StateType b = 0; b < 4; ++b) {
                EXPECT_EQ(rv_model->getOriginalRateMatrix()[4 * a + b], rv_model->getMutationMatrixEntry(a,b,i));
                EXPECT_EQ(rv_model->getOriginalRateMatrix()[4 * a + b], rv_model->getTransposedMutationMatrixEntry(b,a,i));
                if(a == b) {
                    EXPECT_EQ(rv_model->getOriginalRateMatrix()[4 * a + b], rv_model->getDiagonalMutationMatrixEntry(a, i));
                }
            }
        }
    }
}

TEST(RateVariation, merge_N_O)
{
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Alignment aln = loadAln5K();
    Model model(aln.ref_seq.size(), true, false, cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);

    ModelDNARateVariation* rv_model = (ModelDNARateVariation*) tree.model;
    rv_model->setAllMatricesToDefault();

    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    const PositionType end_pos = 8623;
    RealNumType lower_plength = -1;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,
        8.3817209383128356075479820086471249851456377655268e-06,
        0.72401575181023847260775028189527802169322967529297,
        1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.3675361745020691,0.0000068532680661245309,
        0.63243117236202639,0.000025799867838435163};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge[i]);
    }
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    seqregion1.plength_observation2root = 0.311; // plength_observation2root
    // does NOT affect the likelihood of new merged region
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 2);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge[i]);
    }
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    seqregion1.plength_observation2node = 0; // plength_observation2node = 0
    // does NOT affect the likelihood of new merged region
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 3);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge[i]);
    }
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    seqregion1.plength_observation2node = 0.1252; // plength_observation2node > 0
    // affects the likelihood of new merged region
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge1{0.37599425541307963,0.0099624992407331657,
        0.55990605200203047,0.054137193344156891};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge1[i]);
    }
    // ----- Test 4 -----
    
    // ----- Test 5 -----
    lower_plength = 0; // lower_plength = 0; plength_observation2node > 0
    merge_N_O<4>(lower_plength, seqregion1, tree.model,
              end_pos, merged_regions);
    EXPECT_EQ(merged_regions.size(), 5);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge1[i]);
    }
}

TEST(RateVariation, merge_O_ORACGT)
{
    Alignment aln = loadAln5K();
    Model model(aln.ref_seq.size(), true, false, cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);

    ModelDNARateVariation* rv_model = (ModelDNARateVariation*) tree.model;
    rv_model->setAllMatricesToDefault();

    std::unique_ptr<Params>& params = tree.params;
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    
    // ----- Test 1 -----
    SeqRegions merged_regions;
    auto new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value{0.27595781923293211113090706021466758102178573608398,
        8.3817209383128356075479820086471249851456377655268e-06,
        0.72401575181023847260775028189527802169322967529297,
        1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value;
    SeqRegion seqregion1(TYPE_O, 243, -1, -1, std::move(new_lh));
    new_lh = cmaple::make_unique<SeqRegion::LHType>();
    SeqRegion::LHType new_lh_value1{0.27595781923293211113090706021466758102178573608398,
        8.3817209383128356075479820086471249851456377655268e-06,
        0.72401575181023847260775028189527802169322967529297,
        1.8047235891267821277644811672757896303664892911911e-05};
    (*new_lh) = new_lh_value1;
    SeqRegion seqregion2(TYPE_O, 243, -1, 1e-3, std::move(new_lh));
    RealNumType total_blength_1 = -1;
    RealNumType total_blength_2 = -1;
    const RealNumType threshold_prob = params->threshold_prob;
    const PositionType end_pos = 3213;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_EQ(merged_regions.size(), 1);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);
    SeqRegion::LHType new_lh_value_merge{0.12684687976595482,1.1702018350525203E-10,
        0.87315311957450514,5.425200212298543E-10};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge[i]);
    }
    // ----- Test 1 -----
    
    // ----- Test 2 -----
    total_blength_1 = -1;
    total_blength_2 = 1e-5;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_TRUE(merged_regions.size() == 2);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    /*EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_EQ(merged_regions.back().plength_observation2node, 0);*/
    SeqRegion::LHType new_lh_value_merge1{0.12684809781766063,1.3059895886789165E-10,
        0.87315190141833243,6.3340794416577515E-10};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge1[i]);
    }
    // ----- Test 2 -----
    
    // ----- Test 3 -----
    total_blength_1 = 0.01;
    total_blength_2 = -1;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    EXPECT_FALSE(merged_regions.size() > 3);
    EXPECT_FALSE(merged_regions.back().type != TYPE_O);
    // EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_FALSE(merged_regions.back().plength_observation2node != 0);
    SeqRegion::LHType new_lh_value_merge2{0.12842468075362262,1.1709109654051516E-8,
        0.87157516057518813,1.4696207956334386E-7};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge2[i]);
    }
    // ----- Test 3 -----
    
    // ----- Test 4 -----
    total_blength_1 = 0.011321;
    total_blength_2 = 1e-3;
    merge_O_ORACGT<4>(seqregion1, seqregion2, total_blength_1, total_blength_2, end_pos,
                      threshold_prob, tree.model, tree.aln, merged_regions);
    // EXPECT_EQ(merged_regions.size(), 4);
    EXPECT_EQ(merged_regions.back().type, TYPE_O);
    // EXPECT_EQ(merged_regions.back().plength_observation2root, 0);
    EXPECT_TRUE(merged_regions.back().plength_observation2node == 0);
    SeqRegion::LHType new_lh_value_merge3{0.12875794255340739,1.6716816801441426E-7,
        0.87123893271963804,0.0000029575587865296985};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*merged_regions.back().likelihood)[i], new_lh_value_merge3[i]);
    }
}



TEST(RateVariation, computeAbsoluteLhAtRoot)
{
    Alignment aln = loadAln5K();
    Model model(aln.ref_seq.size(), true, false, cmaple::ModelBase::GTR);
    std::unique_ptr<Params> params = ParamsBuilder().build();
    Tree tree(&aln, &model);

    ModelDNARateVariation* rv_model = (ModelDNARateVariation*) tree.model;
    rv_model->setAllMatricesToDefault();
    
    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    std::unique_ptr<SeqRegions> seqregions3 = nullptr;
    
    // dummy variables
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    std::unique_ptr<SeqRegions> seqregions_1 = aln.data[0]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_2 = aln.data[10]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_3 = aln.data[100]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    
    seqregions_1->mergeTwoLowers<4>(seqregions1, 1e-5, *seqregions_2, 123e-3, \
            tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_2->mergeTwoLowers<4>(seqregions2, 214e-4, *seqregions_3, 13e-8,
            tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_3->mergeTwoLowers<4>(seqregions3, 5421e-8, *seqregions_1, 1073e-6,
            tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    // ----- Test 1 on a more complex seqregions -----
    EXPECT_DOUBLE_EQ(seqregions_1->computeAbsoluteLhAtRoot<4>(tree.model,
                tree.cumulative_base), -40549.6785849070511176250874996185302734375);
    // ----- Test 1 on a more complex seqregions -----
    
    // ----- Test 2 on a more complex seqregions -----
    EXPECT_DOUBLE_EQ(seqregions_2->computeAbsoluteLhAtRoot<4>(tree.model,
                tree.cumulative_base), -40549.19068355686613358557224273681640625);
    // ----- Test 2 on a more complex seqregions -----
    
    // ----- Test 3 on a more complex seqregions -----
    EXPECT_DOUBLE_EQ(seqregions_3->computeAbsoluteLhAtRoot<4>(tree.model,
                tree.cumulative_base), -40548.41030627492000348865985870361328125);
    // ----- Test 3 on a more complex seqregions -----
}

TEST(RateVariation, computeTotalLhAtRoot)
{
    Alignment aln = loadAln5K();
    Model model(aln.ref_seq.size(), true, false, cmaple::ModelBase::GTR);
    Tree tree(&aln, &model);
    std::unique_ptr<Params>& params = tree.params;

    ModelDNARateVariation* rv_model = (ModelDNARateVariation*) tree.model;
    rv_model->setAllMatricesToDefault();

    std::unique_ptr<SeqRegions> seqregions1 = nullptr;
    std::unique_ptr<SeqRegions> seqregions2 = nullptr;
    std::unique_ptr<SeqRegions> seqregions3 = nullptr;
    std::unique_ptr<SeqRegions> seqregions_total_lh = nullptr;
    
    // Generate complex seqregions
    const PositionType seq_length = aln.ref_seq.size();
    const StateType num_states = aln.num_states;
    std::unique_ptr<SeqRegions> seqregions_1 = aln.data[20]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_2 = aln.data[200]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    std::unique_ptr<SeqRegions> seqregions_3 = aln.data[2000]
        .getLowerLhVector(aln.ref_seq.size(), aln.num_states, aln.getSeqType());
    
    seqregions_1->mergeTwoLowers<4>(seqregions1, 1073e-6, *seqregions_2, 13e-8,
                    tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_2->mergeTwoLowers<4>(seqregions2, 123e-3, *seqregions_3, 5421e-8,
                    tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    seqregions_3->mergeTwoLowers<4>(seqregions3, 214e-4, *seqregions_1, 1e-5,
                tree.aln, tree.model, tree.cumulative_rate, params->threshold_prob);
    
    // ----- Test 1 on a more complex seqregions -----
    seqregions1->computeTotalLhAtRoot<4>(seqregions_total_lh, tree.model);
    EXPECT_EQ(seqregions_total_lh->size(), 22);
    SeqRegion::LHType lh_value1_10{7.2185298165439574223086912192976286051226963991212e-10,
        0.00012093455079658274576269449962495627914904616773129,
        5.7783014021321174033127313583984435707563420692168e-09,
        0.99987905894904904879894047553534619510173797607422};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*seqregions_total_lh->at(10).likelihood)[i], lh_value1_10[i]);
    }
    SeqRegion::LHType lh_value1_12{0.00012109700727777851427414274043670161518093664199114,
        3.7920254781824222995438774532709486075887639344728e-09,
        0.99987887895812166405562493309844285249710083007812,
        2.0242575103081267528512744062821337998059334495338e-08};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*seqregions_total_lh->at(12).likelihood)[i], lh_value1_12[i]);
    }
    SeqRegion::LHType lh_value1_14{7.2185298165439574223086912192976286051226963991212e-10,
        0.00012093455079658274576269449962495627914904616773129,
        5.7783014021321174033127313583984435707563420692168e-09,
        0.99987905894904904879894047553534619510173797607422};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*seqregions_total_lh->at(14).likelihood)[i], lh_value1_14[i]);
    }
    SeqRegion::LHType lh_value1_20{0.999878969078505264178602374158799648284912109375,
        3.79202547818242147236326490024327373618007186451e-09,
        0.00012100688689399335864343987267943703045602887868881,
        2.0242575103081260911067843638599939026789797935635e-08};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*seqregions_total_lh->at(20).likelihood)[i], lh_value1_20[i]);
    }
    // ----- Test 1 on a more complex seqregions -----
    
    // ----- Test 2 on a more complex seqregions -----
    seqregions2->computeTotalLhAtRoot<4>(seqregions_total_lh, tree.model, 1e-5);
    EXPECT_EQ(seqregions_total_lh->size(), 26);
    SeqRegion::LHType lh_value2_17{1.0426241904738277496250383261089389463904808508232e-06,
        0.00036249067899242780116039752691392550332238897681236,
        6.3013501138813336668181852573411561024840921163559e-06,
        0.99963016534670312562838034864398650825023651123047};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*seqregions_total_lh->at(17).likelihood)[i], lh_value2_17[i]);
    }
    SeqRegion::LHType lh_value2_21{0.00042526525345440798929128045635650323674781247973442,
        2.4910168754262691925165963680033343052855343557894e-06,
        0.99955743551688391868026428710436448454856872558594,
        1.4808212786277890710097057680449950112233636900783e-05};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*seqregions_total_lh->at(21).likelihood)[i], lh_value2_21[i]);
    }
    // ----- Test 2 on a more complex seqregions -----
    
    // ----- Test 3 on a more complex seqregions -----
    seqregions3->computeTotalLhAtRoot<4>(seqregions_total_lh, tree.model, 0.03762);
    EXPECT_EQ(seqregions_total_lh->size(), 27);
    SeqRegion::LHType lh_value3_10{0.0036579520570985094192473230378936932538636028766632,
        0.93983378575187570547200266446452587842941284179688,
        0.0036637083625175163349718676641941783600486814975739,
        0.052844553828508410153741436943164444528520107269287};
    for(int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ((*seqregions_total_lh->at(10).likelihood)[i], lh_value3_10[i]);
    }
}
