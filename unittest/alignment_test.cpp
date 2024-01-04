#include "gtest/gtest.h"
#include "../alignment/alignment.h"
using namespace cmaple;

/*
 Test read() and generateRef()
 */
TEST(Alignment, readSequences)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Alignment aln;
    
    // Test readFasta()
    aln.read(example_dir + "input.fa");
    EXPECT_EQ(aln.data.size(), 10);
    EXPECT_EQ(aln.data[9].seq_name, "T6"); // after sorting
    
    // Test generateRef()
    EXPECT_EQ(aln.ref_seq.size(), 20);
    EXPECT_EQ(aln.ref_seq[0], 0);
    EXPECT_EQ(aln.ref_seq[7], 2);
    EXPECT_EQ(aln.ref_seq[17], 0);
    
    // Test readPhylip()
    aln.read(example_dir + "input.phy");
    EXPECT_EQ(aln.data.size(), 10);
    EXPECT_EQ(aln.data[9].seq_name, "T6"); // after sorting
    
    // Test generateRef()
    EXPECT_EQ(aln.ref_seq.size(), 20);
    EXPECT_EQ(aln.ref_seq[0], 0);
    EXPECT_EQ(aln.ref_seq[10], 3);
    EXPECT_EQ(aln.ref_seq[18], 1);
    
    // Test generateRef() with empty input -> already moved to private section
    // EXPECT_EXIT(aln.generateRef(sequences), ::testing::ExitedWithCode(2), ".*");
    
    // Test input not found
    EXPECT_THROW(aln.read(example_dir + "notfound"), std::ios_base::failure);
}

/*
 Test readRefSeq(const std::string& ref_path)
 */
TEST(Alignment, readRef)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Alignment aln;
    aln.setSeqType(cmaple::SeqRegion::SEQ_DNA);
    std::string ref_seq = aln.readRefSeq(example_dir + "ref.fa", "REF");
    
    EXPECT_EQ(ref_seq.length(), 20);
    EXPECT_EQ(ref_seq[0], 'A');
    EXPECT_EQ(ref_seq[2], 'T');
    EXPECT_EQ(ref_seq[7], 'G');
    EXPECT_EQ(ref_seq[11], 'A');
    
    // Test readRef() with an empty seq_name
    EXPECT_THROW(aln.readRefSeq(example_dir + "ref.fa", ""), std::invalid_argument);
    
    //  Test readRef() from not-found file
    EXPECT_THROW(aln.readRefSeq(example_dir + "notfound", "REF"), std::ios_base::failure);
}

/*
 extractMutations() -> was moved to private
 Test extractMutations(StrVector &sequences, StrVector &seq_names,
 std::string ref_sequence, std::ofstream &out, bool only_extract_diff)
 Also test outputMutation()
 */
/*TEST(Alignment, extractMutations)
{
    // Don't need to test empty inputs (sequences, seq_names, ref_sequence);
    // we have ASSERT to check them in extractMutations()
    // ----- test on input.fa; generate ref_sequence -----
    Alignment aln;
    StrVector sequences, seq_names;
    std::string ref_sequence;
    
    // output file
    std::string diff_file_path(example_dir + "output.maple");
    
    // read sequences
    std::string file_path = example_dir + "input.fa";
    aln.readSequences(file_path.c_str(), sequences, seq_names);
    // detect the type of the input sequences
    aln.setSeqType(aln.detectSequenceType(sequences));
    
    // generate ref_sequence
    ref_sequence = aln.generateRef(sequences);

    // open the output file
    std::ofstream out = std::ofstream(diff_file_path);
    
    // write reference sequence to the output file
    out << ">" << REF_NAME << std::endl;
    out << ref_sequence << std::endl;
    
    // extract and write mutations of each sequence to file
    aln.extractMutations(sequences, seq_names, ref_sequence, out, true);
    out.close();
    
    // test the output data
    EXPECT_EQ(aln.data.size(), 0); // because only_extract_diff = true;
    EXPECT_EQ(aln.ref_seq.size(), 20);
    // ----- test on input.fa; generate ref_sequence-----
    
    // ----- test on input.fa; read ref_sequence -----
    std::string ref_file_path = example_dir + "ref.fa";
    ref_sequence = aln.readRef(ref_file_path);
    
    // open the output file
    out = std::ofstream(diff_file_path);
    
    // write reference sequence to the output file
    out << ">" << REF_NAME << std::endl;
    out << ref_sequence << std::endl;
    
    // extract and write mutations of each sequence to file
    aln.extractMutations(sequences, seq_names, ref_sequence, out, false);
    out.close();
    
    // test the output data
    EXPECT_EQ(aln.data.size(), 10);
    EXPECT_EQ(aln.data[1].seq_name, "T2");
    EXPECT_EQ(aln.data[6].size(), 2);
    EXPECT_EQ(aln.data[7][0].type, TYPE_DEL);
    EXPECT_EQ(aln.data[7][1].getLength(), 2);
    EXPECT_EQ(aln.data[9][0].position, 7);
    EXPECT_EQ(aln.ref_seq.size(), 20);
    // ----- test on input.fa; read ref_sequence -----
    
    // ----- Test to invalid data , violating assert(str_sequences.size() == seq_names.size()
    // && str_sequences.size() > 0 && out)
 #ifdef DEBUG
    // str_sequences.size() != seq_names.size()
    sequences.pop_back();
    EXPECT_DEATH(aln.extractMutations(sequences, seq_names, ref_sequence, out, false), ".*");
    // str_sequences.size() == 0
    sequences.clear();
    EXPECT_DEATH(aln.extractMutations(sequences, seq_names, ref_sequence, out, false), ".*");
 #endif
    
    // ----- test on input_full.phy; generate ref_sequence -----
    aln.data.clear();
    sequences.clear();
    seq_names.clear();
    ref_sequence = "";
    
    // read sequences
    file_path = example_dir + "input_full.phy";
    aln.readSequences(file_path.c_str(), sequences, seq_names);
    
    // generate ref_sequence
    ref_sequence = aln.generateRef(sequences);
    
    // open the output file
    out = std::ofstream(file_path);
    
    // write reference sequence to the output file
    out << ">" << REF_NAME << std::endl;
    out << ref_sequence << std::endl;
    
    // extract and write mutations of each sequence to file
    aln.extractMutations(sequences, seq_names, ref_sequence, out, false);
    out.close();
    
    // test the output data
    EXPECT_EQ(aln.data.size(), 1000);
    EXPECT_EQ(aln.data[3].seq_name, "T881");
    EXPECT_EQ(aln.data[220].size(), 284);
    EXPECT_EQ(aln.data[390][7].type, 2);
    EXPECT_EQ(aln.data[477][186].getLength(), 1);
    EXPECT_EQ(aln.data[796][356].position, 28897);
    EXPECT_EQ(aln.ref_seq.size(), 29903);
    // ----- test on input_full.phy; generate ref_sequence -----
}*/

/*
 Test readMapleFile(const std::string& diff_path, const std::string& ref_path)
 */
TEST(Alignment, readMapleFile)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Alignment aln;
    
    // ----- test on test_100.maple -----
    aln.read(example_dir + "test_100.maple");
    
    // test the output data
    EXPECT_EQ(aln.data.size(), 100);
    EXPECT_EQ(aln.data[4].seq_name, "EPI_ISL_530098");
    EXPECT_EQ(aln.data[18].size(), 8);
    EXPECT_EQ(aln.data[90][3].type, 1);
    EXPECT_EQ(aln.data[53][7].getLength(), 1);
    EXPECT_EQ(aln.data[46][7].position, 29871);
    EXPECT_EQ(aln.data[38].size(), 8);
    EXPECT_EQ(aln.data[49][5].type, 3);
    EXPECT_EQ(aln.data[66][10].getLength(), 1);
    EXPECT_EQ(aln.data[75][3].position, 3036);
    EXPECT_EQ(aln.data[30].size(), 8);
    EXPECT_EQ(aln.data[43][5].type, 1);
    EXPECT_EQ(aln.data[51][7].getLength(), 1);
    EXPECT_EQ(aln.data[66][4].position, 7797);

    EXPECT_EQ(aln.ref_seq.size(), 29891);
    EXPECT_EQ(aln.ref_seq[8], 3);
    EXPECT_EQ(aln.ref_seq[467], 0);
    EXPECT_EQ(aln.ref_seq[1593], 1);
    // ----- test on test_100.maple -----
    
    // ----- test on test_5K.maple, load ref_seq from test_100.maple -----
    //std::string ref_seq = aln.readRefSeq(example_dir + "ref_test_100.maple", "REF");
    aln.read(example_dir + "test_5K.maple");
    
    // test the output data
    EXPECT_EQ(aln.data.size(), 5000);
    EXPECT_EQ(aln.data[454].seq_name, "3521");
    EXPECT_EQ(aln.data[1328].size(), 12);
    EXPECT_EQ(aln.data[943][8].type, 3);
    EXPECT_EQ(aln.data[953][9].getLength(), 1);
    EXPECT_EQ(aln.data[76][5].position, 23402);
    EXPECT_EQ(aln.data[1543].size(), 13);
    EXPECT_EQ(aln.data[2435][8].type, 3);
    EXPECT_EQ(aln.data[4864][17].getLength(), 1);
    EXPECT_EQ(aln.data[3854][10].position, 23402);
    EXPECT_EQ(aln.data[2454].size(), 14);
    EXPECT_EQ(aln.data[4423][14].type, 0);
    EXPECT_EQ(aln.data[643][5].getLength(), 1);
    EXPECT_EQ(aln.data[59][2].position, 3036);

    EXPECT_EQ(aln.ref_seq.size(), 29891);
    EXPECT_EQ(aln.ref_seq[8], 3);
    EXPECT_EQ(aln.ref_seq[467], 0);
    EXPECT_EQ(aln.ref_seq[1593], 1);
    // ----- test on test_5K.maple, load ref_seq from test_100.maple -----
    
    // ----- Test read() with an empty input
    EXPECT_THROW(aln.read(""), std::invalid_argument);
    
    //  Test readMapleFile() from not-found file
    EXPECT_THROW(aln.read(example_dir + "notfound"), std::ios_base::failure);
    
    // ----- Test readMapleFile() with wrong format file
    EXPECT_THROW(aln.read(example_dir + "input.fa", "", cmaple::Alignment::IN_MAPLE), std::invalid_argument);
}

/*
 Test write()
 */
TEST(Alignment, write)
{
    // detect the path to the example directory
    std::string example_dir = "../../example/";
    if (!fileExists(example_dir + "example.maple"))
        example_dir = "../example/";
    
    Alignment aln;
    
    // extract MAPLE file
    std::string ref_seq = aln.readRefSeq(example_dir + "ref.fa", "REF");
    aln.read(example_dir + "input.phy", ref_seq);
    aln.write(example_dir + "input.phy.maple", cmaple::Alignment::IN_MAPLE, true);
    
    // reset data
    aln.data.clear();
    aln.ref_seq.clear();
    
    // read MAPLE file (for testing)
    aln.read(example_dir + "input.phy.maple");
    
    // test the output data
    EXPECT_EQ(aln.data.size(), 10);
    EXPECT_EQ(aln.data[0].seq_name, "T3");
    EXPECT_EQ(aln.data[1].size(), 2);
    EXPECT_EQ(aln.data[2][0].type, 10);
    EXPECT_EQ(aln.data[3][0].getLength(), 1);
    EXPECT_EQ(aln.data[4][1].position, 17);
    EXPECT_EQ(aln.data[5].size(), 3);
    EXPECT_EQ(aln.data[6][1].type, TYPE_DEL);
    EXPECT_EQ(aln.data[7][2].getLength(), 1);
    EXPECT_EQ(aln.data[8][1].position, 12);
    EXPECT_EQ(aln.data[9].size(), 5);
    EXPECT_EQ(aln.data[0][0].type, 10);
    EXPECT_EQ(aln.data[1][1].getLength(), 1);
    EXPECT_EQ(aln.data[2][1].position, 19);

    EXPECT_EQ(aln.ref_seq.size(), 20);
    EXPECT_EQ(aln.ref_seq[8], 3);
    EXPECT_EQ(aln.ref_seq[2], 3);
    EXPECT_EQ(aln.ref_seq[15], 1);
    // ----- test on input.phy with ref file from ref.fa -----
    
    // ----- Test write() with an empty input
    EXPECT_THROW(aln.write(""), std::invalid_argument);
    
    /*// ----- test on input.fa without ref file, specifying MAPLE file path -----
    aln.data.clear();
    aln.ref_seq.clear();
    aln_file_path = example_dir + "input_full.fa";
    params.aln_path = aln_file_path;
    params.ref_path = "";
    
    diff_file_path = aln_file_path + ".output.maple";
    params.maple_path = diff_file_path;
    
    // extract MAPLE file
    aln.extractDiffFile(params);
    
    // reset data
    aln.data.clear();
    aln.ref_seq.clear();
    
    // read MAPLE file (for testing)
    aln.readMapleFile(diff_file_path, ""); // read ref_seq from the MAPLE file
    
    // test the output data
    EXPECT_EQ(aln.data.size(), 1000);
    EXPECT_EQ(aln.data[432].seq_name, "T400");
    EXPECT_EQ(aln.data[745].size(), 127);
    EXPECT_EQ(aln.data[23][12].type, 1);
    EXPECT_EQ(aln.data[8][6].getLength(), 1);
    EXPECT_EQ(aln.data[475][18].position, 4082);
    EXPECT_EQ(aln.data[51].size(), 355);
    EXPECT_EQ(aln.data[16][5].type, 1);
    EXPECT_EQ(aln.data[678][22].getLength(), 1);
    EXPECT_EQ(aln.data[87][14].position, 3216);
    EXPECT_EQ(aln.data[954].size(), 232);
    EXPECT_EQ(aln.data[245][7].type, 2);
    EXPECT_EQ(aln.data[58][4].getLength(), 1);
    EXPECT_EQ(aln.data[47][20].position, 2716);

    EXPECT_EQ(aln.ref_seq.size(), 29903);
    EXPECT_EQ(aln.ref_seq[8424], 0);
    EXPECT_EQ(aln.ref_seq[223], 3);
    EXPECT_EQ(aln.ref_seq[175], 1);
    // ----- with/without MAPLE file path -----*/
}

/*
 Test convertState2Char(StateType state)
 */
TEST(Alignment, convertState2Char)
{
    EXPECT_EQ(Alignment::convertState2Char(-1, cmaple::SeqRegion::SEQ_DNA), '?');
    EXPECT_EQ(Alignment::convertState2Char(0.1, cmaple::SeqRegion::SEQ_DNA), 'A');
    EXPECT_EQ(Alignment::convertState2Char(0.99, cmaple::SeqRegion::SEQ_DNA), 'A');
    EXPECT_EQ(Alignment::convertState2Char(1000, cmaple::SeqRegion::SEQ_DNA), '?');
    
    EXPECT_EQ(Alignment::convertState2Char(TYPE_N, cmaple::SeqRegion::SEQ_DNA), '-');
    EXPECT_EQ(Alignment::convertState2Char(TYPE_DEL, cmaple::SeqRegion::SEQ_DNA), '-');
    EXPECT_EQ(Alignment::convertState2Char(TYPE_INVALID, cmaple::SeqRegion::SEQ_DNA), '?');
    
    EXPECT_EQ(Alignment::convertState2Char(0, cmaple::SeqRegion::SEQ_DNA), 'A');
    EXPECT_EQ(Alignment::convertState2Char(3, cmaple::SeqRegion::SEQ_DNA), 'T');
    
    EXPECT_EQ(Alignment::convertState2Char(1+4+3, cmaple::SeqRegion::SEQ_DNA), 'R');
    EXPECT_EQ(Alignment::convertState2Char(2+8+3, cmaple::SeqRegion::SEQ_DNA), 'Y');
    EXPECT_EQ(Alignment::convertState2Char(2+4+8+3, cmaple::SeqRegion::SEQ_DNA), 'B');
    EXPECT_EQ(Alignment::convertState2Char(1+2+4+3, cmaple::SeqRegion::SEQ_DNA), 'V');
    EXPECT_EQ(Alignment::convertState2Char(1+2+3, cmaple::SeqRegion::SEQ_DNA), 'M');
}

/*
 Test convertChar2State(char state) -> moved to private
 */
/*TEST(Alignment, convertChar2State)
{
    Alignment aln;
    aln.setSeqType(SEQ_DNA);
    
    // ----- Test convertChar2State() with an invalid state
 #ifdef DEBUG
    EXPECT_DEATH(aln.convertChar2State('e'), ".*");
 #endif
    
    // convertChar2State requires input is a capital character
    
    EXPECT_EQ(aln.convertChar2State('-'), TYPE_DEL);
    EXPECT_EQ(aln.convertChar2State('?'), TYPE_N);
    EXPECT_EQ(aln.convertChar2State('.'), TYPE_N);
    EXPECT_EQ(aln.convertChar2State('~'), TYPE_N);
    
    EXPECT_EQ(aln.convertChar2State('C'), 1);
    EXPECT_EQ(aln.convertChar2State('G'), 2);
    EXPECT_EQ(aln.convertChar2State('X'), TYPE_N);
    EXPECT_EQ(aln.convertChar2State('Y'), 2+8+3);
    EXPECT_EQ(aln.convertChar2State('K'), 4+8+3);
    EXPECT_EQ(aln.convertChar2State('B'), 2+4+8+3);
    EXPECT_EQ(aln.convertChar2State('D'), 1+4+8+3);
    EXPECT_EQ(aln.convertChar2State('V'), 1+2+4+3);
}*/
