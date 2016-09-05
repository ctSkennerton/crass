#include <string>

#include "catch.hpp"
#include "libcrispr.h"
#include "ReadHolder.h"

// 0                                                                                                   1                         
// 0         1         2         3         4         5         6         7         8         9         0         1         2     
// 012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
// CACCATGGAAGACCTTCCTAACACCATGGTAGACATTCCTTACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTTCTAA
// rrrrrrrr                                                       rrrrrrrr                                  rrrrrrrr

TEST_CASE("searching for additional repeated kmer in a 126bp read", "[libcrispr]") {
    ReadHolder read("CACCATGGAAGACCTTCCTAACACCATGGTAGACATTCCTTACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTTCTAA","HWI-D00456:77:C70WLANXX:1:1101:10963:2182");
    SECTION("where there should be one additional match with a minimum spacer length of 26"){
        read.startStopsAdd(0, 7);
        read.startStopsAdd(63,70);
        std::string pattern = "CACCATGG";
        scanRight(read, pattern, 26, 24);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos.size() == 6);
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 63);
        REQUIRE(reppos[3] == 70);
        REQUIRE(reppos[4] == 105);
        REQUIRE(reppos[5] == 112);
    }
}

TEST_CASE("check extending repeat with 126bp read", "[libcrispr]") {
    ReadHolder read("CACCATGGAAGACCTTCCTAACACCATGGTAGACATTCCTTACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTCCTAACACCATGGTAGACCTTTCTAA","HWI-D00456:77:C70WLANXX:1:1101:10963:2182");

    SECTION("The search window length is 8 and the min spacer length is 26") {
        read.startStopsAdd(0, 7);
        read.startStopsAdd(63,70);
        read.startStopsAdd(105,112);
        int repeat_length = extendPreRepeat(read, 8, 26);
        REQUIRE(repeat_length == 23);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos.size() == 6);
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 21);
        REQUIRE(reppos[2] == 62);
        REQUIRE(reppos[3] == 84);
        REQUIRE(reppos[4] == 104);
        REQUIRE(reppos[5] == 125);
    }
}

TEST_CASE("searching for additional repeated kmers in 100bp read", "[libcrispr]"){
// read
// 0                                                                                                   1
// 0         1         2         3         4         5         6         7         8         9         0
// 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
// CCCCGCAGGCGCGGGGATGAACCGAGCGAGACATCACCGGCGAGTCGGAGCGCGTTGCGTTCCCCGCAGGCGCGGGGATGAACCGAAGATAAACGCCGGCG
// rrrrrrrrRRRRRRRRRRRRRRRRR                                    rrrrrrrrRRRRRRRRRRRRRRRRR               
//
    ReadHolder read("CCCCGCAGGCGCGGGGATGAACCGAGCGAGACATCACCGGCGAGTCGGAGCGCGTTGCGTTCCCCGCAGGCGCGGGGATGAACCGAAGATAAACGCCGGCG",
                    "HWI-EAS165_0052:2:58:8891:11288#CGATGT/1_C233_C23");
    SECTION("where there should be no additional matches with a minimum spacer length of 21"){
        read.startStopsAdd(0, 7);
        read.startStopsAdd(61,68);
        std::string pattern = "CCCCGCAG";
        scanRight(read, pattern, 21, 24);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos.size() == 4);
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 68);
    }
    SECTION("where there should be no additional matches with a minimum spacer length of 10"){
        read.startStopsAdd(0, 7);
        read.startStopsAdd(61,68);
        std::string pattern = "CCCCGCAG";
        scanRight(read, pattern, 10, 24);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos.size() == 4);
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 68);
    }
// read2
// 0                                                                                                   1
// 0         1         2         3         4         5         6         7         8         9         0
// 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
// TTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTT
// rrrrrrrrRRRRRRRRR                           RrrrrrrrrRRRRRRRRR            RrrrrrrrrRRRRRRRRR          Spacer Length: 21 
// rrrrrrrrRRRRRRR                                RrrrrrrrrRRRRRRR                    RrrrrrrrrRRRRRR    Spacer Length: 24                    
    ReadHolder read2("TTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTT", 
            "SRR438795.13216");

    SECTION("where there should be one additional match with a minimum spacer length of 21"){
        read2.startStopsAdd(0, 7);
        read2.startStopsAdd(45,52);
        std::string pattern = "TTGTTGTT";
        scanRight(read2, pattern, 21, 24);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos.size() == 6);
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 45);
        REQUIRE(reppos[3] == 52);
        REQUIRE(reppos[4] == 75);
        REQUIRE(reppos[5] == 82);
    }

    SECTION("where there should be one additional match with a minimum spacer length of 24"){
        //0,7,48,55,81,88,
        read2.startStopsAdd(0, 7);
        read2.startStopsAdd(48,55);
        std::string pattern = "TTGTTGTT";
        scanRight(read2, pattern, 24, 24);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos.size() == 6);
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 48);
        REQUIRE(reppos[3] == 55);
        REQUIRE(reppos[4] == 81);
        REQUIRE(reppos[5] == 88);
    }

    SECTION("where there should be three additional matches with a minimum spacer length of 10"){
        //0,7,
        //33,40,
        //51,58,
        //69,76,
        //87,94,
        read2.startStopsAdd(0, 7);
        read2.startStopsAdd(33,40);
        std::string pattern = "TTGTTGTT";
        scanRight(read2, pattern, 10, 24);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos.size() == 10);
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 33);
        REQUIRE(reppos[3] == 40);
        REQUIRE(reppos[4] == 51);
        REQUIRE(reppos[5] == 58);
        REQUIRE(reppos[6] == 69);
        REQUIRE(reppos[7] == 76);
        REQUIRE(reppos[8] == 87);
        REQUIRE(reppos[9] == 94);
    }
}
TEST_CASE("check extending repeat with 100bp read", "[libcrispr]") {
    ReadHolder read("CCCCGCAGGCGCGGGGATGAACCGAGCGAGACATCACCGGCGAGTCGGAGCGCGTTGCGTTCCCCGCAGGCGCGGGGATGAACCGAAGATAAACGCCGGCG",
                    "HWI-EAS165_0052:2:58:8891:11288#CGATGT/1_C233_C23");

    SECTION("The search window length is 8 and the min spacer length is 21") {
        read.startStopsAdd(0, 7);
        read.startStopsAdd(61,68);
        int repeat_length = extendPreRepeat(read, 8, 21);
        REQUIRE(repeat_length == 25);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 24);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 85);
    }
    SECTION("The search window length is 6 and the min spacer length is 21") {
        read.startStopsAdd(0, 5);
        read.startStopsAdd(61,66);
        int repeat_length = extendPreRepeat(read, 6, 21);
        REQUIRE(repeat_length == 25);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 24);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 85);
    }
    SECTION("The search window length is 11 and the min spacer length is 21") {
        read.startStopsAdd(0, 10);
        read.startStopsAdd(61,71);
        int repeat_length = extendPreRepeat(read, 11, 21);
        REQUIRE(repeat_length == 25);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 24);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 85);
    }

    ReadHolder read2("TTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTT", 
            "SRR438795.13216");
    SECTION("The search window length is 8 and the min spacer length is 21") {
        read2.startStopsAdd(0, 7);
        read2.startStopsAdd(45,52);
        read2.startStopsAdd(75,82);
        int repeat_length = extendPreRepeat(read2, 8, 21);
        REQUIRE(repeat_length == 18);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 16);
        REQUIRE(reppos[2] == 44);
        REQUIRE(reppos[3] == 61);
        REQUIRE(reppos[4] == 74);
        REQUIRE(reppos[5] == 91);
    }
    SECTION("The search window length is 8 and the min spacer length is 24") {
        //0,7,48,55,81,88,
        read2.startStopsAdd(0, 7);
        read2.startStopsAdd(48,55);
        read2.startStopsAdd(81,88);
        int repeat_length = extendPreRepeat(read2, 8, 24);
        REQUIRE(repeat_length == 18);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 16);
        REQUIRE(reppos[2] == 47);
        REQUIRE(reppos[3] == 64);
        REQUIRE(reppos[4] == 80);
        REQUIRE(reppos[5] == 97);
    }
}
// read2
// 0                                                                                                   1                                                 
// 0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9
// 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
// CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAAC
// rrrrrrrrRRRRRRRRRRRRRRRRRRRRRRRRRRRR                                          rrrrrrrrRRRRRRRRRRRRRRRRRRRRRRRRRRRR                                         rrrrrrrrRRRRRRRRRRRRRRRRRRRRRRRRRRRR
TEST_CASE("searching for additional repeated kmers in 190bp read", "[libcrispr]"){
    ReadHolder read2("CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAAC",
            "NC_019753_1");
    SECTION("where there should be one additional match with a minimum spacer length of 21"){
        read2.startStopsAdd(0, 7);
        read2.startStopsAdd(78,85);
        std::string pattern = "CTTCCACT";
        scanRight(read2, pattern, 21, 24);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 85);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 162);
    }
    SECTION("where there should be one additional match with a minimum spacer length of 10"){
        read2.startStopsAdd(0, 7);
        read2.startStopsAdd(78,85);
        std::string pattern = "CTTCCACT";
        scanRight(read2, pattern, 10, 24);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 7);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 85);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 162);
    }
}
TEST_CASE("check extending repeat with 190bp read", "[libcrispr]") {

    ReadHolder read2("CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAAC",
            "NC_019753_1");

    SECTION("The search window length is 8 and the min spacer length is 21") {
        read2.startStopsAdd(0,7);
        read2.startStopsAdd(78,85);
        read2.startStopsAdd(155,162);
        int repeat_length = extendPreRepeat(read2, 8, 21);
        REQUIRE(repeat_length == 36);
        StartStopList reppos = read2.getStartStopList();

        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 35);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 113);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 190);
    }
    SECTION("The search window length is 6 and the min spacer length is 21") {
        read2.startStopsAdd(0,5);
        read2.startStopsAdd(78,83);
        read2.startStopsAdd(155,160);
        int repeat_length = extendPreRepeat(read2, 6, 21);
        REQUIRE(repeat_length == 36);
        StartStopList reppos = read2.getStartStopList();

        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 35);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 113);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 190);
    }
    SECTION("The search window length is 11 and the min spacer length is 21") {
        read2.startStopsAdd(0,10);
        read2.startStopsAdd(78,88);
        read2.startStopsAdd(155,165);
        int repeat_length = extendPreRepeat(read2, 11, 21);
        REQUIRE(repeat_length == 36);
        StartStopList reppos = read2.getStartStopList();

        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 35);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 113);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 190);
    }
}
// read3
// 0                                                                                                   1                                                                                                   2                                                                                                   3
// 0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0     
// 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
// CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACCAGAAATTACGACAGACGCGCAAGCAGGATCAGCCACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACAGCCAGATAATGATGATAATGCTAGAGGCTGTCCGTT
// RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR                                          RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR                                         RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR                                    RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
TEST_CASE("check extending repeat with 300bp read", "[libcrispr]") {

    ReadHolder read2("CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACCAGAAATTACGACAGACGCGCAAGCAGGATCAGCCACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACAGCCAGATAATGATGATAATGCTAGAGGCTGTCCGTT",
            "NC_019753_1");

    SECTION("The search window length is 8 and the min spacer length is 21") {
        read2.startStopsAdd(0,7);
        read2.startStopsAdd(78,85);
        read2.startStopsAdd(155,162);
        read2.startStopsAdd(227,234);
        int repeat_length = extendPreRepeat(read2, 8, 21);
        REQUIRE(repeat_length == 36);
        StartStopList reppos = read2.getStartStopList();

        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 35);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 113);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 190);
        REQUIRE(reppos[6] == 227);
        REQUIRE(reppos[7] == 262);
    }
    SECTION("The search window length is 6 and the min spacer length is 21") {
        read2.startStopsAdd(0,5);
        read2.startStopsAdd(78,83);
        read2.startStopsAdd(155,160);
        read2.startStopsAdd(227,232);
        int repeat_length = extendPreRepeat(read2, 6, 21);
        REQUIRE(repeat_length == 36);
        StartStopList reppos = read2.getStartStopList();

        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 35);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 113);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 190);
        REQUIRE(reppos[6] == 227);
        REQUIRE(reppos[7] == 262);
    }
    SECTION("The search window length is 11 and the min spacer length is 21") {
        read2.startStopsAdd(0,10);
        read2.startStopsAdd(78,88);
        read2.startStopsAdd(155,165);
        read2.startStopsAdd(227,237);
        int repeat_length = extendPreRepeat(read2, 11, 21);
        REQUIRE(repeat_length == 36);
        StartStopList reppos = read2.getStartStopList();

        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 35);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 113);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 190);
        REQUIRE(reppos[6] == 227);
        REQUIRE(reppos[7] == 262);
    }
}

