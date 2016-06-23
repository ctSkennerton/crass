#include <string>

#include "catch.hpp"
#include "libcrispr.h"
#include "ReadHolder.h"

// read
// 0                                                                                                   1
// 0         1         2         3         4         5         6         7         8         9         0
// 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
// CCCCGCAGGCGCGGGGATGAACCGAGCGAGACATCACCGGCGAGTCGGAGCGCGTTGCGTTCCCCGCAGGCGCGGGGATGAACCGAAGATAAACGCCGGCG
// RRRRRRRRRRRRRRRRRRRRRRRRR                                    RRRRRRRRRRRRRRRRRRRRRRRRR               
TEST_CASE("searching for additional repeated kmers in 100bp read", "[libcrispr]"){
    ReadHolder read("CCCCGCAGGCGCGGGGATGAACCGAGCGAGACATCACCGGCGAGTCGGAGCGCGTTGCGTTCCCCGCAGGCGCGGGGATGAACCGAAGATAAACGCCGGCG",
                    "HWI-EAS165_0052:2:58:8891:11288#CGATGT/1_C233_C23");
    SECTION("where there should be no additional matches with a minimum spacer length of 21"){
        read.startStopsAdd(0, 8);
        read.startStopsAdd(61,69);
        std::string pattern = "CCCCGCAGG";
        scanRight(read, pattern, 21, 24);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 8);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 69);
    }
    SECTION("where there should be no additional matches with a minimum spacer length of 10"){
        read.startStopsAdd(0, 8);
        read.startStopsAdd(61,69);
        std::string pattern = "CCCCGCAGG";
        scanRight(read, pattern, 10, 24);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 8);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 69);
    }
}
TEST_CASE("check extending repeat with 100bp read", "[libcrispr]") {
    ReadHolder read("CCCCGCAGGCGCGGGGATGAACCGAGCGAGACATCACCGGCGAGTCGGAGCGCGTTGCGTTCCCCGCAGGCGCGGGGATGAACCGAAGATAAACGCCGGCG",
                    "HWI-EAS165_0052:2:58:8891:11288#CGATGT/1_C233_C23");

    SECTION("The search window length is 8 and the min spacer length is 21") {
        read.startStopsAdd(0, 8);
        read.startStopsAdd(61,69);
        int repeat_length = extendPreRepeat(read, 8, 21);
        REQUIRE(repeat_length == 25);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 24);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 85);
    }
    SECTION("The search window length is 6 and the min spacer length is 21") {
        read.startStopsAdd(0, 6);
        read.startStopsAdd(61,67);
        int repeat_length = extendPreRepeat(read, 6, 21);
        REQUIRE(repeat_length == 25);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 24);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 85);
    }
    SECTION("The search window length is 11 and the min spacer length is 21") {
        read.startStopsAdd(0, 11);
        read.startStopsAdd(61,72);
        int repeat_length = extendPreRepeat(read, 11, 21);
        REQUIRE(repeat_length == 25);
        StartStopList reppos = read.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 24);
        REQUIRE(reppos[2] == 61);
        REQUIRE(reppos[3] == 85);
    }
}
// read2
// 0                                                                                                   1                                                 
// 0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9
// 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
// CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAAC
// RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR                                          RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR                                         RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
TEST_CASE("searching for additional repeated kmers in 190bp read", "[libcrispr]"){
    ReadHolder read2("CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAAC",
            "NC_019753_1");
    SECTION("where there should be one additional match with a minimum spacer length of 21"){
        read2.startStopsAdd(0, 8);
        read2.startStopsAdd(78,86);
        std::string pattern = "CTTCCACT";
        scanRight(read2, pattern, 21, 24);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 8);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 86);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 163);
    }
    SECTION("where there should be one additional match with a minimum spacer length of 10"){
        read2.startStopsAdd(0, 8);
        read2.startStopsAdd(78,86);
        std::string pattern = "CTTCCACT";
        scanRight(read2, pattern, 10, 24);
        StartStopList reppos = read2.getStartStopList();
        REQUIRE(reppos[0] == 0);
        REQUIRE(reppos[1] == 8);
        REQUIRE(reppos[2] == 78);
        REQUIRE(reppos[3] == 86);
        REQUIRE(reppos[4] == 155);
        REQUIRE(reppos[5] == 163);
    }
}
TEST_CASE("check extending repeat with 190bp read", "[libcrispr]") {

    ReadHolder read2("CTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACATTTTTTATCCAGATTTTTCCCCAAATTTGCAATAATTGCTACTTCCACTAACCATTTCCCCGTAAGGGGACGGAAACTTCTGTAGAGTTATTGTATAAGAACCCCACGTAGAAACGAGCTTCCACTAACCATTTCCCCGTAAGGGGACGGAAAC",
            "NC_019753_1");

    SECTION("The search window length is 8 and the min spacer length is 21") {
        read2.startStopsAdd(0,8);
        read2.startStopsAdd(78,86);
        read2.startStopsAdd(155,163);
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
        read2.startStopsAdd(0,6);
        read2.startStopsAdd(78,84);
        read2.startStopsAdd(155,161);
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
        read2.startStopsAdd(0,11);
        read2.startStopsAdd(78,89);
        read2.startStopsAdd(155,166);
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
        read2.startStopsAdd(0,8);
        read2.startStopsAdd(78,86);
        read2.startStopsAdd(155,163);
        read2.startStopsAdd(227,235);
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
        read2.startStopsAdd(0,6);
        read2.startStopsAdd(78,84);
        read2.startStopsAdd(155,161);
        read2.startStopsAdd(227,233);
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
        read2.startStopsAdd(0,11);
        read2.startStopsAdd(78,89);
        read2.startStopsAdd(155,166);
        read2.startStopsAdd(227,238);
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
