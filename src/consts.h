//
// Created by enrico on 29/12/22.
//

#ifndef USTAR_CONSTS_H
#define USTAR_CONSTS_H

#include <map>

#define VERSION "0.1"

using namespace std;

enum class seeding_method_t{
    FIRST,
    LOWER_MEDIAN_ABUNDANCE,
    LOWER_AVERAGE_ABUNDANCE,
    SIMILAR_ABUNDANCE,
    BIGGER_LENGTH,
    SMALLER_LENGTH,
    MORE_CONNECTED,
    LESS_CONNECTED
};

enum class extending_method_t{
    FIRST,
    SIMILAR_ABUNDANCE,
    BIGGER_LENGTH,
    SMALLER_LENGTH,
    MORE_CONNECTED,
    LESS_CONNECTED
};

enum class encoding_t{
    PLAIN,
    RLE,
    AVG_RLE,
    FLIP,
    FLIP_RLE,
    AVG_FLIP_RLE,
    BINARY,
    BWT
};
// ----------------------------------------
const map<encoding_t, string> encoding_suffixes = {
        {encoding_t::PLAIN, ""},
        {encoding_t::RLE, ".rle"},
        {encoding_t::FLIP, ".flip"},
        {encoding_t::AVG_RLE, ".avg_rle"},
        {encoding_t::FLIP_RLE, ".flip_rle"},
        {encoding_t::AVG_FLIP_RLE, ".avg_flip_rle"},
        {encoding_t::BINARY, ".bin"},
        {encoding_t::BWT, ".bwt"}
};

const map<string, encoding_t> encoding_names = {
        {"plain", encoding_t::PLAIN},
        {"flip", encoding_t::FLIP},
        {"rle", encoding_t::RLE},
        {"avg_rle", encoding_t::AVG_RLE},
        {"flip_rle", encoding_t::FLIP_RLE},
        {"avg_flip_rle", encoding_t::AVG_FLIP_RLE},
        {"bin", encoding_t::BINARY},
        {"bwt", encoding_t::BWT}
};

const map<string, seeding_method_t> seeding_method_names = {
        {"f", seeding_method_t::FIRST},
        {"-ma", seeding_method_t::LOWER_MEDIAN_ABUNDANCE},
        {"-aa", seeding_method_t::LOWER_AVERAGE_ABUNDANCE},
        {"=a", seeding_method_t::SIMILAR_ABUNDANCE},
        {"-l", seeding_method_t::SMALLER_LENGTH},
        {"+l", seeding_method_t::BIGGER_LENGTH},
        {"-c", seeding_method_t::LESS_CONNECTED},
        {"+c", seeding_method_t::MORE_CONNECTED}
};

const map<string, extending_method_t> extending_method_names = {
        {"f", extending_method_t::FIRST},
        {"=a", extending_method_t::SIMILAR_ABUNDANCE},
        {"-l", extending_method_t::SMALLER_LENGTH},
        {"+l", extending_method_t::BIGGER_LENGTH},
        {"-c", extending_method_t::LESS_CONNECTED},
        {"+c", extending_method_t::MORE_CONNECTED}
};
// -----------------------------------------------------
#endif //USTAR_CONSTS_H
