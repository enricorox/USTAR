//
// Created by enrico on 22/12/22.
//

#ifndef USTAR_ENCODER_H
#define USTAR_ENCODER_H

#include "SPSS.h"

class Encoder{
public:
    explicit Encoder(SPSS *spss);

    void to_fasta_file(const string &file_name);

    void to_counts_file(const string &file_name);
};
#endif //USTAR_ENCODER_H
