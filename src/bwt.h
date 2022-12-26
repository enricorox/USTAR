//
// Created by enrico on 26/12/22.
//

#ifndef USTAR_BWT_H
#define USTAR_BWT_H

#include <cstdint>

typedef uint32_t byte;

void bwt_encode(byte *buf_in, byte *buf_out, int size, int *primary_index);

void bwt_decode(byte *buf_in, byte *buf_out, int size, int primary_index);


#endif //USTAR_BWT_H
