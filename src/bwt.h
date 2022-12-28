//
// Created by enrico on 26/12/22.
//

#ifndef USTAR_BWT_H
#define USTAR_BWT_H

#include <cstdint>

void bwt_encode(uint32_t *buf_in, uint32_t *buf_out, unsigned long size, int *primary_index);

void bwt_decode(uint32_t *buf_in, uint32_t *buf_out, int size, int primary_index);


#endif //USTAR_BWT_H
