#!/bin/env python3
"""
This script shows that B.M.'s thesis results is flawed by a bug in the NumPy library, see https://github.com/numpy/numpy/issues/17277
Not only RLE but even Delta Encoding uses bugged loadtxt()

He used NumPy version 1.16.5

Create two different environment with two different version of NumPy.
This will show a lot of warnings about 1-counts that should not be there:
$ virtualenv venv-flawed
$ pip install  pip install numpy==1.16.5
$ . venv-flawed/bin/activate
$ python avg-rle-brunato.py

This is correct, showing no warning at all:
$ virtualenv venv-correct
$ pip install  pip install numpy==1.23.5
$ . venv-correct/bin/activate
$ python avg-rle-brunato.py
"""

# Quasi-copy-pasted Brunato's code below

# CELL 1
# relative path of the dataset
file_spss = "SRR001665_1.a2.unitigs.fa.ust.fa"
file_count = "SRR001665_1.a2.unitigs.fa.ust.counts"

# relative path of avg-permuted counts, RLE is + ".RLC_multi"
file_output = "SRR001665_2.unitigs.fa.ust.counts.permuted"

# CELL 2
# libraries
import numpy as np
from numpy import loadtxt
import time
from collections import Counter

# debug
print(f"numpy version: {np.__version__}")


# mini functions
def subroutine_avg(elem):
    return (sum(elem)) / len(elem)


def subroutine_top_frequent(elem):
    occurence_count = Counter(elem)
    return occurence_count.most_common(1)[0][0]


# CELL 3
# functions for Encoding********************************************************

def RLC_special_char_only_for_multi(abs_path, flag):  # output file: .RLC_multi
    print("avvio funzione RLC_special_char_only_for_multi")
    lines = loadtxt(abs_path, dtype=str, comments="#", delimiter=" ", unpack=False)

    # debug
    print(f"n_kmers(avg) = {len(lines)}")

    file_name = abs_path + ".RLC_multi"
    local_count = -1
    global_count = 0  # debug
    prova = open(file_name, 'w')
    prova.write(str(lines[0]))
    if ((lines[0] == lines[1])):
        prova.write('-')
        local_count = 1
        global_count = 1  # debug
    else:
        if (flag == 1):
            prova.write('\n')
        else:
            prova.write(' ')

    for i in range(1, lines.size, 1):
        if ((lines[i] == lines[i - 1])):
            local_count = local_count + 1
            global_count += 1  # debug
            if (local_count == 2 and i != 1):
                prova.write('-')
        else:
            if (local_count == 1):
                if (flag == 1):
                    prova.write('\n')
                else:
                    prova.write(' ')
            else:
                if (local_count != -1):  # if altrimenti STAMPO SEMpre il local_count del primo count
                    prova.write(str(local_count))
                if (i != 1):  # if i!=1 altrimenti STAMPO SEMpre il local_count del primo count
                    if (flag == 1):
                        prova.write('\n')
                    else:
                        prova.write(' ')
            prova.write(str(lines[i]))
            local_count = 1

            # debug
            global_count += 1
            if str(lines[i]) == "1":
                print(f"WARNING: line {i} has count 1!")

    prova.write(str(local_count))
    prova.close()

    # debug
    print(f"n_kmers(avg+rle) = {global_count}")

    print("terminata funzione RLC_special_char_only_for_multi")
    return


# Functions for Permutations***************************************

def sort_by_average(abs_path_spss, abs_path_count, abs_path_output, flag):  # output file: .avg
    print("Inizio funzione sort_by_average")
    size = len(open(abs_path_spss).readlines())
    size = int(size / 2)  # numero di unitig/nodi
    unitigs = open(abs_path_spss, "r")
    countsRAM = loadtxt(abs_path_count, dtype=int, comments="#", delimiter=" ", unpack=False)
    # creo una lista
    lista = []
    offset = 0
    m = 0
    for i in range(0, size, 1):
        unitig = unitigs.readline()  # metadato >
        unitig = unitigs.readline()  # unitigs max vero e proprio
        m = len(unitig) - 31
        lista.append(countsRAM[offset:offset + m])

        # debug
        if len(countsRAM[offset:offset + m]) == 0:
            print("WARNING: appending a 0 length count list!")

        offset += m
    lista.sort(key=subroutine_avg)
    # scrivo su disco
    name = abs_path_output + ".avg"
    prova = open(name, 'w')
    for i in range(0, len(lista), 1):
        for k in range(0, len(lista[i]), 1):
            prova.write(str(lista[i][k]))
            if (flag == 1):
                prova.write('\n')
            else:
                prova.write(' ')
    prova.close()
    print("Fine funzione sort_by_average")
    return


# CELL 4

sort_by_average(file_spss, file_count, file_output, 1)
RLC_special_char_only_for_multi(file_output + ".avg", 1)
