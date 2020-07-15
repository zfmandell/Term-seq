#!/usr/bin/env python3

"""
this program takes a rev strand .bedgraph file (can be for either Term-seq called 3' ends or RNA-seq coverage)
and returns invsersed .bedgraph file with coverage values all negative, used for IGV viewing
run as : rev_inverse.py (input file) (output file)
"""

import sys

with open(sys.argv[1],'r') as inp:
    with open(sys.argv[2],'w') as outp:
        for line in inp:
            contents = line.strip().split('\t')
            contents[-1] = str(float(contents[-1])*-1.0)
            outp.write('\t'.join(contents))
            outp.write('\n')
