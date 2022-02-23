#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as inp:
    with open('inverse_'+sys.argv[1],'w') as outp:
        for line in inp:
            contents = line.strip().split()
            delta = float(contents[3])
            inversed = delta*-1
            contents[3] = str(inversed)
            outp.write('\t'.join(contents)+'\n')



