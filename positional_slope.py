#!/usr/bin/env python3

"""
runs python3
requires numpy

This program is used to find the Cv at each base position

Cv is calculated as shown below
Cv=np.mean(coverage values in 15bp upstream)-np.mean(depth values in 15bp downstream)

run in a directory with your stranded coverage files (*_fwd.cov)/(*_rev.cov)
make sure only .cov files are those to be run

if you have a bam file, you can get a coverage file using
bedtools  genomecov -d  -ibam file.bam -g ref.fai >file.cov

output files are in bedgraph format placed in same directory as input files
Run as : python positional_slope.py <plus/minus> (+optional -window_size x)
"""

import argparse
import numpy as np
import glob

def read_coverage(wig_fyle):
    #read in a coverage file, returns coverage values as a flat list
    return_list = []
    with open(wig_fyle,"r") as inp:
        for line in inp:
            return_list.append(float(line.strip().split("\t")[2]))
    return return_list

def delta_calc(coverage_list,windows,strand):
    #creates list of delta value, position dependent, does NOT handle circular nature of prokaryotic chromasome
    # if window size = 15, first 15 nt, last 15 nt will be registered as 0
    #only keeps positive Cv (where Term-seq coverage decreases)
    output = []
    counter = 0
    length = len(coverage_list)
    while counter < length:
        if counter < windows:
            output.append(0.0)
            counter += 1
        elif counter >= windows and counter <= (length - windows):
            left = np.mean(coverage_list[counter-windows:counter+1],dtype=np.float64)
            right = np.mean(coverage_list[counter:(counter+windows+1)],dtype=np.float64)
            if strand.lower() == 'plus' or strand.lower() == 'both':
                to_add = left - right
                if to_add >= 0:
                    output.append(to_add)
                else:
                    output.append(0.0)
            elif strand.lower() == 'minus':
                to_add = right - left
                if to_add >= 0:
                    output.append(to_add)
                else:
                    output.append(0.0)
            counter += 1
        elif counter > (length - windows):
            output.append(0.0)
            counter += 1

    return output

def writer(wig_fyle,delta_list):
    #creates output Cv (delta) bedgraph file
    new_name = 'delta_'+str(wig_fyle[:-4])+".bedgraph"
    counter = 0
    with open(wig_fyle,'r') as inp:
        firstline = inp.readline()
        first_cont = firstline.strip().split("\t")
        sample = first_cont[0]
        with open(new_name,"w") as outp:
            for item in delta_list:
                outp.write(str(sample) + "\t"+ str(counter)+"\t"+str(counter+1)+"\t")
                outp.write(str(delta_list[counter]))
                outp.write("\n")
                counter += 1

def main():
    parser = argparse.ArgumentParser(description='calculates the Cv value at each position, outputs bedgraph files')
    parser.add_argument('strand',type=str,help='specify strand: <plus/minus>, no default value, must be inputted')
    parser.add_argument('-window_size',type=int,default=15,help='nucleotides upstream/downstream of x used to calculate delta, default = 15')
    parser.add_argument('-file', default = None, help='Specific .cov file to analyze')

    args = parser.parse_args()

    #for all coverage files in directory, read file into list, calculate delta values, create final bedgraph file
    if args.file != None:
        sub_data = read_coverage(args.file)
        sub_delta = delta_calc(sub_data,args.window_size,args.strand)
        writer(args.file,sub_delta)

    else:
        for fyle in sorted(glob.glob('*.cov')):
            sub_data = read_coverage(fyle)
            sub_delta = delta_calc(sub_data,args.window_size,args.strand)
            writer(fyle,sub_delta)

if __name__ == '__main__':
    main()
