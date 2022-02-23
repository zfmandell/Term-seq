#!/usr/bin/env python3

"""
runs python3
requires BioPython
requires numpy
requires SciPy >= 0.11

This program finds all local maxima Cv values from genome-wide Cv files via SciPy

run in a directory with your stranded Cv (delta) files (delta_*_fwd.bedgraph)/(delta_*_rev.bedgraph)
make sure only .bedgraph files are those to be run

if you have a folder of cov files, you can get a Cv (delta) file using:
python positional_slope.py <plus/minus> (+optional -window_size x)


optionally returns a .fasta file, with sequences x nucleotides upstream of each peak, default length = 50
returns bedgraph file containing all 3' ends that pass threshold for loading onto IGV
Run as : python local_peak.py <plus/minus> (+optional -min_delta x) (+optional -fasta_out <true/false>) (+optional -fasta <.fasta file used for alignment>) (+optional -fasta_window x)


"""

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from scipy.signal import argrelextrema
import operator

def read_coverage(wig_fyle):
    #create list of coverage values
    return_list = []
    with open(wig_fyle,"r") as inp:
        for line in inp:
            return_list.append(float(line.strip().split("\t")[3]))
    return return_list

def read_position(wig_fyle):
    #create list of position values, 1:1 with coverage list
    return_list = []
    with open(wig_fyle,"r") as inp:
        for line in inp:
            return_list.append(int(line.strip().split("\t")[2]))
    return return_list

def read_fasta(genome_fasta):
    #reads in fasta
    record = SeqIO.read(genome_fasta, "fasta")
    return str(record.seq)

def reverse_complement(seq):
    #returns reverse complement of a list of nucleotides
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def local_max(delta_list,positions):
    #finds all local maxima, ranks in order of peak height,thresholds based on min_delta, smoothing optional via savitzky golay
    peak_dict = {}

    peaks = argrelextrema(np.asarray(delta_list), np.greater)[0].tolist()

    for item in peaks:
        peak_dict[positions[item]] = delta_list[item]
    return peak_dict

def writer(genome,peak_dict,fyle,window,min_delta,perc_delta,strand,fasta):
    #creates output bedgraph 3' end file based on threshold and option .fasta file with sequence upstream of each 3' end

    if perc_delta != None:
        peaks_all = [np.log10(x) for x in peak_dict.values()]
        perc = 10**np.percentile(peaks_all,perc_delta)
        sorted_peaks = [x for x in sorted(peak_dict.items(), key=operator.itemgetter(0)) if x[1] >= perc]
    else:
        sorted_peaks = [x for x in sorted(peak_dict.items(), key=operator.itemgetter(0)) if x[1] >= min_delta]

    if fasta.lower() == 'true':
        new_name_fasta = str(window)+"nt_"+str(fyle[:-9])+".fasta"
        with open(new_name_fasta,"w") as outp:
            for item in sorted_peaks:
                if strand.lower() == 'plus':
                    if item[0] <= window:
                        outp.write(">position:"+str(item[0])+'strand:+_Cv:'+str(item[1])+"\n")
                        edge = window - item[0]
                        gap = len(genome)-edge-1
                        seq = ''.join(genome[gap:])+''.join(genome[:item[0]])
                        outp.write(str(seq)+"\n")
                    else:
                        outp.write(">position:"+str(item[0])+'strand:+_Cv:'+str(item[1])+"\n")
                        seq = ''.join(genome[(item[0]-window-1):item[0]])
                        outp.write(str(seq)+"\n")

                elif strand.lower() == 'minus':
                    if item[0] <= len(genome)-(window):
                        outp.write(">position:"+str(item[0])+'strand:-_Cv:'+str(item[1])+"\n")
                        seq = genome[item[0]-1:item[0]+window]
                        rc_seq = reverse_complement(seq)
                        outp.write(str(rc_seq)+"\n")
                    else:
                        outp.write(">position:"+str(item[0])+'strand:-_Cv:'+str(item[1])+"\n")
                        edge = len(genome)-(item[0]-1)
                        gap = window - edge
                        seq = str(genome[item[0]-2:])+str(genome[:gap])
                        rc_seq = reverse_complement(seq)
                        outp.write(str(rc_seq)+"\n")


    new_name_bedgraph = "final_"+str(fyle)
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        contents = firstline.split("\t")
        sample = contents[0]
    with open(new_name_bedgraph,"w") as outp:
        for item in sorted_peaks:
            outp.write(str(sample)+"\t"+str(int(item[0]-1))+"\t"+str(int(item[0]))+"\t"+str(item[1])+"\n")

def main():
    parser = argparse.ArgumentParser(description='determines all Cv local maxima, outputs bedgraph file and optional .fasta file')
    parser.add_argument('strand',type=str,help='specify strand: <plus/minus>')
    parser.add_argument('-min_delta',type=int,default=10,help='min delta for thresholding Cv values, default value = 10')
    parser.add_argument('-perc_delta',type=float,default= None,help='min Cv percentile (based on log10 transformed set of values) for thresholding Cv values')
    parser.add_argument('-fasta_out',type=str,default='false',help='whether program outputs .fasta file  containing sequences upstream of all called 3-OH ends <true/false> default = false')
    parser.add_argument('-fasta',type=str,help='.fasta file used for alignment, only needed if -fasta_out = true')
    parser.add_argument('-fasta_window',type=int,default=50,help='length of upstream sequences in .fasta output file, default = 50, , only needed if -fasta_out = true')
    parser.add_argument('-file', default = None, help='Specific .bedgraph file to analyze')
    args = parser.parse_args()

    #for all Cv (delta) files in a directory, find local maxima, write to final .bedgraph file and optional .fasta file

    if args.file != None:
        if args.fasta_out.lower() == 'true':
            genome = read_fasta(args.fasta)
            sub_data = read_coverage(args.file)
            positions = read_position(args.file)
            peaks = local_max(sub_data,positions)
            writer(genome,peaks,args.file,args.fasta_window,args.min_delta,args.perc_delta,args.strand,args.fasta_out)
        else:
            sub_data = read_coverage(args.file)
            positions = read_position(args.file)
            peaks = local_max(sub_data,positions)
            writer('false',peaks,args.file,'false',args.min_delta,args.perc_delta,args.strand,'false')

    else:
        if args.fasta_out.lower() == 'true':
            genome = read_fasta(args.fasta)
            for fyle in sorted(glob.glob('*.bedgraph')):
                sub_data = read_coverage(fyle)
                positions = read_position(fyle)
                peaks = local_max(sub_data,positions)
                writer(genome,peaks,fyle,args.fasta_window,args.min_delta,args.perc_delta,args.strand,args.fasta_out)
        else:
            for fyle in sorted(glob.glob('*.bedgraph')):
                sub_data = read_coverage(fyle)
                positions = read_position(fyle)
                peaks = local_max(sub_data,positions)
                writer('false',peaks,fyle,'false',args.min_delta,args.perc_delta,args.strand,'false')

if __name__ == '__main__':
    main()
