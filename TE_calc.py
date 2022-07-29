#!/usr/bin/env python3

import argparse
import numpy as np
from operator import itemgetter
import codecs

def TE_calc(coord,cov,strand,upstream,downstream):
	
	if strand == '+':
		up = np.median(cov[coord-upstream:coord])
		down = np.median(cov[coord+1:coord+downstream])

		if up >= 10:
			TE = round(((up-down)/up)*100,2)
			if TE < 0.0:
				TE = 0.0 
			return TE
		else:
			return None

	else:
		up = np.median(cov[coord-1:coord+(upstream-1)])
		down = np.median(cov[coord-(downstream+1):coord-1])

		if up >= 10:
			TE = round(((up-down)/up)*100,2)
			if TE < 0.0:
				TE = 0.0 
			return TE
		else:
			return None
	
def read_coverage(cov_fyle):
	#create list of coverage values, position dependent
	return_list = []
	with open(cov_fyle,"r") as inp:
		for line in inp:
			return_list.append(float(line.strip().split('\t')[2]))
	return return_list

def read_ends(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			coord = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
			return_dict[coord] = strand
	return return_dict

def main():
	parser = argparse.ArgumentParser(description='for all 3-OH ends in a CSV, calculates TE in 2 biological conditions. Only outputs ends with %T above threshold in the WT strain.')
	parser.add_argument('ends',type=str,help='CSV of 3-OH ends')
	parser.add_argument('threshold',type=int,help='TE threshold, must be positive integer')
	parser.add_argument('WT_fwd',type=str,help='fwd cov file of WT strain [.cov]')
	parser.add_argument('WT_rev',type=str,help='rev cov file of WT strain [.cov]')
	parser.add_argument('Mutant_fwd',type=str,help='fwd cov file of Mutant strain [.cov]')
	parser.add_argument('Mutant_rev',type=str,help='rev cov file of Mutant strain [.cov]')
	parser.add_argument('upstream',type=int,help='upstream distance (nt) for TE calc, must be positive intenger')
	parser.add_argument('downstream',type=int,help='downstream distance (nt) for TE calc, must be positive intenger')
	parser.add_argument('output',type=str,help='')

	args = parser.parse_args()

	ends = read_ends(args.ends)
	WT_fwd = read_coverage(args.WT_fwd)
	WT_rev = read_coverage(args.WT_rev)
	Mutant_fwd = read_coverage(args.Mutant_fwd)
	Mutant_rev = read_coverage(args.Mutant_rev)

	final_list = []
	for key,value in ends.items():
		if value == '+':
			TE_WT = TE_calc(key,WT_fwd,'+',args.upstream,args.downstream)
			TE_Mutant = TE_calc(key,Mutant_fwd,'+',args.upstream,args.downstream)
			
			if TE_WT != None and TE_Mutant != None:
				if TE_WT >= args.threshold:
					final_list.append([key,value,TE_WT,TE_Mutant,round((TE_WT-TE_Mutant),0)])
		else:
			TE_WT = TE_calc(key,WT_rev,'-',args.upstream,args.downstream)
			TE_Mutant = TE_calc(key,Mutant_rev,'-',args.upstream,args.downstream)
			
			if TE_WT != None and TE_Mutant != None:
				if TE_WT >= args.threshold:
					final_list.append([key,value,TE_WT,TE_Mutant,round((TE_WT-TE_Mutant),0)])

	final_list_sorted = sorted(final_list, key=itemgetter(0))

	with open(args.output,'w',encoding='utf-8',newline='') as outp:
		outp.write('POT,Strand,%T WT,%T Mutant, d%T\n')
		for item in final_list_sorted:
			outp.write(','.join(list(map(str,item)))+'\n')

		
if __name__ == '__main__':
	main()




	