#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import numpy
import os,subprocess
from pyRSD_CoEv._version import __version__

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		add_help=False,
		description=('Segments were then merged if successive 10 kb windows were shared between same strain combinations.'))
	
	requiredParser= parser.add_argument_group('requiredParser')
	
	requiredParser.add_argument('--input','-i',
		help='The rsd value file.',
		nargs=1,
		metavar='input file',
		type=str,
		required=True)
	
	requiredParser.add_argument('--outpre','-o',
		help='output file prefix',
		nargs=1,
		metavar='output file prefix',
		type=str,
		required=True)
	

	optionParser=parser.add_argument_group('optionParser')

	optionParser.add_argument('--shuffle','-s',
		help='Number of shuffle. The default value for this parameter is 100000',
		required=False,
		default=100000,
		type=int)
	
	optionParser.add_argument('--help','-h',
		action="help",
		help="show this help message and exit")
	
	#optionParser.add_argument('--version',action='version',version='%(prog)s {}'.format(__version__))
	return parser

def mergeRsdWindows(inputfile=None,mergeout=None):
	'''
	'''
	all_combine=[]
	needpvalue={}
	mergedinfo=[]
	
	count=0
	f=open(inputfile,'r')
	o=open(mergeout,'w')
	
	for line in f:
		newline=line.strip().split('\t')
		all_combine.append(newline[3])
		chrom=newline[0]
		if count==0:
			win=int(newline[2])-int(newline[1])+1
		count+=1

	old_start=0
	old_end=0

	i=0
	while i<len(all_combine):
		new_combine=all_combine[i]
		if i==0:
			old_combine=new_combine
			old_start=i
			old_end=i+1
			length=1
			i+=1
		elif i>0:
			if new_combine==old_combine:
				old_end=i+1
				length+=1
				i=i+1
			else:
				if length>=2 and old_combine!='gap' and (i+2)<len(all_combine) and all_combine[i+1]==old_combine and all_combine[i+2]==old_combine:
					old_end=i+2+1
					new_combine=old_combine
					length+=3
					i=i+3
				else:
					if length>=2 and old_combine!='gap':
						o.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(chrom,old_start*win+1,(old_end)*win,old_combine,length))
						needpvalue.setdefault(old_combine,{})
						needpvalue[old_combine].setdefault(length,0)
						mergedinfo.append([chrom,old_start*win+1,(old_end)*win,old_combine,length])
					old_start=i
					old_end=i+1
					old_combine=new_combine
					length=1
					i=i+1
	if length>=2 and old_combine!='gap':
		needpvalue.setdefault(old_combine,{})
		needpvalue[old_combine].setdefault(length,0)
		mergedinfo.append([chrom,old_start*win+1,(old_end)*win,old_combine,length])
		o.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(chrom,old_start*win+1,(old_end)*win,old_combine,length))
	f.close()
	o.close()
	return all_combine,needpvalue,mergedinfo

def smooth(all_combine=[],needpvalue={},mergedinfo=[],shuffle=100000,output=None):
	'''
	estimated the probability of 
	occurrence of long genomic regions 
	that were unique to a given strain or to a number of strain
	combinations using a permutation test (n=10 million)
	
	all_combine:type list. contain all 10kb combinations. 
	
	needpvalue:type dict. merged successive 10 kb windows combinations, value:length
	
	shuffle: permutation number

	needpvalue={strain_combin:{length:0}}
	'''
	#contain rsdvaluefile strain combination
	
	# b_combin=[]
	# rsdfile=open(rsdvaluefile,'r')
	# for line in rsdfile:
	# 	newline=line.strip().split('\t')
	# 	b_combin.append(newline[3])
	# rsdfile.close()
	
	filelength=len(all_combine)
	


	#mergedfile=open(mergedrsdfile,'r')
	# needpvalue={}
	# mergedinfo=[]
	
	# for line in mergedfile:
	# 	newline=line.strip().split('\t')
		
	# 	length=newline[4]
	# 	strain_combin=newline[3]
		
	# 	if strain_combin in needpvalue:
	# 		if length not in needpvalue[strain_combin]:
	# 			needpvalue[strain_combin][length]=0
	# 	else:
	# 		needpvalue[strain_combin]={}
	# 		needpvalue[strain_combin][length]=0
	# 	#needpvalue.setdefault(strain_combin,{})
		
		
	# 	mergedinfo.append(line.strip())
	# mergedfile.close()
	
	

	out=open(output,'w')
	for shuffle_num in range(0,shuffle):
		longest={}
		for c in needpvalue.keys():
			longest[c]=[0]

		mylist=list(numpy.random.permutation(filelength))
		mycombs=[]
		for k in mylist:
			mycombs.append(all_combine[k])

		i=0
		while i<len(mycombs):
			new_combine=mycombs[i]
			if i==0:
				old_combine=new_combine
				#old_start=i
				#old_end=i+1
				length=1
				i+=1
			elif i!=0:
				if new_combine==old_combine:
					#old_end=i+1
					length+=1
					i=i+1
				elif new_combine!=old_combine:
					if length>=2 and old_combine!='gap' and (i+2)<len(mycombs) and mycombs[i+1]==old_combine and mycombs[i+2]==old_combine:
						#old_end=i+2+1
						new_combine=old_combine
						length+=3
						i=i+3
					else:
						if length>=2 and old_combine!='gap' and (old_combine in needpvalue):
							longest[old_combine].append(length)
						#old_start=i
						#old_end=i+1
						old_combine=new_combine
						length=1
						i=i+1
		if length>=2 and old_combine!='gap' and (old_combine in needpvalue) :
			longest[old_combine].append(length)

		newlongest={}
		for c in longest:
			maxval=max(longest[c])
			newlongest[c]=maxval
				
		for c in needpvalue:
			maxval=newlongest[c]
			for d in needpvalue[c]:
				if maxval>=int(d):
					needpvalue[c][d]+=1
	
	for i in mergedinfo:
		pvalue=float(needpvalue[i[-2]][i[-1]])/float(shuffle)
		out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(i[0],i[1],i[2],i[3],len(i[3].split(',')),i[4],pvalue))
	out.close()

def correctPvalue(inputfile,out):
	R = subprocess.getoutput('which Rscript')

	print (R)
	if 'no Rscript in' in R:
		sys.exit(f'Cannot find R in your $PATH, you could use -r to to specify the R path')
	else:
		path=os.path.dirname(os.path.abspath(__file__))
		cmd='Rscript %s %s %s %s'%(os.path.join(path,'fdr.R'),0,inputfile,out)
		os.system(cmd)
def main():
	parse=parse_arguments()
	args=parse.parse_args()

	print('***start merge window***')
	mergefile=args.outpre[0]+'_merged.txt'
	all_combin,needpvalue,mergedinfo=mergeRsdWindows(args.input[0],mergefile)
	print('***finish merge window***')

	print('***start shuffle file***')
	shufflefile=args.outpre[0]+'_merged_shuffled.txt'
	smooth(all_combin,needpvalue,mergedinfo,args.shuffle,shufflefile)
	print('***finish shuffle file***')

	print('***strat correct P value***')
	correctfile=args.outpre[0]+'_merged_shuffled_fdr.txt'
	correctPvalue(shufflefile,correctfile)
	print('***finish correct P value***')


if __name__ == '__main__':
	main()




	





