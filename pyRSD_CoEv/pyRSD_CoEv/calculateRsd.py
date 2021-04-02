#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import sys
import os
import re
import argparse
import itertools
from multiprocessing import Pool
from pyRSD_CoEv._version import __version__

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		add_help=False,
		description=('Calculate relative snp density in nonoverlapping windows.'))
	
	requiredParser= parser.add_argument_group('requiredParser')
	
	requiredParser.add_argument('--input','-i',
		help='The input file you need processed.',
		nargs=1,
		metavar='input file',
		type=argparse.FileType('r'),
		required=True)
	
	requiredParser.add_argument('--chrom','-c',
		help='The name of the chromosome. The name of the chromosome should be the same as that in your VCF file',
		nargs=1,
		metavar='chromosome name',
		type=str,
		required=True)
	
	requiredParser.add_argument('--sample_num','-s',
		help='The number of individual',
		nargs=1,
		metavar='individual number',
		type=int,
		required=True)

	optionParser=parser.add_argument_group('optionParser')
	optionParser.add_argument('--threads','-t',
		help='Number of threads. The minimum value for this parameter is 2',
		required=False,
		default=4,
		type=int)
	
	optionParser.add_argument('--cutoff','-u',
		help='SNV exists in at least a certain number of individuals. The default value is 1 ',
		metavar='cutoff number',
		type=int,
		default=1,
		required=False)

	optionParser.add_argument('--windowsize','-w',
		help='The sizes of nonoverlapping window. The default value is 10000',
		nargs=1,
		metavar='windowsize',
		type=int,
		default=10000,
		required=False)
	
	optionParser.add_argument('--help','-h',
		action="help",
		help="show this help message and exit")
	
	#optionParser.add_argument('--version',action='version',version='%(prog)s {}'.format(__version__))
	return parser

def readFai(fai):
	chr_len={}
	fai=open(fai,'r')
	for line in fai:
		line=line.strip().split('\t')
		chr_len[line[0]]=int(line[1])
	return chr_len

def rsdChisquare(row=None,i_index=None):
	total=sum(row)

	#each strain RSD value
	RSD_list=[i/total*100 for i in row]

	#Chi-square statistics were calculated for all possible combinations of strains
	#Only the combination corresponding to the minimum Chi-square value will keep
	smallest_value=100
	for i in range(1,len(row)+1):
		for subset in itertools.combinations(i_index,i):
			RSD = []
			subset_len=len(subset)
			Eval=float(100)/float(subset_len)
			each_total=0
			for k in subset:
				each_total+=((RSD_list[k]-Eval)**2)/Eval
				RSD.append([int(k),RSD_list[k]])
			if each_total<=smallest_value:
				smallest_value=each_total
				min_RSD=RSD

	min_RSD.sort()
	newmin_RSD=[]
	newmin_com=[]
	for each in min_RSD:
		newmin_RSD.append(str(each[1]))
		newmin_com.append(str(each[0]))

	return ','.join(newmin_RSD),newmin_com

def run_rsdChisquare(array=None,chr_name=None,cutoff=None,windowsize=10000,chrlength=None):
	start=array[0,0]
	end=array[array.shape[0]-1,0]
	num=array.shape[1]-1
	stuff=[i for i in range(num)]

	outfile_name='tpm_%s_%s_%s_RSD.txt'%(chr_name,start,end)
	outfile=open(outfile_name,'w')
	
	array_list=array.tolist()

	for row in array_list:
		if row[1:].count(0)<=num-cutoff:
			out_RSD,out_com=rsdChisquare(row[1:],stuff)
			if (row[0]+1)*windowsize<chrlength:
				start=row[0]*windowsize+1
				end=(row[0]+1)*windowsize
			if (row[0]+1)*windowsize>=chrlength:
				start=row[0]*windowsize+1
				end=chrlength
			outfile.write('%s\t%s\t%s\t%s\t%s\n'%(chr_name,start,end,','.join(out_com),out_RSD))
		else:
			if (row[0]+1)*windowsize<chrlength:
				start=row[0]*windowsize+1
				end=(row[0]+1)*windowsize
			if (row[0]+1)*windowsize>=chrlength:
				start=row[0]*windowsize+1
				end=chrlength
			outfile.write('%s\t%s\t%s\t%s\t%s\n'%(chr_name,start,end,'gap','novalue'))
	outfile.close()

def countSnv(inputfile=None,chr_name=None,sample_num=None,chrlength=None,windowsize=10000):
	'''
	:param inputfile: The outputfile of 01fliter.py
	:param chr_name: Chromosome
	:param outfile: Output file
	:return:
	'''
	#储存10kb的snv个数
	#chr_bin=readFai(fai)
	if chrlength%windowsize==0:
		binnum=chrlength//windowsize
	else:
		binnum=chrlength//windowsize+1
	allcount=np.zeros((binnum,int(sample_num)),dtype=int)

	i=0
	input=open(inputfile,'r')
	
	if windowsize/1000>=1:
		outfile='%s_snvnum_%skb.txt'%(chr_name,windowsize/1000)
	else:
		outfile='%s_snvnum_%sbp.txt'%(chr_name,windowsize)
	
	for line in input:
		line=line.strip().split('\t')
		if line[0]==chr_name:
			gt=line[3].split(' ')
			if len(gt)!=int(sample_num):
				sys.exit("please check your sample number and inputfile sample number!")
			else:
				row_index=int(line[1])//windowsize
				for i in range(0,len(gt)):
					if gt[i]=='11':
						allcount[row_index,i]+=1
					
	input.close()
	np.savetxt(outfile,allcount,fmt='%d')
	return allcount
	
def main():
	parse=parse_arguments()
	args=parse.parse_args()
	
	#Check whether the chromosome size file exists
	pwd_dir=os.path.abspath(os.path.split(args.input[0].name)[0])
	fai_file=os.path.join(pwd_dir,'chrom_sizes.txt')
	if not os.path.exists(fai_file):
		sys.exit('No chromosome size information file (chrom_sizes.txt) was found in the {0} '.format(pwd_dir))
	
	#mkdir output dir
	out_dir=pwd_dir+'/02caculatersd'
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	os.chdir(out_dir)
	
	allchr_length=readFai(fai_file)
	allcount=countSnv(inputfile=args.input[0].name,chr_name=args.chrom[0],sample_num=args.sample_num[0],chrlength=allchr_length[args.chrom[0]],windowsize=args.windowsize)


	#切割文件计算RSD值
	if args.threads<2:
		print('The minimum value of --thread is 2')
		args.threads=2
	newallcount=np.insert(allcount,0,range(allcount.shape[0]),axis=1)
	tpm=np.array_split(newallcount,args.threads)
	new_tpm=[]
	file_index=[]
	#run_rsdChisquare(array=None,chr_name=None,cutoff=None,windowsize=10000,chrlength=None)
	for i in range(args.threads):
		new_tpm.append((tpm[i],args.chrom[0],args.cutoff,args.windowsize,allchr_length[args.chrom[0]]))
		file_index.append((tpm[i][0,0],tpm[i][tpm[i].shape[0]-1,0]))

	pool = Pool(processes=args.threads)
	pool.starmap(run_rsdChisquare, new_tpm)
	pool.close()
	pool.join()

	#合并切割文件
	output=open('%s_RSD.txt'%(args.chrom[0]),'wb')
	for i in file_index:
		f = open('tpm_%s_%s_%s_RSD.txt'%(args.chrom[0],i[0],i[1]), 'rb')
		output.write(f.read())
		f.close()
		os.remove('tpm_%s_%s_%s_RSD.txt'%(args.chrom[0],i[0],i[1]))
	output.close()


if __name__ == '__main__':
	main()
	

