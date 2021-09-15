#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import argparse
from scipy.stats.stats import pearsonr
import networkx as nx
from pyRSD_CoEv._version import __version__

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		add_help=False,
		description=(''))
	
	requiredParser= parser.add_argument_group('requiredParser')
	
	'''requiredParser.add_argument('--inputpath','-i',
		help='The path to the inputfile.NOTICE: the path should only contain inputfile',
		nargs=1,
		metavar='input file path',
		type=str,
		required=True)'''
	
	
	requiredParser.add_argument('--sampleindex','-s',
		help='strains name file.',
		nargs=1,
		metavar='strains name file',
		type=str,
		required=True)

	requiredParser.add_argument('--annotype','-a',
		help='The software you used to predicted functional consequences of variants. Support VEP and ANNOVAR output format. Type vep or ann for this parameter.',
		nargs=1,
		metavar='vep or ann',
		type=str,
		required=True)

	requiredParser.add_argument('--output','-o',
		help='output file. ',
		nargs=1,
		metavar='output file',
		type=str,
		required=True)

	requiredParser.add_argument('--varianttype','-v',
		help='variant type',
		nargs=1,
		metavar='variant type: indel or snv',
		type=str,
		required=True)	

	requiredParser.add_argument('--pca','-p',
		help='pca file.',
		nargs=1,
		metavar='pca file',
		type=str,
		required=True)

	optionParser=parser.add_argument_group('optionParser')
	
	optionParser.add_argument('--help','-h',
		action="help",
		help="show this help message and exit")
	
	return parser
def read_indeldata(inputfile,filetype):
	transcript2frameshift={}
	frameshift_type=['frameshift insertion','frameshift deletion','frameshift block substitution','frameshift_variant']

	f = open(inputfile, 'r')
	for line in f:
		if line.startswith('#'):
			pass
		else:
			line.rstrip()
			parts = line.split('\t')
			if filetype=='vep':
				transcriptname = [parts[4]]
				conseq=parts[6].split(',')
			if filetype=='ann':
				transcriptname = [i.split(':')[1] for i in parts[2].split(',')[:-1]]
				conseq=parts[1]

			if len(set(conseq)&set(frameshift_type))==1:
				for i in transcriptname:
					transcript2frameshift.setdefault(i,0)
					transcript2frameshift[i]+=1
	f.close()
	return transcript2frameshift

def buldindelMtrix(indexfile,filetype,output):
	files=open(indexfile,'r')
	out=open(output,'w')
	
	all_transcript2rate=[]
	union_names={}
	all_samplename=[]
	
	for i in files:
		i=i.strip().split('\t')
		file=i[1]
		tpm=read_indeldata(file,filetype)
		all_transcript2rate.append(tpm)
		all_samplename.append(i[0])
		for k in tpm:
			union_names.setdefault(k,0)
	files.close()
	
	for c in union_names:
		rates=[]
		for dictor in all_transcript2rate:
			if c in dictor:
				rate=dictor[c]
			else:
				rate=0
			rates.append(str(rate))
		out.write(c+'\t'+'\t'.join(rates)+'\n')
	out.close()




def read_snvdata(inputfile,filetype):
	transcript2rate={}
	transcript2NSC={}
	transcript2SC={}
	
	nsctype = ['nonsynonymous SNV', 'stopgain SNV', 'stoploss SNV', 'stopgain', 'stoploss', 'stop_lost','stop_gained', 'missense_variant', 'start_lost']
	sctype = ['synonymous SNV', 'synonymous_variant']
	
	f = open(inputfile, 'r')
	
	for line in f:
		if line.startswith('#'):
			pass
		else:
			line.rstrip()
			parts = line.split('\t')
			if filetype=='vep':
				transcriptname = [parts[4]]
				conseq=parts[6].split(',')
			if filetype=='ann':
				transcriptname = [i.split(':')[1] for i in parts[2].split(',')[:-1]]
				conseq=parts[1]

			if len(set(conseq)&set(nsctype))==1:
				for i in transcriptname:
					transcript2NSC.setdefault(i,0)
					transcript2NSC[i]+=1
			if len(set(conseq)&set(sctype))==1:
				for i in transcriptname:
					transcript2SC.setdefault(i,0)
					transcript2SC[i]+=1
	f.close()
	
	names_union=list(set(transcript2SC.keys()).union(set(transcript2NSC.keys())))
	for c in names_union:
		if c in transcript2NSC:
			NSC=transcript2NSC[c]
		else:
			NSC=0

		if c in transcript2SC:
			SC=transcript2SC[c]
		else:
			SC=0
		
		if SC!=0:
			#Include 0/X, X/X
			rate=float(NSC)/float(SC)
		else:
			rate=float(NSC)
		transcript2rate[c]=rate
	return transcript2rate

def buldsnvMtrix(indexfile,filetype,output):
	'''
	indexfile:
	starin_name1	starin1_vepfile
	starin_name2	strain2_vepfile
	'''
	files=open(indexfile,'r')
	out=open(output,'w')
	
	all_transcript2rate=[]
	all_samplename=[]
	union_names={}
	for i in files:
		i=i.strip().split('\t')
		file=i[1]
		tpm=read_snvdata(file,filetype)
		all_transcript2rate.append(tpm)
		all_samplename.append(i[0])

		for k in tpm:
			union_names.setdefault(k,0)

	for c in union_names:
		has_NSC=0
		has_NA=0
		for dictor in all_transcript2rate:
			if c in dictor:
				if (dictor[c]!='NA') and (dictor[c]!=0):
					has_NSC=1
				if dictor[c]=='NA':
					has_NA=1
		if has_NSC==1 and has_NA==0:
			rates=[]
			for dictor in all_transcript2rate:
				if (c in dictor):
					rate=dictor[c]
				else:
					rate=0
				rates.append(str(rate))
			out.write(c+'\t'+'\t'.join(rates)+'\n')
	out.close()

def lmregression(ratiofile,pcaevc,output):
	R = subprocess.getoutput('which Rscript')
	#outpath=os.path.dirname(os.path.abspath(ratiofile))
	print (R)
	if 'no Rscript in' in R:
		sys.exit(f'Cannot find R in your $PATH, you could use -r to to specify the R path')
	else:
		path=os.path.dirname(os.path.abspath(__file__))
		cmd='Rscript %s %s %s %s'%(os.path.join(path,'regression.R'),ratiofile,pcaevc,output)
		os.system(cmd)

def getcor(regrefile,output):
	values=[]
	names=[]
	o=open(output,'w')
	refile=open(regrefile,'r')
	for line in refile:
		line=line.rstrip()
		parts=line.split('\t')
		each=[float(i) for i in parts[1:]]
		values.append(each)
		names.append(parts[0])

	k=0
	while k<(len(names)-1):
		for j in range(k+1,len(names)):
			try:
				corr_result=pearsonr(values[k],values[j])
				correlation_coeff=corr_result[0]
				if float(corr_result[0])>0:
					pvalue=float(corr_result[1])/2.0
				else:
					pvalue=1-float(corr_result[1])/2.0
				if str(correlation_coeff)!='nan':
					o.write(names[k]+'\t'+names[j]+'\t'+str(correlation_coeff)+'\t'+str(pvalue)+'\n')
			except:
				print(names[k],names[j],'invalid value')
		k+=1
	o.close()
	refile.close()

def getfdr(correlationfile,out):
	R = subprocess.getoutput('which Rscript')
	path=os.path.dirname(os.path.abspath(__file__))
	cmd='Rscript %s %s %s %s'%(os.path.join(path,'fdr.R'),"1",correlationfile,out)
	os.system(cmd)

def getnetwork(input=None,r=0.96,pvalue=9.33e-15,fdr=0.001,out=None):
	f=open(input,'r')
	o=open(out,'w')
	edges =[]
	for line in f:
		if not line.startswith('#'):
			parts=line.strip().split('\t')
			if float(parts[4])<=fdr and abs(float(parts[2]))>=r and float(parts[3])<=pvalue:
				edges.append((parts[0],parts[1]))

	G = nx.Graph()
	G.add_edges_from(edges)
	G_cluster=list(nx.connected_components(G))
	count=1
	
	for i in G_cluster:
		o.write('cluster_%s\t%s\t%s\n'%(count,','.join(list(i)),len(i)))
		count+=1
	f.close()
	o.close()

def main():
	parse=parse_arguments()
	args=parse.parse_args()
	
	print("***Start build matrix***")
	outputmatrix=args.output[0]+'_matrix.txt'
	if args.annotype[0] not in ['vep','ann']:
		sys.exit(f'-f only support vep or ann')
	if args.varianttype[0] not in ['snv','indel']:
		sys.exit(f'-v only support snv or indel')
	
	'''if args.varianttype[0]=='snv':
		buldsnvMtrix(args.sampleindex[0],args.annotype[0],outputmatrix)
	if args.varianttype[0]=='indel':
		buldindelMtrix(args.sampleindex[0],args.annotype[0],outputmatrix)
	
	print("STEP2: Regression")
	#outputreg=args.output[0]+'_regre.txt'
	outputreg=args.output[0]+'_re.txt'
	print(outputreg)
	lmregression(outputmatrix,args.pca[0],outputreg)
	
	print("STEP3: Pairwise correlation")
	outputcor=args.output[0]+'_re_cor.txt'
	getcor(outputreg,outputcor)

	print("STEP4: Adjust p value")'''
	outfdr=args.output[0]+'_re_cor_fdr.txt'
	#getfdr(outputcor,outfdr)

	print("STEP5: construct gene cluster")
	outnet=args.output[0]+'_net.txt'
	getnetwork(input=outfdr,r=0.96,pvalue=9.33e-15,fdr=0.001,out=outnet)

	
if __name__ == '__main__':
	main()
