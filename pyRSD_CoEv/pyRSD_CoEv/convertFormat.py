import sys
import os
import re
import argparse
#import logging
from pyRSD_CoEv._version import __version__
def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		add_help=False,
		description=('Fliter out the heterozygotes and missing value. Create a RSDer input format file.'))
	requiredParser= parser.add_argument_group('requiredParser')
	requiredParser.add_argument('--vcf','-v',
		help='The vcf file you need processed.',
		nargs=1,
		metavar='vcf file',
		type=argparse.FileType('r'),
		required=True)
	optionParser=parser.add_argument_group('optionParser')
	optionParser.add_argument('--help','-h',
		action="help",
		help="show this help message and exit")
	optionParser.add_argument('--version',action='version',version='%(prog)s {}'.format(__version__))
	return parser

def keepSnp(ref,alt,gttype):
	issnp=0
	heter=0
	hasmiss=0
	hasN=0
	homo=0
	if len(ref)==1 and len(alt)==1:
		issnp=1
	if '01' in gttype or '10' in gttype:
		heter=1
	if '.' in ''.join(gttype):
		hasmiss=1
	if '11' in gttype:
		homo=1
	if 'N' in ''.join(gttype) or 'N' in ref or 'N' in alt:
		hasN=1
	if homo==1 and issnp==1 and heter==0 and hasmiss==0 and hasN==0:
		return True


def readVcf(vcffile):
	'''
	:param vcffile: vcf format input file
	:return: chrom size dict and flitered line
	'''
	all_chr_len={}
	f=open(vcffile,'r')
	
	dir_path=os.path.abspath(os.path.dirname(vcffile))
	outpath=os.path.join(dir_path,'01fliteredvcf')
	
	if not os.path.exists(outpath):
		os.mkdir(outpath)

	outfile_name=os.path.basename(vcffile).replace('.vcf','')+'_fliter.txt'
	o=open(os.path.join(outpath,outfile_name),'w')
	chrfile=open(os.path.join(dir_path,'chrom_sizes.txt'),'w')

	for line in f:
		line=line.strip()
		if line.startswith('#'):
			if line.startswith('##contig'):
				pattern1=re.compile('##contig=<ID=(.+),length')
				pattern2=re.compile('length=(.+)>')
				chr_ID=pattern1.findall(line)[0]
				chr_len=int(pattern2.findall(line)[0])
				all_chr_len[chr_ID]=chr_len
				chrfile.write('%s\t%s\n'%(chr_ID,chr_len))
			if line.startswith('#CHROM'):
				ind_names=line.split('\t')[8:]
		else:
			newline=line.split('\t')
			chrom=newline[0]
			pos=newline[1]
			gt=[i.split(':')[0].replace('/','') for i in newline[9:]]
			ref=newline[3]
			alt=newline[4]
			if keepSnp(ref,alt,gt)==True:
				o.write('%s\t%s\t%s\\%s\t%s\n'%(chrom,pos,ref,alt,' '.join(gt)))
	f.close()
	o.close()
	chrfile.close()

def main():
	parse=parse_arguments()
	args=parse.parse_args()
	if not args.vcf[0].name.endswith('.vcf'):
		sys.exit('Please defin the extension of inputfile. VCF format is supported.')
	else:
		readVcf(args.vcf[0].name)

if __name__ == '__main__':
	main()
