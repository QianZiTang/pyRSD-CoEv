#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import argparse
def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=('annotation genes located in fdr<0.001 genome region identified by RSDer'))
    requiredParser = parser.add_argument_group('requiredParser')

    requiredParser.add_argument('--gtffile','-g',
                                help='GTF format gene annotation file.',
                                nargs=1,
                                metavar='GTF file',
                                type=str,
                                required=True)

    requiredParser.add_argument('--inputfile','-i',
                                help='input file.',
                                nargs=1,
                                metavar='input file',
                                type=str,
                                required=True)

    requiredParser.add_argument('--samplefile','s',
                                help='sample file',
                                nargs=1,
                                metavar='sample file',
                                type=str,
                                required=True)

    requiredParser.add_argument('--fai','-f',
                                help='chrom length file',
                                nargs=1,
                                metavar='chrom length file',
                                type=str,
                                required=True)


    requiredParser.add_argument('--output','-o',
                                help='output file',
                                nargs=1,
                                metavar='output file',
                                type=str,
                                required=True)

    optionParser = parser.add_argument_group('optionParser')

    optionParser.add_argument('--help', '-h',
                              action="help",
                              help="show this help message and exit")
def extractGene(gtffile):
    f=open(gtffile,'r')

    all_gene={}
    for line in f:
        newline=line.strip().split('\t')
        eachparts = str(parts[8]).split(';')
        gene_id = (eachparts[0]).split('"')[1]
        start = int(newline[3])
        end = int(newline[4])
        strand = str(newline[6])
        if 'chr' not in newline[0]:
            chrom='chr'+newline[0]
        else:
            chrom =newline[0]

        if newline[1] == 'protein_coding' and newline[2] == 'gene':
            all_gene.setdefault(chrom,[])
            all_gene[chrom].append([start,end,gene_id,strand])
    f.close()
    return all_gene

def annotationGene(inputfile,samplefile,fai,genelist,outfile):
    inputfile1=open(inputfile,'r')
    samplefile=open(samplefile,'r')
    fai=open(fai,'r')
    out=open(outfile,'w')

    sampleindex=[line.strip() for line in samplefile]
    samplefile.close()

    chr_len={}
    for line in fai:
        newline=line.strip().split('\t')
        chr_len[newline[0]]=int(newline[1])
    fai.close()

    for line in inputfile:
        newline=line.strip().split('\t')
        start=int(newline[1])-5000
        end=int(newline[2])+5000
        if 'chr' not in newline[0]:
            chrom = 'chr' + newline[0]
        else:
            chrom = newline[0]

        start= max(1, start)
        end=min(end,chr_len[chrom])

        tpm=[]
        for l in genelist[chrom]:
            gen_start=l[0]
            gene_end=l[1]
            gene_id=l[chrom][2]

            if gen_start<=end and start<=gene_end:
                overlap=min(end,gene_end)-max(start,gen_start)+1
                tpm.append(gene_id)
        combin=newline[3]
        newcombin=[sampleindex[int(i)] for i in combin]

        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(chrom,newline[1],newline[2],';'.join(newcombin),';'.join(tpm),'\t'.join(newline[4:])))
    out.close()

def main():
    parse = parse_arguments()
    args = parse.parse_args()

    genelist=extractGene(args.gtffile[0])
    annotationGene(args.inputfile[0], args.samplefile[0], args.fai[0], genelist, args.outfile[0])







