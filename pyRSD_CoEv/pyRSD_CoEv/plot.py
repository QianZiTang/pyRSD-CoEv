#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator 

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		add_help=False,
		description=('Plot Relative snp density.'))
	
	requiredParser= parser.add_argument_group('requiredParser')
	
	requiredParser.add_argument('--input','-i',
		help='The rsd value file.',
		nargs=1,
		metavar='input file',
		type=str,
		required=True)

	requiredParser.add_argument('--indexfile','-d',
		help='The index file.',
		nargs=1,
		metavar='index file',
		type=str,
		required=True)

	requiredParser.add_argument('--strain','-s',
		help='The strain name.',
		nargs=1,
		metavar='strain name',
		type=str,
		required=True)
	
	requiredParser.add_argument('--outpre','-o',
		help='output file prefix',
		nargs=1,
		metavar='output file prefix',
		type=str,
		required=True)
	

	optionParser=parser.add_argument_group('optionParser')
	
	optionParser.add_argument('--help','-h',
		action="help",
		help="show this help message and exit")
	
	#optionParser.add_argument('--version',action='version',version='%(prog)s {}'.format(__version__))
	return parser

def getindex(indexfile):
    index={}
    f1=open(indexfile,'r')
    for line in f1:
        line=line.strip().split('\t')
        index[line[1]]=line[0]
    f1.close()
    return index

def plotman(inputdata,out):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    x_labels = []
    x_labels_pos = []
    #colors=['#333333','#999999']
    colors = ['#E24E42', '#008F95']
    
    #画散点图
    for num,(name, group) in enumerate(inputdata):
        plotcolor=num% len(colors)
        group.plot(x='indexnum',y='rsdvalue2',kind='scatter',color=colors[plotcolor], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['indexnum'].iloc[-1] - (group['indexnum'].iloc[-1] - group['indexnum'].iloc[0]) / 2))

    
    ax.set_xticks(x_labels_pos) 
    #添加x轴刻度上的值
    ax.set_xticklabels(x_labels)
    
    #设置y轴刻度范围以及显示的值
    #y_major_locator=MultipleLocator(50)
    #ax.yaxis.set_major_locator(y_major_locator)
    
    #设置图例并且设置图例的字体及大小
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontname('Times New Roman')
    
    #旋转90度
    plt.xticks(fontsize=10, rotation=60)
    plt.yticks(fontsize=10)
    
    #设置y轴和x轴的说明
    font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 10,}
    
    plt.xlabel('Chromosome',font)
    plt.ylabel('Relative SNV density(RSD)',font)
    plt.savefig(out,format='png',dpi=300,bbox_inches ="tight")

def getvalue(x,y,z):
    if x=='gap':
        return 0
    else:
        t=x.split(',')
        if z not in t:
        #if len(set(z)&set(t))==0:
            return 0
        else:
            i=t.index(z)
            v=y.split(',')[i]
            return v
def main():
	parse=parse_arguments()
	args=parse.parse_args()

	strainindex=getindex(args.indexfile[0])
	
	data=pd.read_table(args.input[0],sep='\t',header=None,names=['chromosome','start_bin','end_bin','strains_com','rsd_value'])
	data['strain_index']=strainindex[args.strain[0]]

	data['rsdvalue2'] = data[['strains_com','rsd_value','strain_index']].apply(lambda x:getvalue(x.strains_com,x.rsd_value,x.strain_index), axis = 1)
	data['chromosome']=data['chromosome'].str.replace('chr','')
	data.sort_values(by=['chromosome','start_bin','end_bin'])
	data['chromosome'] = data['chromosome'].astype('category')
	data['indexnum'] = range(len(data))
	data_grouped = data.groupby(('chromosome'))

	out=args.outpre[0]+'_plot.png'
	plotman(data_grouped,out)

if __name__ == '__main__':
	main()