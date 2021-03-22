# pyRSD-CoEv
An python package for selection sweep detection and co-evolutionary gene cluster identification

## Installation
1. python libraries
numpy>= 1.18.1
scipy>=1.4.1
networkx>=2.4
itertools
multiprocessing
subprocess
1. R packages
fdrtool



## 2. Identification of co-evolve gene cluster

## 2.1 Preparation
In order to use the coevgeneCluster function, you have to provide:
+ Variant effect prediction files. VEP or ANNOVAR output format are supported.

+ Sample index file
This file contain two colum separated by tab. The first cloumn is strains name, The second column is the location information of variant effect prediction files of this strain.
```
#example of sample index file
strain1	path/to/strain1_variant_effect_prediction_file
strain2	path/to/strain2_variant_effect_prediction_file
```
+ PCA file
First 15 principal components eigenvector values file. You can use smartpca to generate this file
```
#example of principal components eigenvector values file
eigvals	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15
strain1	-0.0420	-0.0697	0.0083	0.2046	0.1675	0.7422	0.4327	0.0497	0.3474	0.1036	0.0616	0.0155	0.0845	0.0372	0.0007
strain2	-0.0825	-0.1165	0.2635	-0.0447	-0.1126	-0.1236	0.0674	-0.0288	-0.0185	0.1594	0.1436	0.8599	-0.0110	0.2235	-0.0030
```
**！！ATTENTION**
The first column of sample index file a must be the same as the first column of file PCA

## 2.2 Usage
`coevgeneCluster --sampleindex < sampleindex file> --annotype < vep or ann> --output <prefix output file name> --varianttype <snv or indel> --pca <pca file> --rvalue <default 0.96> --pvalue <default 9.33e-10> --fdr <default 0.001>
## Required arguments
`--sampleindex,-s` 

Type: file. 
Sample index file

`--annotype,-a` 

Type: str. 
The software you used to predicted functional consequences of variants. Support VEP and ANNOVAR output format. Type 'vep' or 'ann' for this parameter.

`--output,-o`
Type: str.
Prefix of output file.

`--varianttype,-v`
Type: str.
The type of your variant data, snv or indel.

`--pca,-p`
PCA file
## Optional arguments
`--rvalue`
The cutoff of Pearson’s correlation coefficient， default value is 0.96.

`--pvalue`
The cutoff of p value, default value is 9.33e-10.

`--fdr`
The cutoff of fdr, default value is 0.001.
## 2.3 Output file
prefix_output_matrix.txt
(N+1)*M matrix. n+1:strain numer
prefix_output__re.txt

prefix_output_cor.txt
prefix_output_cor_fdr.txt
prefix_output_net.txt
