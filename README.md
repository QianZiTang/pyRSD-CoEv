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
