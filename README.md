# pyRSD-CoEv
An python package for selection sweep detection and co-evolutionary gene cluster identification


## 00 Installation
## Requirement
+ numpy>= 1.18.1
+ scipy>=1.4.1
+ networkx>=2.4
+ itertools
+ multiprocessing
+ subprocess
+ fdrtool(R package)

```
git https://github.com/QianZiTang/pyRSD-CoEv.git
cd pyRSD-CoEv
python setup.py install
```
## 1. Identification of selective sweeps using relative homozygous SNV density

### 1.1 Tools for pre-processing variant call format SNV data：` convertFormat`

The tool is designed for filtering polymorphic sites of missing or heterozygote genotypes in at least one individual. This step generates two result files (chrom_sizes.txt and input_filename_filter.txt) in the input file path. chrom_sizes.txt contains chromosome length information. input_filename_filter.txt is a tab delimited file which contains information about polymorphic sites. The first column indicates chromosome name, the second indicates the position in the chromosome, and the third column is REF/ALT (reference allele/alternate allele). The remaining columns include the genotypical information for each individual. input_filename_filter.txt is used as input file for` calculateRSD` function. user can download example file in https://github.com/QianZiTang/pyRSD-CoEv/releases/download/v1.0-alpha/example.rar

- Usage
`convertFormat –vcf <vcf_filename>`
- Required arguments
`--vcf, -v <vcf_filename>` This argument defines the VCF file to be processed. VCF format v4.1 is supported.
- Example usage
`convertFormat --vcf /example/example.vcf`
- Output file
1) /example/01fliteredvcf/input_filename_filter.txt
```
chr15		33506	G/A	00 00 11 00 00 00 00 00 00 11 11 00 00 00 00 00 00 00 11 00 00 00
chr15	33901	G/A	00 00 00 11 11 00 11 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
```
This file contains information about polymorphic sites shown in 4 columns. The first and second columns represent genomic coordinates of the variant. The third column is the same as the REF/ALT column in the inputted VCF file. The fourth column contains genotypical information for each strain separated by space. This file should be used as input file for the `calculateRSD` function. 

 2) /example/01fliteredvcf /chrom_sizes.txt
This file contains chromosome length information. If there is no header information in your VCF file, please replace it with your own chromosome length file at the end of this step and name it as ‘chrom_sizes.txt’
 3) /example/01fliteredvcf /index.txt
Index information of strains.

###  1.2 Tools for calculating relative snp density in nonoverlapping windows : `calculateRsd`
- Usage
`calculateRSD --input <input_filename> --chrom <chromosome> --sample_num <integer>`
- Required arguments
`--input, -i <input_filename>`
convertFormat result file (input_filename_filter.txt) can be used as inputfile for this function.
`--chrom, -c <chromosome>`
Chromosome name. eg chr1
`--sample_num, -s <integer>`
The number of individuals.
- Optional arguments
`--threads, -t < number of threads, the default value is 2>`
`--cutoff, -u < The minimum number of strains with SNV. The default value is 1>`
`--windowsize, -w <default value is 10000>`
- Example Usage
`
calculateRsd -i /example/01fliteredvcf/input_filename_filter.txt --chrom chr15 --sample_num 22 -t 6
`
- outputfile
1) /example/02caculatersd/chr15_RSD.txt
This file contains information about relative homozygous SNV density (RSD) for each non-overlapping 10-kb window. The first column is the chromosome name. The second and third columns are the start and end coordinates of the window. The fourth column is the combination of strains, and the fifth column is the  RSD values of strains.

 ** This function generates intermediate files (/example/02caculatersd/tpm_chr*_*_*_RSD.txt) when it is running. After running, the intermediate file will be automatically deleted and the final result file (/example/02caculatersd/chr*_RSD.txt) will be generated.

### 1.3	Tools for merging and smoothing `mergeANDsmooth`
- Usage
`mergeANDsmooth --input <input file> --output <output prefix>`
- Required arguments
` --input, -i <input file> `
` --output, -o <output prefix>`
- Example usage
`mergeANDsmooth -i /example/02caculatersd/chr15_RSD.txt –o /example/02caculatersd/chr15_RSD`
- outputfile
/example/02caculatersd/chr15_RSD/chr15_RSD_merged.txt
/example/02caculatersd/chr15_RSD/chr15_RSD_merged_shuffled.txt
/example/02caculatersd/chr15_RSD/chr15_RSD_merged_shuffled_fdr.txt

The final output file of this function. Only merged genomic regions with p value<1e-5 and FDR<0.001 will be kept. The first and fourth columns are same as the output file of calculateRsd. The second column and the third column are the new start and end coordinates after merging, and the sixth column is the number of merged 10 kb windows. The fourth column is the combination of strains, and the fifth column is the number of strains in this combination. Columns 7 to 9 are p, Q-val and FDR, respectively.

### 1.4 Tools for identifying coevolutionary gene clusters `coevgeneCluster`
In order to use the `coevgeneCluster` function, you have to provide:
- Variant effect prediction files. VEP or ANNOVAR output formats are supported.
- Sample index file. This file contains two columns separated by tab. The first column is the strain name, the second column is the folder path of variant effect prediction files for this strain.
```
#example of sample index file
strain1	path/to/strain1_variant_effect_prediction_file
strain2	path/to/strain2_variant_effect_prediction_file
```
- PCA file. Eigenvector values for the first 15 principal components. You can use smartpca to generate this file
```
#example of eigenvector values for principal components
eigvals	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15
strain1	-0.0420	-0.0697	0.0083	0.2046	0.1675	0.7422	0.4327	0.0497	0.3474	0.1036	0.0616	0.0155	0.0845	0.0372	0.0007
strain2	-0.0825	-0.1165	0.2635	-0.0447	-0.1126	-0.1236	0.0674	-0.0288	-0.0185	0.1594	0.1436	0.8599	-0.0110	0.2235	-0.0030
```
- Example Usage
`coevCluster --sampleindex <sampleindex file> --annotype <vep or ann> --output < prefix_output_matrix.txt> –varianttype <snv or indel> --pca <pca file> --rvalue <default 0.96> --pvalue <default 9.33e-10> --fdr <default 0.001>`

- outputfile
1) prefix_output_matrix.txt
(N+1)*M matrix. 
N+1: the first column is transcript, N is the number of strains. M: the number of transcripts showing NSC or frameshift
2) prefix_output_re.txt
corrected (N+1)*M matrix
3) prefix_output_cor.txt
Transcript pairs Pearson’s correlation and p value
4) prefix_output_cor_fdr.txt
Transcript pairs with Pearson’s correlation, p value and FDR.
5) prefix_output_net.txt
Transcript cluster file. Transcript pairs with FDR less than 0.1% will be selected to construct the network.

### 1.5	Tools for visualizing results
- The relative SNV density (RSD) in each 10kb window in all chromosomes
Example Usage:
`Plot --input <inputfile> --indexfile <index file> --strain <strain name> --output <output file name>`

- The average linkage disequilibrium (LD) between SNV pairs  
Example Usage:
`Rscript plotLD.R <first inputfile> <second inputfile> <outputname>`
`<first inputfile>`
linkage disequilibrium (LD) values of genomic background
`<second inputfile>`
linkage disequilibrium (LD) values of PASS regions
`<outputname>`
output name
- Distribution of Tajima’s D in 10 kb windows
Example Usage:
`Rscript plotTaj.R <first inputfile> <second inputfile> <outputname>`
`<first inputfile>` 
Tajima’s D value of genomic background
`<second inputfile>`
Tajima’s D value of PASS regions
`<outputname>`
output name
