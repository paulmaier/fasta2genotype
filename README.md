# fasta2genotype

This program takes a fasta file listing all sequence haplotypes of all individuals at all loci as well as a list of individuals/populations and list of white loci then outputs data in one of eight formats:

(1) migrate-n, (2) Arlequin, (3) DIYabc, (4) LFMM, (5) Phylip, (6) G-Phocs, or (7) Treemix
(8) Additionally, the data can be coded as unique sequence integers (haplotypes) in Structure, Genepop, SamBada, Bayescan, Arlequin, GenAlEx format, or summarized as allele frequencies by population.

### Execute program in the following way:

`python fasta2genotype.py [fasta file] [whitelist file] [population file] [VCF file] [output name]`

### Quality filtering options:
* In addition to coverage filtering, several other quality control measures can be selected.
* Monomorphic loci can be removed.
* Loci suspected of being paralogs (assembled as one ‘Frankenstein’ locus, but including alleles from different loci) can be removed. Loci surpassing the desired threshold value of heterozygosity are tested for Hardy-Weinberg genotype proportions. If they fail (based on a chi-squared test), and have higher heterozygosity than expected under Hardy-Weinberg, they are removed.
* Alleles under a locus-wide frequency can be removed
* Alleles represented under a given frequency of populations can be removed (e.g. 4 pops of 16, frequency of 0.25)
* Loci can be removed from populations where those populations fall under a missing data threshold for those loci
* Loci under a given locus-wide missing data threshold can be removed.
* Individuals missing a threshold frequency of loci can be removed.
* Restriction enzyme sites can be removed if requested, for single- or double-digest setups. Simply select this option and provide the 5' and/or 3’ sequence(s). There can be multiple sequences for either 5’ or 3’ end. For example, you might have double-digest RAD data with READ 1 and READ 2 in the same fasta file. You can tell the program to remove either ‘TGCAGG’ or ‘CGG’ depending on which is found on the 5’ end of a given sequence. Regardless of restriction enzymes or adapters used, it might help to examine the fasta file before choosing this option. Your sequences might be reversed or reverse-complemented.
* If you are exporting an Phylip alignment file, there are several options:
1. SNPs or sequences: concatenate just the polymorphic sites (SNPs), or full sequences, including invariable sites.
2. The alignment can be summarized at the level of haplotypes, individuals, or populations, using IUPAC ambiguity codes.
3. You can select a subset of loci that are “phylogenetically informative” (PI), meaning fixed for alternate alleles at 2+ taxa, or simply “fixed” (at 1+ taxa).
4. If you choose the option for “PI”/”fixed” loci, you will probably want to output SNPs rather than full sequences, since you probably want only PI or fixed sites. However, if you choose “PI”/”fixed” loci and “full sequences,” then the program will output complete sequences containing at least 1 PI or fixed SNP.
5. The alignment can be separated by locus with “!” symbols.
6. A header of tab-delimited locus names can be added.

### Citation:

#### Please cite in the following way:
Maier P.A., Vandergast A.G., Ostoja S.M., Aguilar A., Bohonak A.J. (2019). Pleistocene glacial cycles drove lineage diversification and fusion in the Yosemite toad (*Anaxyrus canorus*). Evolution, in press. https://www.doi.org/10.1111/evo.13868
