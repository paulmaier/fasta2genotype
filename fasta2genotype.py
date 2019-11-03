#!/usr/bin/env python
#
# (c) Paul Maier
# July 21, 2017
# fasta2genotype.py
# V 1.10
# Written for Python 2.7.10
#
# This program takes a fasta file listing sequence haplotypes of all individuals at all loci
# as well as a list of individuals/populations and list of variable loci then outputs data in 
# one of eight formats:
# (1) migrate-n, (2) Arlequin, (3) DIYabc, (4) LFMM, (5) Phylip, (6) G-Phocs, or (8) Treemix
# (8) Additionally, the data can be coded as unique sequence integers (haplotypes)
#     in Structure/Genepop/SamBada/Bayescan/Arlequin/GenAlEx format
#     or summarized as allele frequencies by population
#
# Execute program in the following way:
# python fasta2genotype.py [fasta file] [whitelist file] [population file] [VCF file] [output name]
#
#

import sys,re,csv,collections,itertools
from decimal import *
import numpy as np
from scipy import stats


print """
###################################################################
###                                                             ###
###       Fasta2Genotype | Data Conversion | Version 1.10       ###
###                                                             ###
###                        Cite as follows:                     ###
###                                                             ###
###   Maier P.A., Vandergast A.G., Ostoja S.M., Aguilar A.,     ###
###   Bohonak A.J. (2019). Pleistocene glacial cycles drove     ###
###   lineage diversification and fusion in the Yosemite toad   ###
###   (Anaxyrus canorus). Evolution, in press.                  ###
###   https://www.doi.org/10.1111/evo.13868                     ###
###                                                             ###
###################################################################
"""


if len(sys.argv) != 6:
        print "     ** Error: improper number of arguments. Please see manual for instructions. **"
        print "fasta2genotype.py [fasta file] [whitelist file] [population file] [VCF file] [output name]"
        exit(1)

outname = str(sys.argv[5])
outfile = outname + ".out"
outfile_loci = outname + "_loci.out"
outfile_pops = outname + "_pops.out"

while True:
        try: choice = int(raw_input("Output type? [1] Migrate [2] Arlequin [3] DIYABC [4] LFMM [5] Phylip [6] G-Phocs [7] Treemix [8] Haplotype: "))
        except ValueError: print "     ** Warning: Not a valid option. **"
        else:
                if 1 <= choice < 9: break
                else: print "     ** Warning: Not a valid option. **"

if choice == 2 or choice == 3:
        title = raw_input("Title of project? : ")

if choice == 5:
        while True:
                try: haplo = int(raw_input("Use SNPs or full sequences for alignment? [1] SNPs [2] Full Sequences : "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= haplo < 3: break
                        else: print "     ** Warning: Not a valid option. **"

if choice == 5:
        while True:
                try: haplotypes = int(raw_input("Type of sequences for alignment? [1] Haploid [2] Diploid [3] Population: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= haplotypes < 4: break
                        else: print "     ** Warning: Not a valid option. **"

if choice == 5:
        while True:
                try: phylo_inform = int(raw_input("Keep only phylogenetically informative (PI) loci, fixed loci, or all loci? [1] PI [2] Fixed [3] All: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= phylo_inform < 4: break
                        else: print "     ** Warning: Not a valid option. **"

if choice == 5:
        while True:
                try: breakpoints = int(raw_input("Flag break points between loci with '!' symbol? [1] Yes [2] No : "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= breakpoints < 3: break
                        else: print "     ** Warning: Not a valid option. **"

if choice == 5:
        while True:
                try: locheaders = int(raw_input("Insert locus name headers in first row? [1] Yes [2] No : "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= locheaders < 3: break
                        else: print "     ** Warning: Not a valid option. **"

if choice == 7:
        while True:
                try: one_snp = int(raw_input("How many SNPs to keep per locus? [1] Only one [2] All : "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= one_snp < 3: break
                        else: print "     ** Warning: Not a valid option. **"

if choice == 8:
        title = ""
        FourOrSix = 0
        while True:
                try: HaploChoice = int(raw_input("Specific output type? [1] Structure [2] Genepop [3] AlleleFreqency [4] SamBada [5] Bayescan [6] Arlequin [7] GenAlEx : "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= HaploChoice < 8: break
                        else: print "     ** Warning: Not a valid option. **"
        if HaploChoice == 2:
                while True:
                        try: FourOrSix = int(raw_input("Genepop in four [1] or six [2] digit format? "))
                        except ValueError: print "     ** Warning: Not a valid option. **"
                        else:
                                if 1 <= FourOrSix < 3: break
                                else: print "     ** Warning: Not a valid option. **"
        if HaploChoice == 2 or HaploChoice == 6:
                title = raw_input("Title of project? : ")

outtype = 2
if sys.argv[2] != 'NA' and sys.argv[2] != 'na' and sys.argv[2] != 'Na' and sys.argv[2] != 'nA':
        while True:
                try: outtype = int(raw_input("Loci to use? [1] Whitelist [2] All: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 1 <= outtype < 3: break
                        else: print "     ** Warning: Not a valid option. **"

while True:
        try: clipcutsite = int(raw_input("Remove restriction enzyme or adapter sequences? These may bias data. [1] Yes [2] No: "))
        except ValueError: print "     ** Warning: Not a valid option. **"
        else:
                if 1 <= clipcutsite < 3: break
                else: print "     ** Warning: Not a valid option. **"

cutsite1 = ""
cutsite2 = ""
if clipcutsite == 1:
        while True:
                cutsite1 = raw_input("Beginning (5') sequence(s) to remove? (If multiple use spaces, if none leave blank): ")
                if re.match("^[ATCGatcg ]*$", cutsite1): break
                else: print "     ** Warning: Not a valid option. **"
        while True:
                cutsite2 = raw_input("Ending (3') sequence(s) to remove? (If multiple use spaces, if none leave blank): ")
                if re.match("^[ATCGatcg ]*$", cutsite2): break
                else: print "     ** Warning: Not a valid option. **"
	cutsite1 = cutsite1.upper()
	cutsite2 = cutsite2.upper()	
	cutsite1 = cutsite1.split()
	cutsite2 = cutsite2.split()

UseCoverage = 0
CoverageCutoff = 0
monomorphic_filter2 = 2
if sys.argv[4] != 'NA' and sys.argv[4] != 'na' and sys.argv[4] != 'Na' and sys.argv[4] != 'nA':
        UseCoverage = 1
        while True:
                try: CoverageCutoff = int(raw_input("Coverage Cutoff (number reads for locus)? Use '0' to ignore coverage: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 0 < CoverageCutoff: monomorphic_filter2 = 1
                        if 0 <= CoverageCutoff: break
                        else: print "     ** Warning: Not a valid option. **"

while True:
        try: monomorphic_filter = int(raw_input("Remove monomorphic loci? [1] Yes [2] No: "))
        except ValueError: print "     ** Warning: Not a valid option. **"
        else:
                if 1 <= monomorphic_filter < 3: break
                else: print "     ** Warning: Not a valid option. **"

heterocutoff = 0
while True:
        try: hetero_filter = int(raw_input("Remove loci with excess heterozygosity? This can remove paralogs. [1] Yes [2] No: "))
        except ValueError: print "     ** Warning: Not a valid option. **"
        else:
                if 1 <= hetero_filter < 3: break
                else: print "     ** Warning: Not a valid option. **"

if hetero_filter == 1:
        while True:
                try: heterocutoff = float(raw_input("Maximum heterozygosity cutoff for removing loci out of Hardy-Weinberg? "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 0 <= heterocutoff <= 1: break
                        else: print "     ** Warning: Not a valid option. **"

while True:
        try: allele_filter = int(raw_input("Filter for allele frequency? False alleles might bias data. [1] Yes [2] No: "))
        except ValueError: print "     ** Warning: Not a valid option. **"
        else:
                if 1 <= allele_filter < 3: break
                else: print "     ** Warning: Not a valid option. **"

allele_threshold = 0
allele_pop_threshold = 0
if allele_filter == 1:
        while True:
                try: allele_threshold = float(raw_input("Allele frequency threshold for removal across all individuals? Use '0' to ignore this: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 0 <= allele_threshold <= 1: break
                        else: print "     ** Warning: Not a valid option. **"
        while True:
                try: allele_pop_threshold = float(raw_input("Frequency of populations containing allele for removal across all individuals? Use '0' to ignore this: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 0 <= allele_pop_threshold <= 1: break
                        else: print "     ** Warning: Not a valid option. **"

while True:
        try: missing_data_filter = int(raw_input("Filter for missing genotypes? These might bias data. [1] Yes [2] No: "))
        except ValueError: print "     ** Warning: Not a valid option. **"
        else:
                if 1 <= missing_data_filter < 3: break
                else: print "     ** Warning: Not a valid option. **"

if missing_data_filter == 1:
        while True:
                try: locus_threshold = float(raw_input("Locus frequency threshold for locus removal across all individuals? Use '0' to ignore this: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 0 <= locus_threshold <= 1: break
                        else: print "     ** Warning: Not a valid option. **"
        while True:
                try: locus_pop_threshold = float(raw_input("Population frequency threshold for locus removal across each population? Use '0' to ignore this: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 0 <= locus_pop_threshold <= 1: break
                        else: print "     ** Warning: Not a valid option. **"
        while True:
                try: ind_threshold = float(raw_input("Individual frequency threshold for individual removal across all loci? Use '0' to ignore this: "))
                except ValueError: print "     ** Warning: Not a valid option. **"
                else:
                        if 0 <= ind_threshold <= 1: break
                        else: print "     ** Warning: Not a valid option. **"

print " "
print "**************************************************************************************************************"
print "***                                       ... BEGINNING CONVERSION ...                                     ***"
print "**************************************************************************************************************"
print " "




# Create dictionary of populations and individuals
def Pops(populations, seqsdict):
        print "Cataloging populations..."
        try:
		pops = csv.DictReader(open(populations,"U"), delimiter="\t", quotechar='"', dialect="excel-tab")
        except IOError:
                print "Error: File does not appear to exist. Check the file and the directory path."
                exit(1)

        popsdict = {} #Structure: {Population : {SampleID : IndividualID} }
        for i in pops:
                if i[pops.fieldnames[2]] in popsdict.keys():
                        popsdict[i[pops.fieldnames[2]]][i[pops.fieldnames[0]]] = i[pops.fieldnames[1]]
		else:
			popsdict[i[pops.fieldnames[2]]] = {i[pops.fieldnames[0]]:i[pops.fieldnames[1]]}

	inds = []
	for k in popsdict.iterkeys():
                inds.extend(popsdict[k].keys())
        inds2 = seqsdict.keys()
        diff = np.setdiff1d(inds, inds2)

        for i in diff:
                if i in seqsdict.keys(): del seqsdict[i]
                for k in popsdict.iterkeys():
                        if i in popsdict[k].keys(): del popsdict[k][i]
                
        num_pops = len(popsdict)
        
	print "Counting gene copies..."
        gene_copies = {} #Structure: {PopulationID : DiploidGeneCopies}
        for i, j in popsdict.iteritems():
                gene_copies[i] = 2*len(j)
                
        return popsdict, num_pops, gene_copies, 




# Create dictionary of Coverage by individual X locus
def LocusCoverage(VCFfile):
        print "Calculating loci coverage..."
        try:
                with open(VCFfile,"U") as f:
                        cov = csv.reader(f, delimiter="\t")
                        d = list(cov)

        except IOError:
                print "Error: File does not appear to exist. Check the file and the directory path."
                exit(1)
        covdict = {} # Structure: {IndividualID : {LocusID : Coverage} }
        rindex = 0
        for row in d:
                if rindex < 9:
                        rindex += 1
                        continue
                cindex = 0
                locus = int(d[rindex][2])
                for column in row:
                        ind = d[8][cindex]
                        if cindex < 9:
                                cindex += 1
                                continue
                        if ind in covdict.keys():
                                covnum = int(re.sub(r'\S+:(\d+):\S+', r'\1', d[rindex][cindex]))
                                covdict[ind][locus] = covnum
                        else:
                                covnum = int(re.sub(r'\S+:(\d+):\S+', r'\1', d[rindex][cindex]))
                                covdict[ind] = {locus : covnum}
                        cindex += 1
                rindex += 1
        return covdict




# Create dictionary of individuals, sequences, and alleles
def Seqs (outtype, clipcutsite, cutsite1, cutsite2, CoverageCutoff, covdict):
        print "Cataloging loci..."

        whitelist = []
        if outtype == 1:
                try:
                        fin=open(sys.argv[2],"U")
                        for line in fin:
                                whitelist.append(str(line.strip()))
                        fin.close()
                except IOError:
                        print "Error: File does not appear to exist. Check the file and the directory path."
                        exit(1)

        if CoverageCutoff > 0:
                printtext = "     Removing genotypes below coverage threshold of %s..."
                printval = (CoverageCutoff)
                print (printtext % printval)
        newfasta = open(sys.argv[1],"U") #This is original fasta if outtype == 2
        seqsdict = {} #Structure: {SampleID : {LocusID : {AlleleID : DNAsequence} } }
        carrot = ">"

        try: # Make temporary dictionary to associate SampleID and IndividualID using population file
                # This is used for coverage cutoff option
                pops = csv.DictReader(open(sys.argv[3],"U"), delimiter="\t", quotechar='"', dialect="excel-tab")
        except IOError:
                print "     ** Error: File does not appear to exist. Check the file and the directory path. **"
                exit(1)
        popsdicttemp = {}
        for i in pops:
                popsdicttemp[i[pops.fieldnames[0]]] = i[pops.fieldnames[1]]

        if clipcutsite == 1 and cutsite1 != []:
                cliplist1 = ', '.join(cutsite1)
                if len(cutsite1) == 1: print ("     Clipping sequence %s from left side..." %cliplist1)
                if len(cutsite1) > 1: print ("     Clipping whichever sequence %s is found on left side..." %cliplist1)
        if clipcutsite == 1 and cutsite2 != []:
                cliplist2 = ', '.join(cutsite2)
                if len(cutsite2) == 1: print ("     Clipping sequence %s from right side..." %cliplist2)
                if len(cutsite2) > 1: print ("     Clipping whichever sequence %s is found on right side..." %cliplist2)

        duplicated_alleles = 0 # Count instances of duplicated gene copies in homozygotes
        three_alleles = 0 # Count instances of three or more alleles in one individual

        print "Counting locus lengths..."
        num_sites = {} #Structure: {LocusID : NumberNucleotides}
        
        for line in newfasta:
                if carrot in line:
                        indnum = re.sub(r'>CLocus_\w+_Sample_(\w+)_Locus_\w+_Allele_\w+', r'\1', line); indnum=indnum.strip()
                        locusnum = re.sub(r'>CLocus_(\w+)_Sample_\w+_Locus_\w+_Allele_\w+', r'\1', line); locusnum=locusnum.strip()
                        allelenum = re.sub(r'>CLocus_\w+_Sample_\w+_Locus_\w+_Allele_(\w+)', r'\1', line); allelenum=allelenum.strip()
                        nextline = newfasta.next(); nextline = nextline.strip()

                        if locusnum in whitelist or len(whitelist) == 0:

                                # Clip off cut sites
                                if clipcutsite == 1:
                                        for i in range(0,len(cutsite1)):
                                                if nextline[0:len(cutsite1[i])] == cutsite1[i]:
                                                        nextline = nextline[len(cutsite1[i]):]
                                                        break
                                        for i in range(0,len(cutsite2)):
                                                if nextline[(len(nextline)-len(cutsite2[i])):] == cutsite2[i]:
                                                        nextline = nextline[:(len(nextline)-len(cutsite2[i]))]
                                                        break

                                if locusnum not in num_sites: num_sites[locusnum] = len(nextline)
                                
                                # Produce seqsdict without coverage cutoff
                                if CoverageCutoff == 0:
                                        if indnum in seqsdict.keys():
                                                if locusnum in seqsdict[indnum].keys():
                                                        if allelenum not in seqsdict[indnum][locusnum].keys():
                                                                if int(allelenum) in range (0,2): seqsdict[indnum][locusnum][allelenum] = nextline
                                                                else: three_alleles += 1
                                                        else:
                                                                duplicated_alleles += 1
                                                else:
                                                        if int(allelenum) in range (0,2): seqsdict[indnum][locusnum] = {allelenum:nextline}
                                                        else: three_alleles += 1
                                        else:
                                                if int(allelenum) in range (0,2): seqsdict[indnum] = {locusnum:{allelenum:nextline}}
                                                else: three_alleles += 1
                                                
                                # Produce seqsdict with coverage cutoff
                                if CoverageCutoff > 0:
                                        if indnum in popsdicttemp.keys():
                                                if int(locusnum) in covdict[popsdicttemp[indnum]].keys():
                                                        if indnum in seqsdict.keys():
                                                                if locusnum in seqsdict[indnum].keys():
                                                                        if allelenum not in seqsdict[indnum][locusnum].keys() and covdict[popsdicttemp[indnum]][int(locusnum)] >= CoverageCutoff:
                                                                                if int(allelenum) in range (0,2): seqsdict[indnum][locusnum][allelenum] = nextline
                                                                                else: three_alleles += 1
                                                                        elif allelenum not in seqsdict[indnum][locusnum].keys() and covdict[popsdicttemp[indnum]][int(locusnum)] < CoverageCutoff:
                                                                                pass
                                                                        else:
                                                                                duplicated_alleles += 1
                                                                else:
                                                                        if covdict[popsdicttemp[indnum]][int(locusnum)] >= CoverageCutoff:
                                                                                if int(allelenum) in range (0,2): seqsdict[indnum][locusnum] = {allelenum:nextline}
                                                                                else: three_alleles += 1
                                                        else:
                                                                if covdict[popsdicttemp[indnum]][int(locusnum)] >= CoverageCutoff:
                                                                        if int(allelenum) in range (0,2): seqsdict[indnum] = {locusnum:{allelenum:nextline}}
                                                                        else: three_alleles += 1

        if duplicated_alleles > 0: print ("     ** Warning: %s homozyogotes had both gene copies in fasta file. Removing duplicate sequence(s). **" %duplicated_alleles)
        if three_alleles > 0: print ("     ** Warning: %s genotypes had 3 or more alleles in fasta file. Keeping only alleles '0' and '1'. **" %three_alleles)

	newfasta.close()
        return seqsdict, num_sites




# Screen for monomorphic loci and excess heterozygosity
def LocusRemoval(seqsdict, popsdict, gene_copies, num_sites, monomorphic_filter, hetero_filter, heterocutoff):
        if monomorphic_filter2 == 1: monomorphic_filter = 1
        #Create dictionary of loci containing dictionaries of allele counts
        allelecount = {} #Build structure: {LocusID : {Sequence : Count} }
        print "Counting alleles for each locus..."
        for x in sorted(seqsdict.iterkeys()): #Cycle through all individuals
                for p in sorted(seqsdict[x].iterkeys()): # Cycle through all loci
                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through both alleles
                                if p in allelecount.keys():
                                        if str(seqsdict[x][p][a]) not in allelecount[p].keys(): #If allele not added yet, add and set count to 1
                                                
                                                if len(seqsdict[x][p]) == 1:
                                                        allelecount[p][str(seqsdict[x][p][a])] = 2
                                                else:
                                                        allelecount[p][str(seqsdict[x][p][a])] = 1
                                        else:
                                                if len(seqsdict[x][p]) == 1:
                                                        allelecount[p][str(seqsdict[x][p][a])] += 2
                                                else:
                                                        allelecount[p][str(seqsdict[x][p][a])] += 1
                                else:
                                        if len(seqsdict[x][p]) == 1:
                                                allelecount[p] = {str(seqsdict[x][p][a]):2}
                                        else:
                                                allelecount[p] = {str(seqsdict[x][p][a]):1}

        if hetero_filter == 1:
        		#Create dictionary of loci containing dictionaries of observed heterozygote and homozygote counts
        		genocount = {} #Build structure: {LocusID : {Genotype : Count} }
        		print "Identifying loci with excess heterozygosity..."
        		print "     Calculating observed heterozygosity and homozygosity..."
        		for x in sorted(seqsdict.iterkeys()): #Cycle through all individuals
                		for p in sorted(seqsdict[x].iterkeys()): # Cycle through all loci
                                        #if len([k for k,v in seqsdict[x][p].items() if list(seqsdict[x][p].values()).count(v)==1]) > 1: #More than 1 unique allele
                                        if int(len(seqsdict[x][p])) >= 2: # Heterozygote
                                                a1 = min(seqsdict[x][p]["0"], seqsdict[x][p]["1"])
						a2 = max(seqsdict[x][p]["0"], seqsdict[x][p]["1"])
                                                if p in genocount.keys():
                                        		if str(a1+"/"+a2) not in genocount[p].keys(): #If genotype not added yet, add and set count to 1
                                                		genocount[p][str(a1+"/"+a2)] = 1
                                        		else:
                                                		genocount[p][str(a1+"/"+a2)] += 1
                                		else:
                                        		genocount[p] = {str(a1+"/"+a2):1}
                                        elif int(len(seqsdict[x][p])) == 1: # Homozygote
                                                a0 = seqsdict[x][p].values()[0]
                                                if p in genocount.keys():
                                        		if str(a0+"/"+a0) not in genocount[p].keys(): #If genotype not added yet, add and set count to 1
                                                		genocount[p][str(a0+"/"+a0)] = 1
                                        		else:
                                                		genocount[p][str(a0+"/"+a0)] += 1
                                		else:
                                        		genocount[p] = {str(a0+"/"+a0):1}

        		#Create dictionary of loci containing dictionaries of expected heterozygote and homozygote counts
        		genoexpect = {} #Build structure: {LocusID : {Genotype : Count} }
        		print "     Calculating expected heterozygosity and homozygosity..."
                	for p in sorted(allelecount.iterkeys()): # Cycle through all loci
				for a in sorted(allelecount[p].iterkeys()):
                                        for b in sorted(allelecount[p].iterkeys()):
						if a < b:
                                                        a1 = min(a,b); a2 = max(a,b)
                                                        exp_het = 2 * allelecount[p][a1]/float(sum(allelecount[p].values())) * allelecount[p][a2]/float(sum(allelecount[p].values())) * sum(allelecount[p].values())/2.0
                                                        if p in genoexpect.keys():
                                                                genoexpect[p][str(a1+"/"+a2)] = exp_het
                                                        else:
                                                                genoexpect[p] = { str(a1+"/"+a2) : exp_het }
                                                elif a == b:
                                                        exp_hom = ((allelecount[p][a]/float(sum(allelecount[p].values())))**2) * sum(allelecount[p].values())/2.0
                                                        if p in genoexpect.keys():
                                                                genoexpect[p][str(a+"/"+a)] = exp_hom
                                                        else:
                                                                genoexpect[p] = { str(a+"/"+a) : exp_hom }
                                
                        #Create dictionary of loci to remove based on excess heterozygosity (1 = remove, 0 = keep)
                        remove_hetero = {} #Build structure: {LocusID : Remove_value }
                        print "     Flagging loci with excess heterozygosity for removal..."
                        for p in sorted(genoexpect.iterkeys()):
                                remove_hetero[p] = 0
                                high_hetero = 0
                                for a in sorted(genocount[p].iterkeys()):
                                        if genocount[p][a]/(float(sum(allelecount[p].values()))/2.0) >= heterocutoff and re.sub(r'(\w+)/\w+', r'\1', a) != re.sub(r'\w+/(\w+)', r'\1', a):
                                                if genocount[p][a] > genoexpect[p][a]:
                                                        high_hetero = 1
                                if high_hetero == 1:
                                        chi = 0
                                        num_genos = len(allelecount[p]) * (len(allelecount[p]) - 1) / 2 + len(allelecount[p])
                                        for a in sorted(genoexpect[p].iterkeys()):
                                                try:
                                                        chi += (genoexpect[p][a] - genocount[p][a])**2/float(genoexpect[p][a])
                                                except KeyError:
                                                        chi += (genoexpect[p][a] - 0)**2/float(genoexpect[p][a])

                                        df = num_genos - len(allelecount[p])

                                        if (1 - stats.chi2.cdf(chi, df)) < 0.05:
                                                remove_hetero[p] = 1

        if monomorphic_filter == 1 or hetero_filter == 1:
                print "     Removing loci..."
                removed_mono = 0
                removed_hetero = 0

                # If coverage filtering is used, remove loci from num_sites that don't exist in seqsdict
                keep_loci = []
                for i in seqsdict.iterkeys():
                        for j in seqsdict[i].iterkeys():
                                if j not in keep_loci: keep_loci.append(j)
                rem = np.setdiff1d(num_sites.keys(), keep_loci)
                for i in rem: del num_sites[i]
                
                for p in sorted(allelecount.iterkeys()):
                        if len(allelecount[p]) == 1 and monomorphic_filter == 1:
                                removed_mono += 1
                                for x in sorted(seqsdict.iterkeys()):
                                        if p in seqsdict[x]:
                                                del seqsdict[x][p]
                                                if p in num_sites:
                                                        del num_sites[p]
                        if hetero_filter == 1 and remove_hetero[p] == 1:
                                try:
                                        removed_hetero += 1
                                        for x in sorted(seqsdict.iterkeys()):
                                                if p in seqsdict[x]:
                                                        del seqsdict[x][p]
                                                        if p in num_sites:
                                                                del num_sites[p]
                                except KeyError:
                                        continue
        if removed_mono > 0: print ("     Removed %s monomorphic loci." % removed_mono)
        if removed_hetero > 0: print ("     Removed %s overly heterozygous loci." % removed_hetero)

        return seqsdict, num_sites




# Screen for false alleles
def AlleleRemoval(seqsdict, popsdict, gene_copies, num_sites, allele_threshold, allele_pop_threshold, allele_filter):
        #Create dictionary of loci containing dictionaries of allele counts
        allelecount = {} #Build structure: {LocusID : {Sequence : Count} }
        print "Counting alleles..."
        for x in sorted(seqsdict.iterkeys()): #Cycle through all individuals
                for p in sorted(seqsdict[x].iterkeys()): # Cycle through all loci
                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through both alleles
                                if p in allelecount.keys():
                                        if str(seqsdict[x][p][a]) not in allelecount[p].keys(): #If allele not added yet, add and set count to 1
                                                allelecount[p][str(seqsdict[x][p][a])] = 1
                                        else:
                                                allelecount[p][str(seqsdict[x][p][a])] += 1
                                else:
                                        allelecount[p] = {str(seqsdict[x][p][a]):1}

        if allele_filter == 1:
                print "     Removing alleles using thresholds..."
                if allele_threshold != 0: printtext = "     Allele must have an overall frequency of %s..."; printval = allele_threshold; print (printtext % printval)
                if allele_pop_threshold != 0: printtext = "     Allele must be present in %s of populations..."; printval = allele_pop_threshold; print (printtext % printval)
                 
                num_inds = 0
                for k in popsdict.iterkeys():
                        num_inds += int(len(popsdict[k]))

                # Locus by locus in allelecount dictionary
                removed_alleles = 0
                for p in sorted(allelecount.iterkeys()):
                        count = 0
                        remove_flag = 0
                        for b in sorted(allelecount[p].iterkeys()):
                                count = allelecount[p][b]
                                if float(count)/(2.0*float(num_inds)) < allele_threshold and allele_threshold != 0: #Remove alleles below threshold
                                        remove_flag = 1
                                        removed_alleles += 1
                                        for x in sorted(seqsdict.iterkeys()):
                                                if p in seqsdict[x].keys():
                                                        for a in sorted(seqsdict[x][p].iterkeys()):
                                                        		if p in seqsdict[x]: #This locus might have been removed for individual already
                                                                		if seqsdict[x][p][a] == b:
                                                                				del seqsdict[x][p]
                                                                				if b in allelecount[p]: del allelecount[p][b]
                                if remove_flag == 0: #Remove alleles not present in enough populations
                                        count = 0
                                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                                flag = 0
                                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                                if p in seqsdict[x].keys():
                                                                        for a in sorted(seqsdict[x][p].itervalues()):
                                                                                if a == b:
                                                                                        flag = 1
                                                if flag == 1: count += 1
                                        if float(count)/(float(gene_copies[k])) < allele_pop_threshold and allele_pop_threshold != 0:
                                                removed_alleles += 1
                                                for x in sorted(seqsdict.iterkeys()):
                                                        if p in seqsdict[x].keys():
                                                                for a in sorted(seqsdict[x][p].iterkeys()):
                                                        				if p in seqsdict[x]: #This locus might have been removed for individual already
                                                                				if seqsdict[x][p][a] == b:
                                                                						del seqsdict[x][p]
                                                                						if b in allelecount[p]: del allelecount[p][b]
                                                                						
        if removed_alleles > 0: print ("     Removed %s alleles below threshold." % removed_alleles)
        
        return seqsdict, num_sites




# Screen for missing data
def MissingData(seqsdict, popsdict, gene_copies, num_sites, locus_threshold, locus_pop_threshold, ind_threshold):
        print "Applying missing data thresholds..."

        if locus_pop_threshold != 0: printtext = "     Locus must have a frequency of %s in each population..."; printval = locus_pop_threshold; print (printtext % printval)
        if locus_pop_threshold != 0: # Remove loci from pops when below pop threshold
                locus_pop_count = {}
                removed_popmissing = 0
                for k in sorted(popsdict.iterkeys()): #Look at one population
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                if p in locus_pop_count.keys():
                                                        if k in locus_pop_count[p].keys(): locus_pop_count[p][k] += 1
                                                        else: locus_pop_count[p][k] = 1
                                                else:
                                                        locus_pop_count[p] = {k : 1}
                for p in locus_pop_count.iterkeys():
                        for k in locus_pop_count[p].iterkeys():
                                if float(locus_pop_count[p][k])/(float(gene_copies[k])/2.0) < locus_pop_threshold:
                                        removed_popmissing += 1
                                        for x in sorted(seqsdict.iterkeys()):
                                                if p in seqsdict[x] and x in popsdict[k].keys():
                                                        del seqsdict[x][p]
                if removed_popmissing > 0: print ("     Removed %s loci below threshold." % removed_popmissing)

        if locus_threshold != 0: printtext = "     Locus must have an overall frequency of %s..."; printval = locus_threshold; print (printtext % printval)                                                       
        if locus_threshold != 0: # Remove loci below threshold                                                       
                num_inds = 0
                for k in popsdict.iterkeys():
                        num_inds += int(len(popsdict[k]))

                locus_count = {} #Structure: { LocusID : count }
                removed_missing = 0
                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                if p in locus_count.keys(): locus_count[p] += 1
                                else: locus_count[p] = 1
                for p in locus_count.iterkeys():
                        if float(locus_count[p])/(float(num_inds)) < locus_threshold:
                                removed_missing += 1
                                for x in sorted(seqsdict.iterkeys()):
                                        if p in seqsdict[x]:
                                                del seqsdict[x][p]
                                                if p in num_sites:
                                                        del num_sites[p]
                if removed_missing > 0: print ("     Removed %s loci below threshold." % removed_missing)

        if ind_threshold != 0: printtext = "     Individual must have %s of total loci..."; printval = ind_threshold; print (printtext % printval)                                                       
        if ind_threshold != 0: # Remove individuals below threshold
                locus_list = []
                removed_inds = 0
                for x in sorted(seqsdict.iterkeys()):
                        for p in sorted(seqsdict[x].iterkeys()):
                                if p not in locus_list: locus_list.append(p)
                new_num_loci = len(locus_list)
                if new_num_loci == 0:
                        print "     ** All loci removed. Check data. **"
                        exit(1)

                locus_ind_count = {}
                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                        locus_ind_count[x] = len(seqsdict[x])

                for x in sorted(seqsdict.iterkeys()):
                        if float(locus_ind_count[x])/float(new_num_loci) < ind_threshold:
                                remove_inds += 1
                                del seqsdict[x]
                if removed_inds > 0: print ("     Removed %s individuals below threshold." % removed_inds)

        return seqsdict, num_sites




def IUPAC(letters):
        code = letters[0]
        for letter in letters:
                if letter in ['A','T','C','G'] and code in ['A','T','C','G','M','R','W','S','Y','K','V','H','D','B']:
                
                        if letter == "C":
                                if code == "C":
                                        code = "C"                       
                                if code == "G":
                                        code = "S"
                                if code == "A":
                                        code = "M"
                                if code == "T":
                                        code = "Y"
                                if code == "M":
                                        code = "M"                       
                                if code == "R":
                                        code = "V"
                                if code == "W":
                                        code = "H"
                                if code == "S":
                                        code = "S"
                                if code == "Y":
                                        code = "Y"                       
                                if code == "K":
                                        code = "B"
                                if code == "V":
                                        code = "V"
                                if code == "H":
                                        code = "H"
                                if code == "D":
                                        code = "N"
                                if code == "B":
                                        code = "B"
                                        
                                
                        if letter == "G":
                                if code == "G":
                                        code = "G"
                                if code == "C":
                                        code = "S"
                                if code == "A":
                                        code = "R"
                                if code == "T":
                                        code = "K"
                                if code == "M":
                                        code = "V"                       
                                if code == "R":
                                        code = "R"
                                if code == "W":
                                        code = "D"
                                if code == "S":
                                        code = "S"
                                if code == "Y":
                                        code = "B"                       
                                if code == "K":
                                        code = "K"
                                if code == "V":
                                        code = "V"
                                if code == "H":
                                        code = "N"
                                if code == "D":
                                        code = "D"
                                if code == "B":
                                        code = "B"

                        if letter == "A":
                                if code == "A":
                                        code = "A"
                                if code == "C":
                                        code = "M"
                                if code == "G":
                                        code = "R"
                                if code == "T":
                                        code = "W"
                                if code == "M":
                                        code = "M"                       
                                if code == "R":
                                        code = "R"
                                if code == "W":
                                        code = "W"
                                if code == "S":
                                        code = "V"
                                if code == "Y":
                                        code = "H"                       
                                if code == "K":
                                        code = "D"
                                if code == "V":
                                        code = "V"
                                if code == "H":
                                        code = "H"
                                if code == "D":
                                        code = "D"
                                if code == "B":
                                        code = "N"

                        if letter == "T":
                                if code == "T":
                                        code = "T"
                                if code == "C":
                                        code = "Y"
                                if code == "G":
                                        code = "K"
                                if code == "A":
                                        code = "W"
                                if code == "M":
                                        code = "H"                       
                                if code == "R":
                                        code = "D"
                                if code == "W":
                                        code = "W"
                                if code == "S":
                                        code = "B"
                                if code == "Y":
                                        code = "Y"                       
                                if code == "K":
                                        code = "K"
                                if code == "V":
                                        code = "N"
                                if code == "H":
                                        code = "H"
                                if code == "D":
                                        code = "D"
                                if code == "B":
                                        code = "B"

                else:
                        code = "N"
                                
        return code




def IUPAC_fixed(letters):
        fixed = 0
        for letter1 in letters:
                for letter2 in letters:
                        if letter1 < letter2 and letter1 in ['A','T','C','G','M','R','W','S','Y','K','V','H','D','B'] and letter2 in ['A','T','C','G','M','R','W','S','Y','K','V','H','D','B']:
                                if letter1 == "A":
                                        if letter2 in ['T','C','G','S','Y','K','B']:
                                                fixed = 1
                                if letter1 == "T":
                                        if letter2 in ['A','C','G','M','R','S','V']:
                                                fixed = 1
                                if letter1 == "C":
                                        if letter2 in ['A','T','G','R','W','K','D']:
                                                fixed = 1
                                if letter1 == "G":
                                        if letter2 in ['A','T','C','M','W','Y','H']:
                                                fixed = 1
                                if letter1 == "M":
                                        if letter2 in ['T','G','K']:
                                                fixed = 1
                                if letter1 == "R":
                                        if letter2 in ['T','C','Y']:
                                                fixed = 1
                                if letter1 == "W":
                                        if letter2 in ['C','G','S']:
                                                fixed = 1
                                if letter1 == "S":
                                        if letter2 in ['A','T','W']:
                                                fixed = 1
                                if letter1 == "Y":
                                        if letter2 in ['A','G','R']:
                                                fixed = 1
                                if letter1 == "K":
                                        if letter2 in ['A','C','M']:
                                                fixed = 1
                                if letter1 == "V":
                                        if letter2 in ['T']:
                                                fixed = 1
                                if letter1 == "H":
                                        if letter2 in ['G']:
                                                fixed = 1
                                if letter1 == "D":
                                        if letter2 in ['C']:
                                                fixed = 1
                                if letter1 == "B":
                                        if letter2 in ['A']:
                                                fixed = 1
                               
        return fixed




# Output Migrate file
def Fasta2Migrate(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        print "Outputting migrate-n file..."
        try:
                OrderedLoci = []
                par=""
                fout=open(outfile,"w")

                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        
                fout.write(str(num_pops) + '\t' + str(len(OrderedLoci)) + "\n")
                for key in sorted(num_sites.iterkeys()):
                        fout.write("%s\t" % num_sites[key])
                fout.write('\n')
                
                

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        fout.write(str(gene_copies[k]) + '\t' + 'Pop_' + k + '\n')
                        for i in OrderedLoci: #Look at one locus
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                if int(len(str(popsdict[k][x])))>9: print "     ** Error, Ind ID > 9 characters. **"; exit(1)
                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write ?s
                                                        ind = str(popsdict[k][x])
                                                        fout.write((ind+'a').ljust(10)); z=0
                                                        while z < int(num_sites[i]):
                                                                fout.write("?"); z+=1
                                                        fout.write('\n')
                                                        ind = str(popsdict[k][x])
                                                        fout.write((ind+'b').ljust(10)); z=0
                                                        while z < int(num_sites[i]):
                                                                fout.write("?"); z+=1
                                                        fout.write('\n')

                                                else:
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                if p == i:      #If this is the right locus
                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                ind = str(popsdict[k][x])
                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                        if int(a) == 0: par = 'a'
                                                                                        elif int(a) == 1: par = 'b'
                                                                                        fout.write((ind+par).ljust(10)+str(seqsdict[x][p][a])+'\n')
                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                        fout.write((ind+'a').ljust(10)+str(seqsdict[x][p][a])+'\n'+(ind+'b').ljust(10)+str(seqsdict[x][p][a])+'\n')
                print "*** DONE! ***"
                fout.close()

        except IOError:
                print "     ** Error: Problems outputting file. Check the directory path. **"
                return 0

        return 1




# Output Arlequin file
def Fasta2Arlequin(num_pops, popsdict, seqsdict, gene_copies, num_sites, title):
        try:
                OrderedLoci = []
                print "Outputting Arelequin file..."
                fout=open(outfile,"w")

                fout.write("[Profile]\n\n\t\"" + title + "\"\n\n\t\tNbSamples=" + str(num_pops))
                fout.write("\n\t\tGenotypicData=1\n\t\tGameticPhase=0\n\t\tDataType=DNA\n\t\t")
                fout.write("LocusSeparator=TAB\n\t\tMissingData=\"?\"\n\n\n[Data]\n\n\t[[Samples]]\n\n")

                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        fout.write("\t\tSampleName= \"Pop_" + str(k) + "\"\n\t\tSampleSize=" + str(int(gene_copies[k])/2) + "\n\t\tSampleData={\n")
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        count = 0
                                        while count < 2:
                                                ind = str(popsdict[k][x])
                                                if count == 0: fout.write(ind+'\t1\t')
                                                if count == 1: fout.write('\t\t')
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write ?s
                                                                z=0
                                                                while z < int(num_sites[i]):
                                                                        fout.write("?"); z+=1
                                                                fout.write('\t')

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                count += 1
                                                fout.write('\n')
                        fout.write("}\n")
                print "*** DONE! ***"
                fout.close()
        except IOError:
                print "     ** Error: Problems outputting file. Check the directory path. **"
                return 0

        return 1




# Output DIYABC file
def Fasta2DIYABC(num_pops, popsdict, seqsdict, gene_copies, num_sites, title):
        print "Outputting DIYabc file..."
        try:
                OrderedLoci = []
                fout=open(outfile,"w")

                fout.write(title + " <NM=NF>\n")

                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        fout.write(key + "\t<A>\n")

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        fout.write("Pop\n")
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        ind = str(popsdict[k][x])
                                        fout.write(ind+'\t,\t')
                                        for i in OrderedLoci: #Look at one locus
                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write empty brackets
                                                        fout.write("<[][]>\t")

                                                else:
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                if p == i:      #If this is the right locus
                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                        if (int(a)==0):
                                                                                                fout.write("<[" + str(seqsdict[x][p][a]) + "]")
                                                                                        elif (int(a)==1):
                                                                                                fout.write("[" + str(seqsdict[x][p][a]) + "]>\t")

                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                        fout.write("<[" + str(seqsdict[x][p][a]) + "][" + str(seqsdict[x][p][a]) + "]>\t")
                                        fout.write('\n')
                print "*** DONE! ***"
                fout.close()
        except IOError:
                print "     ** Error: Problems outputting file. Check the directory path. **"
                return 0

        return 1




# Output LFMM file
def Fasta2LFMM(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        print "Outputting LFMM (ped) file..."

        try:

                whitelist = []

                with open(sys.argv[2],"U") as f:
                        for line in f:
                                whitelist.append(str(line.strip()))

                with open(sys.argv[4],"U") as f:
                        vcf = csv.reader(f, delimiter="\t")
                        d = list(vcf)
                snpname = str(outname+".snp")
                fout=open(snpname,"w")

                print "     Building dictionary of SNPs..."
                snpdict = {} # Structure: {IndividualID : {SnpID : [Allele1, Allele2] } }
                rindex = 0
                lastlocus = 0
                for row in d:
                        if rindex < 9:
                                rindex += 1
                                continue
                        cindex = 0
                        locus = int(d[rindex][2])
                        if str(locus) not in whitelist:
                                rindex += 1
                                continue
                        if locus == lastlocus:allele += 1
                        else: allele = 1
                        lastlocus = locus
                        fout.write(str(d[rindex][0])+'\t'+str(d[rindex][2])+'_snp'+str(allele)+'\t0\t'+str(d[rindex][1])+'\t'+str(d[rindex][3])+'\t'+str(d[rindex][4])+'\n')
                        for column in row:
                                ind = d[8][cindex]
                                if cindex < 9:
                                        cindex += 1
                                        continue
                                if ind in snpdict.keys():
                                        a1 = str(re.sub(r'(\S)/\S:\d+:\S+', r'\1', d[rindex][cindex]))
                                        a2 = str(re.sub(r'\S/(\S):\d+:\S+', r'\1', d[rindex][cindex]))
                                        if a1 == '.': allele1 = 0
                                        else:
                                                if int(a1) == 1:
                                                        allele1 = str(d[rindex][4])
                                                elif int(a1) == 0:
                                                        allele1 = str(d[rindex][3])
                                        if a2 == '.': allele2 = 0
                                        else:
                                                if int(a2) == 1:
                                                        allele2 = str(d[rindex][4])
                                                elif int(a2) == 0:
                                                        allele2 = str(d[rindex][3])
                                        snpdict[ind][str(str(locus)+"_snp"+str(allele))] = [allele1, allele2]
                                else:
                                        a1 = str(re.sub(r'(\S)/\S:\d+:\S+', r'\1', d[rindex][cindex]))
                                        a2 = str(re.sub(r'\S/(\S):\d+:\S+', r'\1', d[rindex][cindex]))
                                        if a1 == '.': allele1 = 0
                                        else: allele1 = int(a1)
                                        if a2 == '.': allele2 = 0
                                        else: allele2 = int(a2)
                                        snpdict[ind] = {str(str(locus)+"_snp"+str(allele)) : [allele1, allele2]}
                                cindex += 1
                        rindex += 1
                fout.close()

                print "     Writing file..."
                fout=open(outfile,"w")

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        ind = str(popsdict[k][x])
                                        fout.write(str(ind)+'\t'+str(k)+'\t0\t0\t0\t0')
                                        for i in sorted(snpdict[ind].iterkeys()): #Look at one locus
                                                for a in sorted(snpdict[ind][i]): #Cycle through 1 or 2 alleles
                                                        fout.write('\t'+str(a))
                                        fout.write("\n")
                print "*** DONE! ***"
                fout.close()


        except IOError:
                print "     ** Error: Problems outputting file. Make sure the VCF file is included. Check the directory path. **"
                return 0

        return 1




# Output Phylip file
def Fasta2Phylip(num_pops, popsdict, seqsdict, gene_copies, num_sites, haplo, haplotypes, phylo_inform, breakpoints, locheaders):
        print "Outputting Phylip file..."

        try:

                if haplo == 1:
                        
                        print "     Cataloging unique sequences..."
                        globalseqsdict = {} #Structure: {LocusID : {AlleleID : DNAsequence} }
                        for ind in seqsdict.keys():
                                for locus in seqsdict[ind].keys():
                                        for allele_id, allele_seq in seqsdict[ind][locus].iteritems():
                                                if locus not in globalseqsdict.keys():
                                                        globalseqsdict[locus] = { 0 : seqsdict[ind][locus][allele_id] }
                                                        
                                                else:
                                                        if allele_seq not in globalseqsdict[locus].values():
                                                                new_allele = max(globalseqsdict[locus].keys()) + 1
                                                                globalseqsdict[locus][new_allele] = allele_seq
                        
                        print "     Building dictionary of haplotypes..."
                        haplodict = {} #Structure: {LocusID : {AlleleID : Haplotype} }
                        for locus in globalseqsdict.keys():
                                SNP_positions = []
                                for letter, letter_str in enumerate(globalseqsdict[locus][0]):
                                        default = globalseqsdict[locus][0][letter]
                                        for allele in globalseqsdict[locus].keys():
                                                if globalseqsdict[locus][allele][letter] != default and letter not in SNP_positions:
                                                        SNP_positions.append(letter)

                                for allele in globalseqsdict[locus].keys():
                                        hap = ""
                                        for i,j in enumerate(SNP_positions):
                                                hap += globalseqsdict[locus][allele][SNP_positions[i]]

                                        if locus not in haplodict.keys():
                                                haplodict[locus] = { 0 : hap }
                                        else:
                                                haplodict[locus][allele] = hap

                        print "     Cataloging haplotypes..."
                        newseqsdict = seqsdict #Structure: {SampleID : {LocusID : {AlleleID : Haplotype} } }
                        for ind in seqsdict.keys():
                                for locus in seqsdict[ind].keys():
                                        for allele_id, allele_seq in seqsdict[ind][locus].iteritems():
                                                newseqsdict[ind][locus][allele_id] = haplodict[locus][globalseqsdict[locus].keys()[globalseqsdict[locus].values().index(allele_seq)]]

                        print "     Counting haplotype SNP sites..."
                        new_num_sites = {} #Structure: {LocusID : NumberNucleotides}
                        for locus in haplodict.keys():
                                new_num_sites[locus] = len(haplodict[locus][0])

                if haplo == 2:
                        newseqsdict = seqsdict
                        new_num_sites = num_sites

                if haplotypes == 1:
                        nind = sum(len(v) for v in popsdict.itervalues()) * 2
                elif haplotypes == 2:
                        nind = sum(len(v) for v in popsdict.itervalues())
                elif haplotypes == 3:
                        nind = len(popsdict)

                OrderedLoci = []
                for key in sorted(new_num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        
        	if haplotypes == 3:
                        # Dictionary with of sequences by population
                        print "     Creating dictionary of population sequences..."
                        popseqsdict = {} #Structure: {PopID : {LocusID : [Sequences] } }
                        for pop in popsdict.iterkeys(): #Look at one population
                                for ind in newseqsdict.iterkeys(): #look at one individual
                                        if ind in popsdict[pop].iterkeys():
                                                for locus in newseqsdict[ind].iterkeys(): #Look at one locus
                                                        for allele, seq in newseqsdict[ind][locus].iteritems():
                                                                if pop not in popseqsdict.iterkeys():
                                                                        popseqsdict[pop] = {locus : [seq] }
                                                                else:
                                                                        if locus not in popseqsdict[pop].iterkeys():
                                                                                popseqsdict[pop][locus] = [seq]
                                                                        else:
                                                                                if seq not in popseqsdict[pop][locus]:
                                                                                        popseqsdict[pop][locus].append(seq)

                if phylo_inform == 1 or phylo_inform == 2:
                        # Remove SNPs (if haplo == 1) or loci (if haplo == 2) that aren't phylogenetically informative
                        fixed_sites = {} #Structure: {LocusID : [PositionIDs] }
                        print "     Removing loci that are not phylogenetically informative..."
                        if haplotypes == 1:
                                print "     Loci must have SNPs that are present for alternative alleles at 2+ haplotypes..."
                                for i in OrderedLoci: #Look at one locus
                                        for n in range(0,new_num_sites[i]): #Look at one bp position
                                                all_letter = []
                                                for k in sorted(newseqsdict.iterkeys()): #Look at one individual
                                                        if i in newseqsdict[k].keys(): #If individual has this locus
                                                                for a in newseqsdict[k][i].iterkeys(): #Look at one allele
                                                                        all_letter.append(newseqsdict[k][i][a][n])

                                                all_letter = collections.Counter(all_letter) #Count DNA character states across all haplotypes
                                                if phylo_inform == 1: all_letter = collections.Counter(y for y in all_letter.elements() if all_letter[y] >= 2) #Remove singletons (P uninformative)
                                                if len(all_letter) > 1: fixed = 1 #This position has fixed differences that are informative if fixed = 1
                                                elif len(all_letter) <= 1: fixed = 0
                                                
                                                if fixed == 1:
                                                        if i not in fixed_sites.keys():
                                                                fixed_sites[i] = [n]
                                                        else:
                                                                fixed_sites[i].append(n)

                                if haplo == 1:
                                        # Remove SNPs that aren't fixed and PI
                                        for k in sorted(newseqsdict.iterkeys()):
                                                for i in sorted(newseqsdict[k].iterkeys()):
                                                        if i in fixed_sites.keys():
                                                                for a in newseqsdict[k][i].iterkeys():
                                                                        seq = ""
                                                                        for n,m in enumerate(fixed_sites[i]):
                                                                                seq += newseqsdict[k][i][a][m]
                                                                        newseqsdict[k][i][a] = seq
                                                        else:
                                                                del newseqsdict[k][i]
                                        for i in new_num_sites.keys():
                                                if i in fixed_sites.keys():
                                                        new_num_sites[i] = len(fixed_sites[i])
                                                else:
                                                        del new_num_sites[i]
                                if haplo == 2:
                                        # Remove loci that don't have at least one fixed/PI SNP
                                        for k in sorted(newseqsdict.iterkeys()):
                                                for i in sorted(newseqsdict[k].iterkeys()):
                                                        if i not in fixed_sites.keys():
                                                                del newseqsdict[k][i]
                                        for i in new_num_sites.keys():
                                                if i not in fixed_sites.keys():
                                                        del new_num_sites[i]
                        if haplotypes == 2:
                                print "     Loci must have SNPs that are fixed for alternative alleles at 2+ individuals..."
                                for i in OrderedLoci: #Look at one locus
                                        for n in range(0,new_num_sites[i]): #Look at one bp position
                                                all_letter = []
                                                for k in sorted(newseqsdict.iterkeys()): #Look at one individual
                                                        if i in newseqsdict[k].keys(): #If individual has this locus
                                                                letters = []
                                                                for a in newseqsdict[k][i].iterkeys(): #Look at one allele
                                                                        if newseqsdict[k][i][a][n] not in letters:
                                                                                letters.append(newseqsdict[k][i][a][n])
                                                                letter = IUPAC(letters)
                                                                all_letter.append(letter)

                                                all_letter = collections.Counter(all_letter) #Count DNA character states across all individuals
                                                if phylo_inform == 1: all_letter = collections.Counter(y for y in all_letter.elements() if all_letter[y] >= 2) #Remove singletons (P uninformative)
                                                fixed = IUPAC_fixed(all_letter) #This position has fixed differences that are informative if fixed = 1

                                                if fixed == 1:
                                                        if i not in fixed_sites.keys():
                                                                fixed_sites[i] = [n]
                                                        else:
                                                                fixed_sites[i].append(n)

                                if haplo == 1:
                                        # Remove SNPs that aren't fixed and PI
                                        for k in sorted(newseqsdict.iterkeys()):
                                                for i in sorted(newseqsdict[k].iterkeys()):
                                                        if i in fixed_sites.keys():
                                                                for a in newseqsdict[k][i].iterkeys():
                                                                        seq = ""
                                                                        for n,m in enumerate(fixed_sites[i]):
                                                                                seq += newseqsdict[k][i][a][m]
                                                                        newseqsdict[k][i][a] = seq
                                                        else:
                                                                del newseqsdict[k][i]
                                        for i in new_num_sites.keys():
                                                if i in fixed_sites.keys():
                                                        new_num_sites[i] = len(fixed_sites[i])
                                                else:
                                                        del new_num_sites[i]
                                if haplo == 2:
                                        # Remove loci that don't have at least one fixed/PI SNP
                                        for k in sorted(newseqsdict.iterkeys()):
                                                for i in sorted(newseqsdict[k].iterkeys()):
                                                        if i not in fixed_sites.keys():
                                                                del newseqsdict[k][i]
                                        for i in new_num_sites.keys():
                                                if i not in fixed_sites.keys():
                                                        del new_num_sites[i]
                        if haplotypes == 3:
                                print "     Loci must have SNPs that are fixed for alternative alleles at 2+ populations..."
                                for i in OrderedLoci: #Look at one locus
                                        for n in range(0,new_num_sites[i]): #Look at one bp position
                                                all_letter = []
                                                for k in sorted(popseqsdict.iterkeys()): #Look at one population
                                                        if i in popseqsdict[k].keys(): #If population has this locus
                                                                letters = []
                                                                for a in range(0,len(popseqsdict[k][i])): #Look at one allele
                                                                        if popseqsdict[k][i][a][n] not in letters:
                                                                                letters.append(popseqsdict[k][i][a][n])
                                                                letter = IUPAC(letters)
                                                                all_letter.append(letter)

                                                all_letter = collections.Counter(all_letter) #Count DNA character states across all populations
                                                if phylo_inform == 1: all_letter = collections.Counter(y for y in all_letter.elements() if all_letter[y] >= 2) #Remove singletons (P uninformative)
                                                fixed = IUPAC_fixed(all_letter) #This position has fixed differences that are informative if fixed = 1
                                                if fixed == 1:
                                                        if i not in fixed_sites.keys():
                                                                fixed_sites[i] = [n]
                                                        else:
                                                                fixed_sites[i].append(n)
                                if haplo == 1:
                                        # Remove SNPs that aren't fixed and PI
                                        for k in sorted(popseqsdict.iterkeys()):
                                                for i in sorted(popseqsdict[k].iterkeys()):
                                                        if i in fixed_sites.keys():
                                                                for a in range(0,len(popseqsdict[k][i])):
                                                                        seq = ""
                                                                        for n,m in enumerate(fixed_sites[i]):
                                                                                seq += popseqsdict[k][i][a][m]
                                                                        popseqsdict[k][i][a] = seq
                                                        else:
                                                                del popseqsdict[k][i]
                                        for i in new_num_sites.keys():
                                                if i in fixed_sites.keys():
                                                        new_num_sites[i] = len(fixed_sites[i])
                                                else:
                                                        del new_num_sites[i]
                                if haplo == 2:
                                        # Remove loci that don't have at least one fixed/PI SNP
                                        for k in sorted(popseqsdict.iterkeys()):
                                                for i in sorted(popseqsdict[k].iterkeys()):
                                                        if i not in fixed_sites.keys():
                                                                del popseqsdict[k][i]
                                        for i in new_num_sites.keys():
                                                if i not in fixed_sites.keys():
                                                        del new_num_sites[i]


                        OrderedLoci = []
                        for key in sorted(new_num_sites.iterkeys()):
                                OrderedLoci.append(key)


                nbp = sum(new_num_sites.itervalues())
                        
                fout=open(outfile,"w")
                fout.write(str(nind) + '\t' + str(nbp) + '\n')
        	print "     Writing file..."
        	
        	if locheaders == 1:
        			fout.write('\t')
        			for i in OrderedLoci:
        					fout.write(i + '\t')
        			fout.write('\n')
                if haplotypes == 1:
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(newseqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                gene_copies = ['0','1']
                                                for gene_copy in gene_copies:
                                                        if gene_copy == '0':
                                                                ind = str(popsdict[k][x]) + 'a'
                                                                fout.write((ind+' ').ljust(10)); z=0
                                                        if gene_copy == '1':
                                                                ind = str(popsdict[k][x]) + 'b'
                                                                fout.write((ind+' ').ljust(10)); z=0
                                                        for i in OrderedLoci: #Look at one locus
                                                                if int(len(str(popsdict[k][x])))>9: print "     ** Error, Ind ID > 9 characters **"; exit(1)
                                                                if i not in newseqsdict[x].keys(): #If individual doesn't have this locus, write Ns
                                                                        z=0
                                                                        while z < int(new_num_sites[i]):
                                                                                fout.write("N"); z+=1
                                                                        if breakpoints == 1:
                                                                                fout.write("!")

                                                                else:
                                                                        for p in sorted(newseqsdict[x].iterkeys()): #Cycle through its loci
                                                                                if p == i:      #If this is the right locus
                                                                                        if int(len(newseqsdict[x][p]))==2:#If heterozygote
                                                                                                fout.write(str(newseqsdict[x][p][gene_copy]))
                                                                                                if breakpoints == 1:
                                                                                                        fout.write("!")
                                                                                        elif int(len(newseqsdict[x][p]))==1:#If homozygote
                                                                                                fout.write(str(newseqsdict[x][p]['0']))
                                                                                                if breakpoints == 1:
                                                                                                        fout.write("!")
                                                        fout.write('\n')
        	if haplotypes == 2:
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(newseqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                fout.write((ind+' ').ljust(10)); z=0
                                                for i in OrderedLoci: #Look at one locus
                                                        if int(len(str(popsdict[k][x])))>9: print "     ** Error, Ind ID > 9 characters **"; exit(1)
                                                        if i not in newseqsdict[x].keys(): #If individual doesn't have this locus, write Ns
                                                                z=0
                                                                while z < int(new_num_sites[i]):
                                                                        fout.write("N"); z+=1
                                                                if breakpoints == 1:
                                                                        fout.write("!")

                                                        else:
                                                                for p in sorted(newseqsdict[x].iterkeys()): #Cycle through its loci
                                                                        if p == i:      #If this is the right locus
                                                                                if int(len(newseqsdict[x][p]))==2:#If heterozygote
                                                                                        seq = ""
                                                                                        for n, m in enumerate(str(newseqsdict[x][p]['0'])):
                                                                                                seq += IUPAC([str(newseqsdict[x][p]['0'][n]), str(newseqsdict[x][p]['1'][n])])
                                                                                        fout.write(seq)
                                                                                        if breakpoints == 1:
                                                                                                fout.write("!")
                                                                                elif int(len(newseqsdict[x][p]))==1:#If homozygote
                                                                                        fout.write(str(newseqsdict[x][p]['0']))
                                                                                        if breakpoints == 1:
                                                                                                fout.write("!")
                                                fout.write('\n')
        	if haplotypes == 3:                        
                        for k in sorted(popseqsdict.iterkeys()): #Look at one population
                                if int(len(str(k)))>9: print "     ** Error, Pop ID > 9 characters **"; exit(1)
                                pop = str(k)
                                fout.write((pop+' ').ljust(10)); z=0
                                for i in OrderedLoci: #Look at one locus
                                        if i not in popseqsdict[k].keys(): #If individual doesn't have this locus, write Ns
                                                z=0
                                                while z < int(new_num_sites[i]):
                                                        fout.write("N"); z+=1
                                                if breakpoints == 1:
                                                        fout.write("!")

                                        else:
                                                for p in sorted(popseqsdict[k].iterkeys()): #Cycle through its loci
                                                        if p == i:      #If this is the right locus
                                                                if int(len(popseqsdict[k][p]))>1:#If polymorphic
                                                                        seq = ""
                                                                        for n, m in enumerate(str(popseqsdict[k][p][0])):
                                                                                letters = []
                                                                                for a in range(0,len(popseqsdict[k][p])):
                                                                                        if popseqsdict[k][p][a] not in letters:
                                                                                                letters.append(str(popseqsdict[k][p][a][n]))
                                                                                seq += IUPAC(letters)
                                                                        fout.write(seq)
                                                                        if breakpoints == 1:
                                                                                fout.write("!")
                                                                elif int(len(popseqsdict[k][p]))==1:#If monomorphic
                                                                        fout.write(str(popseqsdict[k][p][0]))
                                                                        if breakpoints == 1:
                                                                                fout.write("!")
                                fout.write('\n')
                print "*** DONE! ***"
                fout.close()


        except IOError:
                print "     ** Error: Problems outputting file. Check the directory path. **"
                return 0

        return 1




# Output G-Phocs file
def Fasta2GPhocs(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        print "Outputting G-Phocs file..."
        try:
                OrderedLoci = []
                fout=open(outfile,"w")
                
                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        
                numinds = sum(len(q) for q in popsdict.itervalues())
                
                fout.write(str(len(OrderedLoci)) + '\n\n')
				
                for i in OrderedLoci: #Look at one locus
                        fout.write(str(i) + '\t' + str(numinds) + '\t' + str(num_sites[i]) + '\n')
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write Ns
                                                        ind = str(popsdict[k][x])
                                                        fout.write((ind) + '\t'); z=0
                                                        while z < int(num_sites[i]):
                                                                fout.write("N"); z+=1


                                                else:
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                if p == i:      #If this is the right locus
                                                                        ind = str(popsdict[k][x])
                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                seq = ""
                                                                                for n, m in enumerate(str(seqsdict[x][p]['0'])):
                                                                                        seq += IUPAC([str(seqsdict[x][p]['0'][n]), str(seqsdict[x][p]['1'][n])])
                                                                                fout.write(ind + '\t' + seq + '\n')
                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                fout.write(ind + '\t' + str(seqsdict[x][p]['0'])+'\n')
                        fout.write('\n')
                print "*** DONE! ***"
                fout.close()

        except IOError:
                print "     ** Error: Problems outputting file. Check the directory path. **"
                return 0

        return 1




# Output Treemix file
def Fasta2Treemix(num_pops, popsdict, seqsdict, gene_copies, num_sites, one_snp):
        print "Outputting Treemix file..."

        try:
        	print "     Finding unique sequences..."
                # Build dictionary of unique sequences
        	globalseqsdict = {} #Structure: {LocusID : {AlleleID : DNAsequence} }
        	for ind in seqsdict.keys():
        		for locus in seqsdict[ind].keys():
        			for allele_id, allele_seq in seqsdict[ind][locus].iteritems():
                                        if locus not in globalseqsdict.keys():
                                                globalseqsdict[locus] = { 0 : seqsdict[ind][locus][allele_id] }
                                                
                                        else:
                                                if allele_seq not in globalseqsdict[locus].values():
                                                        new_allele = max(globalseqsdict[locus].keys()) + 1
                                                        globalseqsdict[locus][new_allele] = allele_seq

        	print "     Finding SNP positions..."
                # Find sequence positions that are SNPs
                SNPpositions = {} #Structure {LocusID : [Positions] }
        	for locus in globalseqsdict.keys():
                        SNP_positions = []
                        for letter, letter_str in enumerate(globalseqsdict[locus][0]):
                                default = globalseqsdict[locus][0][letter]
                                for allele in globalseqsdict[locus].keys():
                                        if globalseqsdict[locus][allele][letter] != default and letter not in SNP_positions:
                                                SNP_positions.append(letter)
                        SNPpositions[locus] = SNP_positions

        	print "     Cataloging SNPs..."
                # Build dictionary of SNPs
                SNPdict = {} #Structure: {SampleID : {LocusID : {PositionID : {AlleleID : SNP} } } }
                for ind in seqsdict.keys():                               
                        for locus in seqsdict[ind].keys():
                                for position in SNPpositions[locus]:
                                        for allele_id in seqsdict[ind][locus]:
                                                if ind not in SNPdict.iterkeys():
                                                        SNPdict[ind] = {locus : {position: {allele_id : seqsdict[ind][locus][allele_id][position]} } }
                                                else:
                                                        if locus not in SNPdict[ind].iterkeys():
                                                                SNPdict[ind][locus] = {position: {allele_id : seqsdict[ind][locus][allele_id][position]} }
                                                        else:
                                                                if position not in SNPdict[ind][locus].iterkeys():
                                                                        SNPdict[ind][locus][position] = {allele_id : seqsdict[ind][locus][allele_id][position]}
                                                                else:
                                                                        if allele_id not in SNPdict[ind][locus][position].iterkeys():
                                                                                SNPdict[ind][locus][position][allele_id] = seqsdict[ind][locus][allele_id][position]

        	print "     Doing more SNP cataloging..."
                # Dictionary with both SNP letters present at each SNP locus
                SNPletters = {} #Structure: {LocusID : {PositionID : [Letters] } }
                
        	print "     Flagging SNPs that are not biallelic..."
                # Flag SNPs with >2 alleles
                globalSNPalleleCounter = {} #Structure: {LocusID : {PositionID : NumAlleles } }
                for ind in SNPdict.keys():                               
                        for locus in SNPdict[ind].keys():
                                for position in SNPdict[ind][locus].keys():
                                        for allele_id, allele_snp in SNPdict[ind][locus][position].iteritems():
                                                if locus not in globalSNPalleleCounter.iterkeys():
                                                        globalSNPalleleCounter[locus] = {position : 1 }
                                                        SNPletters[locus] = {position : [allele_snp]}
                                                else:
                                                        if position not in globalSNPalleleCounter[locus].iterkeys():
                                                                globalSNPalleleCounter[locus][position] = 1
                                                                SNPletters[locus][position] = [allele_snp]
                                                        else:
                                                                if allele_snp not in SNPletters[locus][position]:
                                                                        SNPletters[locus][position].append(allele_snp)
                                                                        globalSNPalleleCounter[locus][position] += 1


        	print "     Counting population SNP alleles..."
                # Dictionary with count for each SNP allele
                SNPpopCount = {} #Structure: {PopID : {LocusID : {PositionID : {SNP : Count} } } }
                for pop in popsdict.iterkeys(): #Look at one population
                        for ind in SNPdict.iterkeys(): #look at one individual
                                if ind in popsdict[pop].iterkeys():
                                        for locus in SNPdict[ind].iterkeys(): #Look at one locus
                                                for position in SNPdict[ind][locus].iterkeys(): #Look at one SNP position in locus
                                                        for allele_id, allele_snp in SNPdict[ind][locus][position].iteritems(): #Look at one SNP allele
                                                                if pop not in SNPpopCount.iterkeys():
                                                                        SNPpopCount[pop] = {locus : {position: {allele_snp : 1} } }
                                                                else:
                                                                        if locus not in SNPpopCount[pop].iterkeys():
                                                                                SNPpopCount[pop][locus] = {position: {allele_snp : 1} }
                                                                        else:
                                                                                if position not in SNPpopCount[pop][locus].iterkeys():
                                                                                        SNPpopCount[pop][locus][position] = {allele_snp : 1}
                                                                                else:
                                                                                        if allele_snp not in SNPpopCount[pop][locus][position].iterkeys():
                                                                                                SNPpopCount[pop][locus][position][allele_snp] = 1
                                                                                        else:
                                                                                                SNPpopCount[pop][locus][position][allele_snp] += 1

                OrderedLoci = []
                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)

                fout=open(outfile,"w")
                pops = []
                for pop in sorted(SNPpopCount.iterkeys()): #Look at one population
                        fout.write(str(pop)+' ') #Print each population as column header
                        pops.append(pop)
                fout.write('\n')

                numsnps = 0
                for locus in globalSNPalleleCounter:
                        for position in globalSNPalleleCounter[locus]:
                                if globalSNPalleleCounter[locus][position] <= 2:
                                        numsnps += 1

                # Remove monomorphic loci from globalseqsdict, SNPpositions, OrderedLoci
                keeploci = globalSNPalleleCounter.keys()
                rem = np.setdiff1d(SNPpositions.keys(), keeploci)
                for i in rem: del SNPpositions[i]
                rem = np.setdiff1d(globalseqsdict.keys(), keeploci)
                for i in rem: del globalseqsdict[i]
                rem = np.setdiff1d(OrderedLoci, keeploci)
                for i in rem: OrderedLoci.remove(i)
                
                #Print out counts of allele by population
                print "     Writing file..."
                locus = 0
                while locus < len(OrderedLoci):
                        locus_index = OrderedLoci[locus]
                        position = 0
                        if one_snp == 1:
                                num_SNPs_per_locus = 1
                        elif one_snp == 2: 
                                num_SNPs_per_locus = len(SNPpositions[locus_index])
                        while position < num_SNPs_per_locus:
                                position_index = SNPpositions[locus_index][position]
                                if globalSNPalleleCounter[locus_index][position_index] <= 2:
                                        for pop in pops:
                                                if locus_index in SNPpopCount[pop]:
                                                        for SNP in sorted(SNPpopCount[pop][locus_index][position_index].iterkeys()):
                                                                num_alleles = len(SNPpopCount[pop][locus_index][position_index])
                                                                this_SNP_first = 2
                                                                other = 2
                                                                if SNPletters[locus_index][position_index][0] == SNP and SNPletters[locus_index][position_index][0] < SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 1
                                                                        other = 1
                                                                if SNPletters[locus_index][position_index][1] == SNP and SNPletters[locus_index][position_index][0] < SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 0
                                                                        other = 0
                                                                if SNPletters[locus_index][position_index][0] == SNP and SNPletters[locus_index][position_index][0] > SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 0
                                                                        other = 1
                                                                if SNPletters[locus_index][position_index][1] == SNP and SNPletters[locus_index][position_index][0] > SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 1
                                                                        other = 0
                                                                if num_alleles == 1:#Only one SNP allele present in population                                                                       
                                                                        if this_SNP_first == 1:
                                                                                fout.write(str(SNPpopCount[pop][locus_index][position_index][SNP]) + ",0")
                                                                        if this_SNP_first == 0:
                                                                                fout.write("0," + str(SNPpopCount[pop][locus_index][position_index][SNP]))
                                                                        fout.write(' ')
                                                                if num_alleles == 2:#Both SNP alleles in population
                                                                        if this_SNP_first == 1:
                                                                                fout.write(str(SNPpopCount[pop][locus_index][position_index][SNP]) + "," + str(SNPpopCount[pop][locus_index][position_index][SNPletters[locus_index][position_index][other]]))
                                                                        if this_SNP_first == 0:
                                                                                fout.write(str(SNPpopCount[pop][locus_index][position_index][SNPletters[locus_index][position_index][other]]) + "," + str(SNPpopCount[pop][locus_index][position_index][SNP]))
                                                                        fout.write(' ')
                                                                        break
                                                else:
                                                        fout.write("0,0 ")
                                        fout.write("\n")
                                position += 1
                        locus += 1                                      

                print "*** DONE! ***"
                fout.close()

        except IOError:
                print "     ** Error: Problems outputting file. Check the directory path. **"
                return 0

        return 1




# Output file of haplotypes
def Fasta2Haplotype(num_pops, popsdict, seqsdict, gene_copies, num_sites, HaploChoice, title, FourOrSix):
        try:
                OrderedLoci = []
                LociOrdered = {}
                OrderedPops = []
                PopsOrdered = {}

                fout=open(outfile,"w")

                count = 1
                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        LociOrdered[key] = count
                        count += 1

                count = 1
                for key in sorted(popsdict.iterkeys()):
                        OrderedPops.append(key)
                        PopsOrdered[key] = count
                        count += 1

        	print "Cataloging unique sequences..."
                #Create dictionary of loci containing dictionaries of unique sequences/integers
                newseqdict = {} #Build structure: {LocusID : {Sequence : UniqueInteger} }
                for x in sorted(seqsdict.iterkeys()):
                	for p in sorted(seqsdict[x].iterkeys()):
                		for a in sorted(seqsdict[x][p].iterkeys()):
                			if p in newseqdict.keys():
                				if str(seqsdict[x][p][a]) not in newseqdict[p].keys() and "?" not in seqsdict[x][p][a]:
                					locusints = []
                					for i in sorted(newseqdict[p].itervalues()):
                						locusints.append(i)
                						newseqdict[p][str(seqsdict[x][p][a])] = max(locusints) + 1
                                                elif "?" in seqsdict[x][p][a]:
                                                        newseqdict[p][str(seqsdict[x][p][a])] = 0

                			else:
                                                if "?" not in seqsdict[x][p][a]:
                                                        newseqdict[p] = {str(seqsdict[x][p][a]):1}
                                                else:
                                                        newseqdict[p] = {str(seqsdict[x][p][a]):0}

        	print "Creating dictionary of unique integers for each haplotype..."
                #Convert sequences into unique integers by locus
                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through all loci
                		for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                			seqsdict[x][p][a] = newseqdict[str(p)][str(seqsdict[x][p][a])]

                if HaploChoice == 1: # Structure format
                        print "Outputting Structure file..."
                        fout.write('\t')
                        for i in OrderedLoci:
                                fout.write('\t'+str(i))
                        fout.write('\n')
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                count = 0
                                                while count < 2:
                                                        ind = str(popsdict[k][x])
                                                        fout.write(ind+'\t'+k+'\t')
                                                        for i in OrderedLoci: #Look at one locus
                                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write 0s
                                                                        fout.write('0\t')

                                                                else:
                                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                        count += 1
                                                        fout.write("\n")
                        print "*** DONE! ***"
                        fout.close()

                if HaploChoice == 2: #Genepop format
                        print "Outputting Genepop file..."
                        fout.write(title+'\n')
                        for i in OrderedLoci:
                                fout.write(str(i)+'\n')

                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                fout.write('Pop\n')
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                fout.write(ind+' ,  ')
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write 0s
                                                                if FourOrSix == 1:
                                                                        fout.write('0000\t')
                                                                if FourOrSix == 2:
                                                                        fout.write('000000\t')

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        count = 0
                                                                        while count < 2:
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0):
                                                                                                                if FourOrSix == 1:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(2))
                                                                                                                if FourOrSix == 2:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(3))
                                                                                                        if (count==1 and int(a)==1):
                                                                                                                if FourOrSix == 1:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(2)+'\t')
                                                                                                                if FourOrSix == 2:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(3)+'\t')
                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        if FourOrSix == 1:
                                                                                                                fout.write(str(seqsdict[x][p][a]).zfill(2)+str(seqsdict[x][p][a]).zfill(2)+'\t')
                                                                                                        if FourOrSix == 2:
                                                                                                                fout.write(str(seqsdict[x][p][a]).zfill(3)+str(seqsdict[x][p][a]).zfill(3)+'\t')
                                                                                                        count += 1
                                                                                count += 1
                                                fout.write("\n")
                        print "*** DONE! ***"
                        fout.close()

                        fout=open(outfile_pops,"w")
                        for i in sorted(PopsOrdered.iterkeys()):
                                fout.write(str(PopsOrdered[i])+'\t'+str(i)+'\n')
                        fout.close()

                if HaploChoice == 3: #Allele frequency by locus X population
                        print "Outputting allele frequency X population matrix..."
                        fout.write('\t')
                        for p in sorted(newseqdict.iterkeys()):
                                for s in sorted(newseqdict[p].itervalues()):
                                        fout.write(str(p)+'_'+str(s)+'\t') #Print each locus/allele combo as column header
                        fout.write('\n')

                        popfreq = {} #Build Structure: {PopulationID : {LocusID : {AlleleInteger : Count} } }
                                     #First build empty dictionary structure with 0 count for each allele

                        print "     Creating dictionary of population allele frequencies..."
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                if k in sorted(popfreq.keys()):
                                        for p in sorted(newseqdict.iterkeys()): #Look at one locus
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()): #Look at one allele integer
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                else:
                                        popfreq[k] = {}
                                        for p in sorted(newseqdict.iterkeys()):
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass

                        print "     Tabulating population allele frequencies..."
                        #Add counts of each allele
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                for i in OrderedLoci: #Look at one locus
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                count = 0
                                                                while count < 2:
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                        popfreq[k][p][seqsdict[x][p][a]] += 1
                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                popfreq[k][p][seqsdict[x][p][a]] += 2
                                                                                                count += 1
                                                                        count += 1

                        print "     Writing file..."
                        #Print out frequencies of allele by population
                        for k in sorted(popfreq.iterkeys()):
                                fout.write(str(k)+'\t')
                                for p in sorted(popfreq[k].iterkeys()):
                                        total = 0
                                        for n in sorted(popfreq[k][p].itervalues()):
                                                total += n
                                        for n in sorted(popfreq[k][p].iterkeys()):
                                                if total > 0:
                                                        x = str(Decimal(str(popfreq[k][p][n])).quantize(Decimal('0.00001'))/Decimal(str(total)).quantize(Decimal('0.00001')))
                                                        fout.write(str(Decimal(str(x)).quantize(Decimal('0.00001')))+'\t')
                                                else:
                                                        fout.write(str(Decimal(str(popfreq[k][p][n])).quantize(Decimal('0.00001')))+'\t')
                                fout.write('\n')
                        print "*** DONE! ***"
                        fout.close()

                if HaploChoice == 4: #Sambada format
                        print "Outputting SamBada file..."
                        fout.write('\t')
                        for p in sorted(newseqdict.iterkeys()):
                                for s in sorted(newseqdict[p].itervalues()):
                                        fout.write(str(p)+'_'+str(s)+'\t') #Print each locus/allele combo as column header
                        fout.write('\n')

                        allelecount = {} #Build Structure: {IndividualID : {LocusID : {AlleleInteger : Count} } }
                                     #First build empty dictionary structure with 0 count for each allele

                        print "     Creating dictionary of allele counts..."
                        for i in sorted(popsdict.iterkeys()): #Look at one population
                                for k in sorted(popsdict[i].itervalues()): #Look at one individual
                                        if k in sorted(allelecount.keys()):
                                                for p in sorted(newseqdict.iterkeys()): #Look at one locus
                                                        if p in sorted(allelecount[k].iterkeys()):
                                                                for n in sorted(newseqdict[p].itervalues()): #Look at one allele integer
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass
                                                        else:
                                                                allelecount[k][p] = {}
                                                                for n in sorted(newseqdict[p].itervalues()):
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass
                                        else:
                                                allelecount[k] = {}
                                                for p in sorted(newseqdict.iterkeys()):
                                                        if p in sorted(allelecount[k].iterkeys()):
                                                                for n in sorted(newseqdict[p].itervalues()):
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass
                                                        else:
                                                                allelecount[k][p] = {}
                                                                for n in sorted(newseqdict[p].itervalues()):
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass

                        print "     Tabulating allele counts..."
                        #Add counts of each allele
                        for k in sorted(popsdict.iterkeys()): #Look at one individual
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, input -1s
                                                                for a in sorted(allelecount[popsdict[k][x]][i].iterkeys()): #Cycle through 1 or 2 alleles
                                                                        allelecount[popsdict[k][x]][i][a] = -1

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        count = 0
                                                                        while count < 2:
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                allelecount[popsdict[k][x]][p][seqsdict[x][p][a]] += 1
                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        allelecount[popsdict[k][x]][p][seqsdict[x][p][a]] += 2
                                                                                                        count += 1
                                                                                count += 1

                        print "     Writing file..."
                        #Print out counts of allele by individual
                        for k in sorted(allelecount.iterkeys()):
                                count = 0
                                while count < 2:
                                        if count == 0: fout.write(str(k)+'a\t')
                                        if count == 1: fout.write('\n'+str(k)+'b\t')
                                        for p in sorted(allelecount[k].iterkeys()):
                                                for n in sorted(allelecount[k][p].iterkeys()):
                                                        if count == 0:
                                                                if allelecount[k][p][n] == 0:
                                                                        fout.write('0\t')
                                                                if allelecount[k][p][n] == 1:
                                                                        allelecount[k][p][n] = 0
                                                                        fout.write('1\t')
                                                                if allelecount[k][p][n] == 2:
                                                                        allelecount[k][p][n] = 1
                                                                        fout.write('1\t')
                                                                if allelecount[k][p][n] == -1:
                                                                        fout.write('NaN\t')
                                                                if allelecount[k][p][n] > 2:
                                                                        print "     ** Error! Program retained more than 2 of a particular allele for one individual/locus. Check code. **"
                                                                        exit(1)
                                                                        # This shouldn't ever happen but the code hasn't been thoroughly tested.
                                                        elif count == 1:
                                                                if allelecount[k][p][n] == 0:
                                                                        fout.write('0\t')
                                                                if allelecount[k][p][n] == 1:
                                                                        allelecount[k][p][n] = 0
                                                                        fout.write('1\t')
                                                                if allelecount[k][p][n] == -1:
                                                                        fout.write('NaN\t')
                                                                if allelecount[k][p][n] > 1:
                                                                        print "     ** Error! Program retained more than 2 of a particular allele for one individual/locus. Check code. **"
                                                                        exit(1)
                                                                        # This shouldn't ever happen but the code hasn't been thoroughly tested.
                                        count += 1

                                fout.write('\n')
                        print "*** DONE! ***"
                        fout.close()

                if HaploChoice == 5: #Bayescan format
                        print "Outputting Bayescan file..."
                        fout.write('[loci]='+str(len(OrderedLoci))+'\n\n[populations]='+str(num_pops)+'\n')


                        popfreq = {} #Build Structure: {PopulationID : {LocusID : {AlleleInteger : Count} } }
                                     #First build empty dictionary structure with 0 count for each allele

                        print "     Creating dictionary of allele counts..."
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                if k in sorted(popfreq.keys()):
                                        for p in sorted(newseqdict.iterkeys()): #Look at one locus
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()): #Look at one allele integer
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                else:
                                        popfreq[k] = {}
                                        for p in sorted(newseqdict.iterkeys()):
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass

                        print "     Tabulating allele counts..."
                        #Add counts of each allele
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                for i in OrderedLoci: #Look at one locus
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                count = 0
                                                                while count < 2:
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                        popfreq[k][p][seqsdict[x][p][a]] += 1
                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                popfreq[k][p][seqsdict[x][p][a]] += 2
                                                                                                count += 1
                                                                        count += 1

                        print "     Writing file..."
                        #Print out Bayescan format (Locus  \t  NumGeneCopiesPop  \t  NumAllelesAtLocus  \t  Counts  \t  For  \t  Each  \t  Allele...)
                        for k in sorted(PopsOrdered.iterkeys()):
                                fout.write('\n[pop]='+str(PopsOrdered[k])+'\n')
                                for p in sorted(popfreq[k].iterkeys()):
                                        genecopies = 0
                                        for n in sorted(popfreq[k][p].iterkeys()):
                                                genecopies += popfreq[k][p][n]
                                        fout.write(str(LociOrdered[p])+'\t'+str(genecopies)+'\t'+str(len(newseqdict[p])))
                                        for n in sorted(popfreq[k][p].iterkeys()):
                                                fout.write('\t'+str(popfreq[k][p][n]))
                                        fout.write('\n')

                        fout.close()

                        fout=open(outfile_loci,"w")
                        for i in sorted(LociOrdered.iterkeys()):
                                fout.write(str(LociOrdered[i])+'\t'+str(i)+'\n')
                        fout.close()

                        fout=open(outfile_pops,"w")
                        for i in sorted(PopsOrdered.iterkeys()):
                                fout.write(str(PopsOrdered[i])+'\t'+str(i)+'\n')
                        print "*** DONE! ***"
                        fout.close()

                if HaploChoice == 6: #Arlequin format
                        OrderedLoci = []
                        fout=open(outfile,"w")
                        print "Outputting Arlequin file..."

                        fout.write("[Profile]\n\n\t\"" + title + "\"\n\n\t\tNbSamples=" + str(num_pops))
                        fout.write("\n\t\tGenotypicData=1\n\t\tGameticPhase=0\n\t\tDataType=STANDARD\n\t\t")
                        fout.write("LocusSeparator=TAB\n\t\tMissingData=\"?\"\n\n\n[Data]\n\n\t[[Samples]]\n\n")

                        for key in sorted(num_sites.iterkeys()):
                                OrderedLoci.append(key)

                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                fout.write("\t\tSampleName= \"Pop_" + str(k) + "\"\n\t\tSampleSize=" + str(int(gene_copies[k])/2) + "\n\t\tSampleData={\n")
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                count = 0
                                                while count < 2:
                                                        ind = str(popsdict[k][x])
                                                        if count == 0: fout.write(ind+'\t1\t')
                                                        if count == 1: fout.write('\t\t')
                                                        for i in OrderedLoci: #Look at one locus
                                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write ?s
                                                                        fout.write("?\t")

                                                                else:
                                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                        count += 1
                                                        fout.write('\n')
                                fout.write("}\n")
                        print "*** DONE! ***"
                        fout.close()

                if HaploChoice == 7: # GenAlEx format
                        print "Outputting GenAlEx file..."

                        fout=open(outfile,"w")

                        num_inds = 0
                        for k in popsdict.iterkeys():
                                num_inds += int(len(popsdict[k]))

                        fout.write(str(len(OrderedLoci))+'\t'+str(num_inds)+'\t'+str(num_pops))
                        for k in sorted(popsdict.iterkeys()):
                                fout.write('\t'+str(len(popsdict[k])))

                        fout.write('\n\t\t')
                        for k in sorted(popsdict.iterkeys()):
                                fout.write('\t'+str(k))
                        fout.write('\nIndID\tPopID\t')

                        count=0
                        for i in OrderedLoci:
                                if count > 0:
                                        fout.write('\t\t')
                                fout.write(str(i))
                                count += 1
                        fout.write('\n')

                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                fout.write(str(ind)+'\t'+str(k)+'\t')
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write 0s
                                                                fout.write('0\t0\t')

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                fout.write(str(seqsdict[x][p][a])+'\t'+str(seqsdict[x][p][a])+'\t')
                                                fout.write("\n")
                        print "*** DONE! ***"
                        fout.close()

        except IOError:
                print "     ** Error: Problems outputting file. Check the directory path. ** "
                return 0

        return 1




"""
Execute program options here.
"""

if UseCoverage == 1 and CoverageCutoff != 0:
        covdict = LocusCoverage(sys.argv[4])
else:
        covdict = {}

seqsdict, num_sites = Seqs(outtype, clipcutsite, cutsite1, cutsite2, CoverageCutoff, covdict)

popsdict, num_pops, gene_copies = Pops(sys.argv[3], seqsdict)

if monomorphic_filter == 1 or hetero_filter == 1 or monomorphic_filter2 == 1:
        seqsdict, num_sites = LocusRemoval(seqsdict, popsdict, gene_copies, num_sites, monomorphic_filter, hetero_filter, heterocutoff)

if allele_filter == 1:
        seqsdict, num_sites = AlleleRemoval(seqsdict, popsdict, gene_copies, num_sites, allele_threshold, allele_pop_threshold, allele_filter)

if missing_data_filter == 1:
        seqsdict, num_sites = MissingData(seqsdict, popsdict, gene_copies, num_sites, locus_threshold, locus_pop_threshold, ind_threshold)

if choice == 1:
	Job2 = Fasta2Migrate(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 2:
	Job2 = Fasta2Arlequin(num_pops, popsdict, seqsdict, gene_copies, num_sites, title)
if choice == 3:
	Job2 = Fasta2DIYABC(num_pops, popsdict, seqsdict, gene_copies, num_sites, title)
if choice == 4:
	Job2 = Fasta2LFMM(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 5:
	Job2 = Fasta2Phylip(num_pops, popsdict, seqsdict, gene_copies, num_sites, haplo, haplotypes, phylo_inform, breakpoints, locheaders)
if choice == 6:
	Job2 = Fasta2GPhocs(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 7:
	Job2 = Fasta2Treemix(num_pops, popsdict, seqsdict, gene_copies, num_sites, one_snp)
if choice == 8:
	Job2 = Fasta2Haplotype(num_pops, popsdict, seqsdict, gene_copies, num_sites, HaploChoice, title, FourOrSix)
