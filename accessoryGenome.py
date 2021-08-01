#! /usr/bin/env python
import sys
import re
import argparse
import glob
import os
import copy
from collections import defaultdict

"""
Author: Conor Meehan 14/07/2021

Script creates a map of presence/absence of accessory genome from lists of deleted genes and BLAST result files

Input:
--genomes A file that lists the genome names, one per line. NOTE: should be identical to how they are in the filenames of the deleted genes
--deleted A folder of files, each with a filename in the format genome1_genome2_* (where * can be anything). The genome1/2 names match exactly those in the --genomes file and the contents a list of genes present in genome1 and absent from genome2
--blast a BLAST file that contains all pairwise BLASTS between genes in these genomes. This is likely a combined file from individual BLASTs. It is used to check for orthologs

Output:
A tab delimited table where genes are rows and genomes are columns and there is a 1 if that gene is present in that genome and a 0 if absent
Genes are combined if there are orthologs of each other based on the BLAST file. The name then will be genes1_gene2 up to however many genes are combined

To run:
python accessoryGenome.py --genomes GenomeFile --deleted DeletedGenesFolder --blast blastFile

"""	

#parse the inputs
parser = argparse.ArgumentParser(description= 'Script creates a map of presence/absence of accessory genome from lists of deleted genes and BLAST result files')
parser.add_argument('--genomes', required=True, help='A file that lists the genome names, one per line')
parser.add_argument('--deleted', required=True, help='A folder containing files where deleted genes are listed')
parser.add_argument('--blast', required=True, help='BLAST data table (outfmt 6)')

args = parser.parse_args()

#read in the files
try:
	genomesF=open(args.genomes, 'r')
except IOError:
	print("\n Genome names file not found.")
	sys.exit()
try:
	blastF=open(args.blast, 'r')
except IOError:
	print("\n BLAST file not found.")
	sys.exit()

#put the genome namnes into a list
genomes=[]
while 1:
	line=genomesF.readline()
	if not line:
		break
	line=line.rstrip()
	genomes.append(line)
genomesF.close()


#Go through the deleted genes file and create a growing dictionary of gene names with a dictionary of genome name and presence/absence. Assume presence unless listed in a file as deleted
genes={}

#Create a dictionary of genome names and 1's to put as default for a new gene
pres_abs_default={}
for genome in genomes:
	pres_abs_default[genome]=1

#open the files one by one in the folder
print("Creating gene presence/absence map")
files=glob.glob(args.deleted+'/*')
for file in files:
	#get the genome names involved by first cutting off the folder name, then extracting the first genome and then the 2nd
	filename=re.sub(".*\/","",file)
	for genome in genomes:
		if filename.startswith(genome+"_"):
			g1=genome
		if "_"+genome+"_" in filename:
			g2=genome
	
	#open the file and read in all the genes. If gene not seen before, put a 1 for all genomes for that gene and then replace as a 0 for the g2. If not new, change the g2 to a 0
	try:
		F=open(file, 'r')
	except IOError:
		print("\n Deleted genes file "+file+" not found.")
		sys.exit()
	while 1:
		gene=F.readline()
		if not gene:
			break
		gene=gene.rstrip()
		if gene in genes.keys():
			genes[gene][g2]=0
		else:
			genes[gene]=copy.deepcopy(pres_abs_default)
			genes[gene][g2]=0
	F.close()

#put the BLAST result into a dictionary where the key is the gene name and the values is a list of all orthologs
blast=defaultdict(list)
while 1:
	line=blastF.readline()
	if not line:
		break
	#keep only the first 2 sections, which are the names
	sections=line.split("\t")
	if sections[1] not in blast[sections[0]]:
		blast[sections[0]].append(sections[1])
	if sections[0] not in blast[sections[1]]:
		blast[sections[1]].append(sections[0])	
blastF.close()

#go through the genes and if they are orthologs (i.e. blast dictionary lists them together, thus have a significant hit), combine them
print("Combining orthologs (this could take a while). Number of genes to check: "+str(len(genes.keys())))
orthologs=[]
genesInGroup=[]
for gene1 in genes:
	for gene2 in genes:
		#ensure no self hits
		if gene1==gene2:
			continue
		
		#check the blast dictionary for the gene pairing in either direction	
		if gene2 in blast[gene1] or gene1 in blast[gene2]:
				orthologs.append([gene1,gene2])	

#Go through the ortholog pairs and combine to create larger homolog groups where each gene is in a homolog set only once
#code taken from https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements and then modified
homologs = []

while len(orthologs)>0:
    first, *rest = orthologs
    first = set(first)

    lf = -1
    while len(first)>lf:
        lf = len(first)

        rest2 = []
        for r in rest:
            if len(first.intersection(set(r)))>0:
                first |= set(r)
            else:
                rest2.append(r)     
        rest = rest2

    homologs.append(first)
    orthologs = rest
  
#go through the list of genes and if not in a homolog group, add them as their own group:
for gene in genes:
	inGroup=0
	for group in homologs:
		if gene in group:
			inGroup=1
	if inGroup==0:
		homologs.append([gene])	 

print("Combining done, there are "+str(len(homologs))+" homolog groups")

#combine the gene presence/absence patterns for the orthologs into one entry
accessGenome={}
accNum=1
altNum=1
for gene in genes:
	inGroup=0
	for group in homologs: #check all the homolog groups for this gene
		if gene in group: #if found, then check if an entry has already been made in accessGenome
			geneFound=0
			for acc in accessGenome:
				if gene in accessGenome[acc][0]:
					geneFound=1
					if genes[gene]!=accessGenome[acc][1]:
						print("Inconsistency found for "+str(accessGenome[acc][0])+". Added as separate entries; please check your BLAST outputs for non-recipricol hits")
						accessGenome[acc+"_alt"+str(altNum)]=[accessGenome[acc][0],genes[gene]]
						altNum+=1
					break	
			
			if geneFound==0: #homolog group hasnt been added yet so add it and use this genes distribution pattern for the presence/absence (should be the same for all but that is checked above)
				accSorted=list(group)
				accSorted.sort()
				accessGenome[str(accNum)]=[accSorted,genes[gene]]
				accNum+=1

#create an output table of all the data, Put the homolog number in column 1, the lists of all the genes in column 2 and then the 1/0 in te other columns

#create a save file
try:
	save=open("AccessoryGenomeTable.txt", 'w')
except IOError:
	print("\n No room for save file.")
	sys.exit()

#create the header and save to file
header="Homolog group\tGenes\t"+"\t".join(genomes)+"\n"
save.write(header)

for group in accessGenome:
	save.write(group+"\t")
	save.write(", ".join(accessGenome[group][0]))
	for genome in genomes:
		save.write("\t"+str(accessGenome[group][1][genome]))
	save.write("\n")
save.close()
sys.exit()	