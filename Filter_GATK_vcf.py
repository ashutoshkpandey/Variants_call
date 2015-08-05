##This python script takes the vcf output from GATK Unified Genotyper and filter low quality variants based multiple criteria as mentioned below.
## The script has been inspired from a similar script from Heng Li but here I have added many more criteria to filter or flag low confident reads.
## I have implemented some global filters including removal of adjacent indels, snps adjacent to indels and multiple snps or indels within a very short stretch.
## Please reference: A genomic analysis of allele-specific expression in the mouse liver submitted in G3/2015/018309

#!/usr/bin/python

# Local Filters

# FS Strand Bias
# BB Base Qual rank sum test
# EB End Distance Bias
# MB MapQ Bias
# AF1 EM estimate of the site allele frequency of the strongest non-reference allele
# ABS Minimum non-reference supporting reads for the SNP call
# ABI Minimum non-reference supporting reads for the Indel call
# d low depth
# D high depth
# N Reference is N
# QS min RMS mapping quality (SNP only)
# QI min RMS mapping quality (Indel only)
# Qual Variant calling quality 

# Global Filters 

# W Adjacent Indel distance (Indel only) (if two gaps or indels lie within W (default 25 bp), indel with maximum supporting reads will be marked as PASS. In case of a tie, Indel with higher ratio of # of non-reference allele supporting reads to total depth will be preferred.
# w Adjacent SNP distance between Indel and SNP to remove SNPs.
# c Cluster of SNP
# sw SNP cluster window

import sys,os,random
import getopt

def usage():
	print '''usage: varfilter.py [options] input_vcf 

Options: 	 -FS1 Float Phred based p-value for strand bias (Higher phred = bias),SNP
		 -FS2 Float Phred based p-value for strand bias, Indel
		 -MB Float minimum p-value for MapQ bias (Only for SNP)
		 -EB1 Float GATK based value for end distance  bias(SNP)
		 -EB2 Float GATK based value for end distance  bias(Indel)
		 -AF1 Float minimum site allele frequency of the strongest non-reference allele
		 -ABS INT minimum alternate bases needed to call SNPs
		 -ABI INT minimum alternate bases needed to call Indels
                 -d INT minimum read depth
                 -D INT maximum read depth
                 -QS INT minimum RMS mapping quality for SNPs 
		 -QI INT minimum RMS mapping quality for Indels
		 -Qual  INT minimum variant calling quality
		 		
		 -w INT window size for filtering SNPs adjacent to Indel
		 -W INT window size for filtering adjacent Indels
		 -c INT minimum number of SNPs in 'sw' window to call cluster
		 -sw INT SNP window to determine SNP cluster'''

strandbias1 = 75
strandbias2 = 200
mapbias = -12.5
enddistbias1 = -10.0
enddistbias2 = -20.0
nonreffreq = .95
minaltsnp = 3
minaltindel = 3
mindepth = 3
maxdepth = 300
minsnpmapq = 30
minindelmapq = 20
variantqual = 50

gapgapwin = 20
snpgapwin = 10
snpcluster = 3
snpclusterwindow = 10
chromosome = ""

args = []
args = sys.argv[1:]
print args

input = open(args[0])
chromosome = str(args[1])

#rand = random.randrange(0,1000000)

tempfile = open("/scratch/"+str(args[0].split("/")[-1])+"_"+str(chromosome)+".tmp","w")
outputfile = str(args[0]).replace(".vcf","_"+str(chromosome)+".flt.vcf")   

print strandbias1
print strandbias2
print mapbias 
print enddistbias2 
print enddistbias2 
print nonreffreq
print minaltsnp 
print minaltindel 
print mindepth 
print maxdepth
print minsnpmapq 
print minindelmapq 
print variantqual 
print snpgapwin
print gapgapwin
print snpcluster
print snpclusterwindow
print chromosome
print outputfile

def IsIndel(info_list):
	check = []
	check = info_list

	if len(check[3]) > 1 or len(check[4].split(",")[0]) > 1:
        	return "INDEL"
        else:
                return  "SNP"


#Applying local filters

for line in input:
	if line.startswith("#"):
        	tempfile.write(str(line))
                continue

	if not line.startswith(str(chromosome)):
		continue

	#print line

	v = []
	v = line.strip("\n").split("\t")

	tag = ""

	if  v[3].lower() == v[4].lower():
		tempfile.write(str(line.strip("\n"))+"\t"+"NON VARIANT\n")
		continue
	
	info = []
	depth_info = []
	infodict = {}

	vartype = ""
	vartype = IsIndel(v)
	
	info = v[7].split(";")
	depth_info = (v[9].split(":"))[1].split(",")
		
	if v[3] == "N":
		tag = tag+"REFERENCE N;"
	
	if float(v[5]) < variantqual:
		tag = tag+"LOW QUAL VARIANT;"
	
	for i in info:
		if "=" in i:
			temp = []
			temp = i.split("=")
			infodict[temp[0]] = temp[1]

	if "," in infodict["AF"]:
        	tag = tag +"MULTIPLE ALLELES;"
		tempfile.write(str(line).strip("\n")+"\t"+str(tag)+"\n")
		continue		
                        	
	if vartype == "INDEL":
		if float(infodict["AF"]) < nonreffreq:
			tag = tag +"NON HOMOZYGOUS;"
		if int(depth_info[1]) < minaltindel:
			tag = tag+"LOW NON-REF READS;"
		if int(infodict["DP"]) < mindepth:
			tag = tag+"LOW DEPTH;"
		if float(infodict["MQ"]) < minindelmapq:
			tag = tag+"LOW MAPPING QUAL;"
		if int(infodict["DP"]) > maxdepth:
			tag = tag+"TOO MANY READS;"
		if float(infodict["QD"]) <= 2.0:
			tag = tag+"LOW QUALITY BY DEPTH;" 

		if float(infodict["FS"]) >  strandbias2:
			tag = tag+"STRAND BIAS;"
		if "ReadPosRankSum" in infodict:
			if float(infodict["ReadPosRankSum"]) < enddistbias2:
 				tag = tag + "END DIST BIAS;"
				
	if vartype == "SNP":
		if float(infodict["AF"]) < nonreffreq:
                        tag = tag +"NON HOMOZYGOUS;"
		if int(depth_info[1]) < minaltsnp:
                	tag = tag+"LOW NON-REF READS;"
		if float(int(depth_info[0]))/float(int(depth_info[0])+int(depth_info[1])) >= .20:
			tag = tag+"HIGH REF READS RATIO;" 
               	if int(infodict["DP"]) < mindepth:
                	tag = tag+"LOW DEPTH;"
                if float(infodict["MQ"]) < minsnpmapq:
                	tag = tag+"LOW MAPPING QUAL;"
                if int(infodict["DP"]) > maxdepth:
                	tag = tag+"TOO MANY READS;"
		if float(infodict["QD"]) <= 2.0:
                        tag = tag+"LOW QUALITY BY DEPTH;"
		
		if float(infodict["FS"]) >  strandbias1:
                        tag = tag+"STRAND BIAS;"
                if "ReadPosRankSum" in infodict:
                        if float(infodict["ReadPosRankSum"]) < enddistbias1:
                                tag = tag + "END DIST BIAS;"
		if "MQRankSum" in infodict:
                        if float(infodict["MQRankSum"]) < mapbias:
                                tag = tag + "MAP BIAS";
		
	if tag != "":
		tag = "LOCAL FLT;"+tag
	else:
		tag = "PASS"+tag

	tempfile.write(str(line).strip("\n")+"\t"+str(tag)+"\n")

tempfile.close()

print "Done with local filtering"

#Applying global filter

def split_line(line):
	vline = ""
	vline = line
	return  vline.strip("\n").split("\t")

def select_indel(indel1pos,indel2pos,info1,info2):
	info1_temp = []
	info2_temp = []
	ipos1 = indel1pos
	ipos2 = indel2pos
	
	info1_temp = info1.split(":")
        info2_temp = info2.split(":")

	if int(info1_temp[1].split(",")[1]) < int(info2_temp[1].split(",")[1]):
		return "WINDOW FLT; ADJ HQ IN-DEL","PASS"

	elif int(info1_temp[1].split(",")[1]) ==  int(info2_temp[1].split(",")[1]):	
			if info1_temp[2] < info2_temp[2]:
				return "PASS","WINDOW FLT; ADJ HQ IN-DEL"
			else:
				return "WINDOW FLT; ADJ HQ IN-DEL","PASS"
	else:
		return "PASS","WINDOW FLT; ADJ HQ IN-DEL"

inputfile = open("/scratch/"+str(args[0].split("/")[-1])+"_"+str(chromosome)+".tmp", 'r').read()
input = []
input = inputfile.split("\n")

output = open(outputfile,"w")

VCF = []
for a in range(0,len(input)-1):
     
        v = []
        v = input[a].strip("\n").split("\t")	
	VCF.append(v)

#INDEL SELECTION

for i in range(0,len(VCF)):
	
	if VCF[i][0].startswith("#"):
                continue
		
        if IsIndel(VCF[i]) == "INDEL" and VCF[i][-1] == "PASS":
		for j in range(i+1,len(VCF)):
                        if VCF[j][-1] == "PASS" and IsIndel(VCF[j]) == "INDEL":
				if int(VCF[i][1])+int(gapgapwin) > int(VCF[j][1]) and int(VCF[j][1]) > int(VCF[i][1]):
					VCF[i][-1],VCF[j][-1]= select_indel(VCF[i][1],VCF[j][1],VCF[i][9],VCF[j][9])
					break
				else:
					break

# FILTRATION OF SNPs ADJACENT TO INDEL

for m in range(0,len(VCF)):

        if VCF[m][0].startswith("#"):
                continue

	#VCF[m][6] = VCF[m][-1]

        if IsIndel(VCF[m]) == "INDEL" and VCF[m][-1] == "PASS":
		for n in range(m+1,len(VCF)):
                	if VCF[n][-1] == "PASS":
                        	if int(VCF[m][1])+int(snpgapwin) > int(VCF[n][1]) and int(VCF[n][1]) > int(VCF[m][1]):
					VCF[n][-1] = "WINDOW FLT; SNP ADJ HQ IN-DEL"
				else:
					break

		for o in range(m-1,0,-1):
			if VCF[o][-1] == "PASS":
                                if int(VCF[m][1]) < int(VCF[o][1])+int(snpgapwin) and int(VCF[m][1]) > int(VCF[o][1]):
                                        VCF[o][-1] = "WINDOW FLT; SNP ADJ HQ IN-DEL"
                                else:
                                        break

# SNP CLUSTER FILTRATION (3 SNPS FOUND WITHIN 10 BP)

#snpcluster = 3
#snpclusterwindow = 10

for p in range(0,len(VCF)):

        if VCF[p][0].startswith("#"):
                continue
	
	count = 0
	
        if IsIndel(VCF[p]) != "INDEL" and VCF[p][-1] == "PASS":
		count = 1
                for q in range(p+1,len(VCF)):
                        if VCF[q][-1] == "PASS":
				count = count + 1
			if count == snpcluster:
				if int(VCF[q][1]) > int(VCF[p][1]) + snpclusterwindow:
					break
				else:
					for z in range(p,q+1):
						if VCF[z][-1] == "PASS":
							VCF[z][-1] = "WINDOW FLT; SNP CLUSTER"

for x in VCF:
	if x[0].startswith("#"):
		output.write("\t".join(x)+"\n")
	else:
		out = []
		out = x
		out[6] = out[-1]
		if IsIndel(out) == "INDEL":
			out[7] = "INDEL;"+out[7]
		output.write(str("\t".join(out[0:-1]))+"\n")

output.close()

os.system("rm /scratch/"+str(args[0].split("/")[-1])+"_"+str(chromosome)+".tmp")
