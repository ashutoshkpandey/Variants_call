#!/usr/bin/python

# Local Filters

# SB Strand Bias
# BB Base Qual Bias
# EB End Distance Bias
# MB MapQ Bias
# VDB Variant Qual Bias
# AF1 EM estimate of the site allele frequency of the strongest non-reference allele
# AB Minimum non-reference supporting reads for the SNP call
# ABI Minimum non-reference supporting reads for the Indel call
# d low depth
# D high depth
# N Reference is N
# QS min RMS mapping quality (SNP only)
# QI min RMS mapping quality (Indel only)
# Qual Variant calling quality 

# Global Filters 

# W Adjacent Indel distance (Indel only) (if two gaps or indels lie within W (default 25 bp), indel with maximum supporting reads will be marked as PASS. In case of tie, Indel with higher ratio of # of non-reference allele supporting reads to total depth will be preferred.
# w Adjacent SNP distance between Indel and SNP to remove SNPs.


import sys,os
import getopt

def usage():
	print '''usage: varfilter.py [options] input_vcf 

Options: 	 -SB Float minimum p-value for strand bias
		 -BB Float minimum p-value for base quality bias
		 -EB Float minimum p-value for end distance  bias
		 -MB Float minimum p-value for MapQ bias 
		 -VDB Float minimum p-value for variant quality  bias
		 -AF1 Float minimum site allele frequency of the strongest non-reference allele
		 -ABS INT minimum alternate bases needed to call SNPs
		 -ABI INT minimum alternate bases needed to call Indels
                 -d INT minimum read depth
                 -D INT maximum read depth
                 -QS INT minimum RMS mapping quality for SNPs 
		 -QI INT minimum RMS mapping quality for Indels
		 -Qual  INT minimum variant calling quality
		 		
		 -w INT window size for filtering SNPs adjacent to Indel
		 -W  INT window size for filtering adjacent Indels'''

strandbias = .0001
basebias = 0.0
mapbias = 0.0
enddistbias = .0001
variantdistbias = 0.0
nonreffreq = .95
minaltsnp = 3
minaltindel = 3
mindepth = 5
maxdepth = 200
minsnpmapq = 25
minindelmapq = 20
variantqual = 20

gapgapwin = 20
snpgapwin = 10

try:
	options, args = getopt.getopt(sys.argv[1:], 'SB:BB:EB:MB:VDB:ABS:ABI:d:D:QS:QI:Qual:w:W:', [])
except getopt.GetoptError:
	usage()
	sys.exit()

for (oflag, oarg) in options:
	
	if oflag == '-SB': strandbias = float(oarg) 
	if oflag == '-BB': basebias = float(oarg)
	if oflag == '-MB': mapbias = float(oarg)
	if oflag == '-EB': enddistbias = float(oarg)
	if oflag == '-VDB': variantdistbias = float(oarg)
	if oflag == '-AF1': nonreffreq = float(oarg)
	if oflag == '-ABS': minaltsnp = int(oarg)
	if oflag == '-ABI': minaltindel = int(oarg)
	if oflag == '-d': mindepth = int(oarg)
	if oflag == '-D': maxdepth = int(oarg)
	if oflag == '-QS': minsnpmapq = int(oarg)
	if oflag == '-QI': minindelmapq = int(oarg)
	if oflag == '-Qual': variantqual = float(oarg)
	if oflag == '-w': snpgapwin = int(oarg)
	if oflag == '-W': gapgapwin = int(oarg)


if len(args) < 1:
	input = sys.stdin
else:
	input = open(args[0])

tempfile = open("/scratch/"+str(args[0].split("/")[-1])+".iitmp","w")
outputfile = str(args[0]).replace("var.raw.vcf","var.flt.vcf")   

print strandbias 
print basebias 
print mapbias 
print enddistbias 
print variantdistbias 
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
print outputfile


#Applying local filters

for line in input:
	if line.startswith("#"):
        	tempfile.write(str(line))
                continue

	v = []
	v = line.strip("\n").split("\t")

	tag = ""

	if  v[3].lower() == v[4].lower():
		tempfile.write(str(line.strip("\n"))+"\t"+"NON VARIANT\n")
		continue
	
	info = []
	infodict = {}
	vartype = ""

	if v[7].startswith("INDEL"):
		vartype = "INDEL"
	else:
		vartype = "SNP"
	
	info = v[7].split(";")
		
	if v[3] == "N":
		tag = tag+"REFERENCE N;"
	
	if float(v[5]) < variantqual:
		tag = tag+"LOW QUAL VARIANT;"
	
	for i in info:
		if "=" in i:
			temp = []
			temp = i.split("=")
			infodict[temp[0]] = temp[1]

	
	Parameters = ["DP4","DP","VDB","MQ"]
	
	for param in Parameters:
		if param not in infodict:
			infodict[param] = 1 


	if vartype == "INDEL":
		if float(infodict["AF1"]) < nonreffreq or float(infodict["FQ"]) > 0 :
			tag = tag +"NON HOMOZYGOUS;"
		if int(infodict["DP4"].split(",")[2])+int(infodict["DP4"].split(",")[3]) < minaltindel:
			tag = tag+"LOW NON-REF READS;"
		if int(infodict["DP"]) < mindepth:
			tag = tag+"LOW DEPTH;"
		if int(infodict["MQ"]) < minindelmapq:
			tag = tag+"LOW MAPPING QUAL;"
		if int(infodict["DP"]) > maxdepth:
			tag = tag+"TOO MANY READS;"
		if float(infodict["VDB"]) < variantdistbias:
			tag = tag + "VARIANT DIST BIAS;"
		if "PV4" in infodict:
			PV4 = []
			PV4 = infodict["PV4"].split(",")

			if float(PV4[0]) < strandbias:
				tag = tag+"STRAND BIAS;"
			if float(PV4[1]) < basebias:
				tag = tag + "BASE BIAS;"
			if float(PV4[2]) < mapbias:
				tag = tag + "MAP BIAS;"
			if float(PV4[3]) < enddistbias:
 				tag = tag + "END DIST BIAS;"
				
	if vartype == "SNP":
		if float(infodict["AF1"]) < nonreffreq or float(infodict["FQ"]) > 0 :
                        tag = tag +"NON HOMOZYGOUS;"
		if int(infodict["DP4"].split(",")[2])+int(infodict["DP4"].split(",")[3]) < minaltsnp:
                	tag = tag+"LOW NON-REF READS;"
               	if int(infodict["DP"]) < mindepth:
                	tag = tag+"LOW DEPTH;"
                if int(infodict["MQ"]) < minsnpmapq:
                	tag = tag+"LOW MAPPING QUAL;"
                if int(infodict["DP"]) > maxdepth:
                	tag = tag+"TOO MANY READS;"
                if float(infodict["VDB"]) < variantdistbias:
                        tag = tag + "VARIANT DIST BIAS;"
               	if "PV4" in infodict:
               		PV4 = []
                        PV4 = infodict["PV4"].split(",")

                        if float(PV4[0]) < strandbias:
                        	tag = tag+"STRAND BIAS;"
                        if float(PV4[1]) < basebias:
                        	tag = tag + "BASE BIAS;"
                        if float(PV4[2]) < mapbias:
                        	tag = tag + "MAP BIAS;"
                        if float(PV4[3]) < enddistbias:
                                tag = tag + "END DIST BIAS;"

	if tag != "":
		tag = "LOCAL FLT;"+tag
	else:
		tag = "PASS"+tag

	tempfile.write(str(line).strip("\n")+"\t"+str(tag)+"\n")

tempfile.close()

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
	
	info1_temp = info1.split(";")
        info2_temp = info2.split(";")

	info1_dict = {}
	info2_dict = {}

	for i in info1_temp:
                if "=" in i:
                	temp = []
                        temp = i.split("=")
                        info1_dict[temp[0]] = temp[1]	
	
	for i in info2_temp:
                if "=" in i:
                        temp = []
                        temp = i.split("=")
                        info2_dict[temp[0]] = temp[1]

	if int(info1_dict["DP4"].split(",")[2])+int(info1_dict["DP4"].split(",")[3]) < int(info2_dict["DP4"].split(",")[2])+int(info2_dict["DP4"].split(",")[3]):
		return "WINDOW FLT; ADJ HQ INDEL","PASS"

	elif int(info1_dict["DP4"].split(",")[2])+int(info1_dict["DP4"].split(",")[3]) == int(info2_dict["DP4"].split(",")[2])+int(info2_dict["DP4"].split(",")[3]):
			
			if info1_dict["MQ"] > info2_dict["MQ"]:
				return "PASS","WINDOW FLT; ADJ HQ INDEL"

			elif info1_dict["MQ"] == info2_dict["MQ"]:
				if info1_dict["DP"] > info2_dict["MQ"]:
					return "PASS","WINDOW FLT; ADJ HQ INDEL"
				else:
					return "WINDOW FLT; ADJ HQ INDEL","PASS" 
			else:
				return "WINDOW FLT; ADJ HQ INDEL","PASS"
	else:
		return "PASS","WINDOW FLT; ADJ HQ INDEL"

inputfile = open("/scratch/"+str(args[0].split("/")[-1])+".iitmp", 'r').read()
input = []
input = inputfile.split("\n")

output = open(outputfile,"w")

VCF = []
for a in range(0,len(input)-1):
     
        v = []
        v = input[a].strip("\n").split("\t")	
	VCF.append(v)

for i in range(0,len(VCF)):
	
	if VCF[i][0].startswith("#"):
                continue
		
        if VCF[i][7].startswith("INDEL") and VCF[i][-1] == "PASS":
		for j in range(i+1,len(VCF)):
                        if VCF[j][-1] == "PASS" and VCF[j][7].startswith("INDEL"):
				if int(VCF[i][1])+int(gapgapwin) > int(VCF[j][1]) and int(VCF[j][1]) > int(VCF[i][1]):
					VCF[i][-1],VCF[j][-1]= select_indel(VCF[i][1],VCF[j][1],VCF[i][7],VCF[j][7])
					break
				else:
					break


for m in range(0,len(VCF)):

        if VCF[m][0].startswith("#"):
                continue

	#VCF[m][6] = VCF[m][-1]

        if VCF[m][7].startswith("INDEL") and VCF[m][-1] == "PASS":
		for n in range(m+1,len(VCF)):
                	if VCF[n][-1] == "PASS":
                        	if int(VCF[m][1])+int(snpgapwin) > int(VCF[n][1]) and int(VCF[n][1]) > int(VCF[m][1]):
					VCF[n][-1] = "WINDOW FLT; SNP ADJ HQ INDEL"
				else:
					break

		for o in range(m-1,0,-1):
			if VCF[o][-1] == "PASS":
                                if int(VCF[m][1]) < int(VCF[o][1])+int(snpgapwin) and int(VCF[m][1]) > int(VCF[o][1]):
                                        VCF[o][-1] = "WINDOW FLT; SNP ADJ HQ INDEL"
                                else:
                                        break

for x in VCF:
	if x[0].startswith("#"):
		output.write("\t".join(x)+"\n")
	else:
		out = []
		out = x
		out[6] = out[-1]
		output.write(str("\t".join(out[0:-1]))+"\n")

output.close()

os.system("rm /scratch/"+str(args[0].split("/")[-1])+".iitmp")
