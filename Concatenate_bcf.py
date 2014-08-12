import os, glob
import sys,re,fileinput

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 2:	
	print "Usage: Input_directory Chromosome_size(fai file)" 
	sys.exit()
  							
chr_order = ""

for line in fileinput.input([Argument[1]]):
        if line.startswith("#") or not line.strip():
                continue

        array = []
        line = line.rstrip("\n")
        array = line.split("\t")

 	chr_order = chr_order+"  "+Argument[0]+"/"+array[0]+".var.raw.bcf"
	

print ("/home/apandey/bio/samtools-0.1.18/bcftools/bcftools cat "+str(chr_order)+" > "+str(Argument[0])+"/Whole_genome_variants.bcf")
	
os.system("/home/apandey/bio/samtools-0.1.18/bcftools/bcftools cat "+str(chr_order)+" > "+str(Argument[0])+"/Whole_genome_variants.bcf")



