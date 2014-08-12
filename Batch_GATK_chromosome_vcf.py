import os, glob
import sys,re,fileinput

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 2:	
	print "Usage: Input_file Job_directory Chromosome_size(fai file)" 
	sys.exit()
  	

if not os.path.exists(str(Argument[1])):
        os.makedirs(str(Argument[1]))
				
for line in fileinput.input([Argument[2]]):
        if line.startswith("#") or not line.strip():
                continue

        array = []
        line = line.rstrip("\n")
        array = line.split("\t")
        
	jobname = ""
        jobname = str(Argument[1])+"/"+re.sub(r'\s','',str(array[0])+"_vcf.sh")
        jobfile = open(str(jobname),"w")

	jobfile.write("#!/bin/bash\n#PBS -l walltime=240:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -N Filtering_variants_"+str(array[0])+"\n\n/home/apandey/Scripts/NGS_Scripts/Variants_call/Filter_GATK_vcf.py  "+str(Argument[0])+"   "+str(array[0])+"\n")

        jobfile.close()

        print "qsub "+str(jobname)
        os.system("qsub "+str(jobname))
