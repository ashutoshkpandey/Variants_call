import os, glob
import sys,re,fileinput

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 4:	
	print "Usage: List of Bam_files_directories(csv) Job_Script_directory Chromosome_size(file) Output_directory" 
	sys.exit()
  
def dir(patharray):
    Listdir = []
    for infile in patharray:
        Listdir.append(os.path.join(idpath,infile))
    return Listdir    
									
if not os.path.exists(str(Argument[1])):
	os.makedirs(str(Argument[1]))

if not os.path.exists(str(Argument[3])):
        os.makedirs(str(Argument[3]))

dpath = []
dpath = Argument[0].split(",")

saminput = ""
BAM = []

for idpath in dpath:
	Listoffile = []
	Listoffile = dir(os.listdir(idpath))
	#print Listoffile

	for file in Listoffile:
        	if file.endswith(".bam") and not os.path.getsize(file) == 0:
                	if file not in BAM:
				#print file
                                saminput = saminput+"  "+file
				BAM.append(file)

#print saminput

Chrsize = {}

for line in fileinput.input([Argument[2]]):
        if line.startswith("#") or not line.strip():
                continue

        array = []
        line = line.rstrip("\n")
        array = line.split("\t")

        if array[0] not in Chrsize:
		Chrsize[array[0]] = array[0]+":1-"+array[1]


for region in Chrsize:

	jobname = ""
	jobname = str(Argument[1])+"/"+re.sub(r'\s','',str(region)+"_bcf.sh")
	jobfile = open(str(jobname),"w")
	
	jobfile.write("#!/bin/bash\n#PBS -l walltime=240:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -N Samtools_variant_calling_"+str(region)+"\n\n/home/apandey/bio/samtools-0.1.18/samtools mpileup -C 50 -d 1000 -E -Q 20 -q 17 -uf /home/apandey/Reference_Fasta/mm10/mm10_ucsc/mm10_ucsc.fa -r "+str(Chrsize[region])+" "+str(saminput)+"  | /home/apandey/bio/samtools-0.1.18/bcftools/bcftools view -t .0001 -bvcN - > "+Argument[3]+"/"+str(region)+".var.raw.bcf\n")

	jobfile.close()
	
	print "qsub "+str(jobname)
	os.system("qsub "+str(jobname))



