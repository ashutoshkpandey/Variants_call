import sys
import re, fileinput
 
Argument = []
Argument = sys.argv[1:]  # Reading the arguments

"""
Usage: Vcf file, Big file from Pindel (with scores), Outputfile
Opening the file and storing the lines in list
"""

Filepath = Argument[0]
Filepath1 = Argument[1] 



def extract_BP_SCORE(row_temp):
	score = 0
	bp = 0
	rtemp = ""

	rtemp = row_temp
	ret = []

	info_list = []
	info_list = rtemp.split("\t")

	for a in range(0,len(info_list)-1):
		#print info_list[a]
		if info_list[a].startswith("ChrID"):
                        ret.append(info_list[a].split(" ")[1])
		if info_list[a].startswith("BP "):
			ret.append(int(info_list[a].split(" ")[1])+1)
		if info_list[a].startswith("S1"):
                        ret.append(info_list[a].split(" ")[1])	
	
	#print ret
	return ret

Scores = {}

output = open(Argument[2],"w")

for row in fileinput.input([Filepath]):
        array = []
        if not row[0].isdigit():
		continue

        array = extract_BP_SCORE(row)
	
	if array[0]+"\t"+str(array[1]) not in Scores:
		Scores[array[0]+"\t"+str(array[1])] = array[2]
	

for row1 in fileinput.input([Filepath1]):
        array1 = []
        if row1.startswith("#"):
		output.write(str(row1))
                continue

        array1 = row1.split("\t")
	#print int(Scores[array1[0]+"\t"+array1[1]])

	if array1[0]+"\t"+array1[1] in Scores:
		if int(Scores[array1[0]+"\t"+array1[1]]) >= 25:
			array2 = []
			array2 = array1[7].split(";")
			if int(array2[2].lstrip("SVLEN=")) >= 100:
				output.write(str(row1))
	else:
		print array1[0]+"\t"+array1[1] 

output.close()
