import sys
import re, fileinput
 
Argument = []
Argument = sys.argv[1:]  # Reading the arguments

"""
Opening the file and storing the lines in list
"""

Filepath = Argument[0]
Filepath1 = Argument[1] 

count1 = 0
count2 = 0
countc = 0

countcpp = 0
countcfp = 0
countcpf = 0
countcff = 0
countup = 0
countuf = 0
counta = 0
countb = 0 

Variants = {}

output = open(Argument[2],"w")

for row in fileinput.input([Filepath]):
        array = []
        if row.startswith("#"):
        	continue
        array = row.split("\t")
	
	if array[0]+"\t"+array[1] not in Variants:
		Variants[array[0]+"\t"+array[1]] = array[6].strip()
		count1 = count1 + 1
	else:
		counta = counta + 1
		if array[6].strip() == "PASS":
			Variants[array[0]+"\t"+array[1]] = array[6].strip()

Done = {}
STAT = ""

for row1 in fileinput.input([Filepath1]):
        array1 = []
        if row1.startswith("#"):
		output.write(str(row1))
                continue

        array1 = row1.split("\t")
	count2 = count2 + 1

	if  array1[0]+"\t"+array1[1] in Done:
		if STAT != "PASS":
			if array1[6].strip() == "PASS":
				array1[7] = array1[7]+";"+"DUAL_VARIANT_PASS_FAIL"
			else:
				array1[7] = array1[7]+";"+"DUAL_VARIANT_FAIL_FAIL"
		else:
			if array1[6].strip() == "PASS":
                                array1[7] = array1[7]+";"+"DUAL_VARIANT_PASS_PASS"
                        else:
                                array1[7] = array1[7]+";"+"DUAL_VARIANT_FAIL_PASS"
		countb = countb + 1
	

	FLAG = ""
	TEMP =[]
	TEMP = (array1[6].rstrip(";")).split(";")

	if TEMP[0] == "PASS":
		FLAG = "PASS"
	else:	
		for x in TEMP[1:]:
			FLAG = FLAG + str((x.strip(" ")).replace (" ", "_")) +","
        	FLAG = str("-".join(TEMP[0].split(" ")))+"="+str(FLAG).rstrip(",")
	
	#print FLAG
	STAT = FLAG

	if array1[0]+"\t"+array1[1] in Variants:
		countc = countc + 1

		if Variants[array1[0]+"\t"+array1[1]] == "PASS":
			if FLAG == "PASS":
				array1[7] = array1[7]+";"+str(FLAG)+";"+"COMMON_PASS_PASS"
				array1[6] = "."
				countcpp = countcpp + 1
			else:
				array1[7] = array1[7]+";"+str(FLAG)+";"+"COMMON_FAIL_PASS"
				array1[6] = "."
				countcfp = countcfp + 1
		else:
			if FLAG == "PASS":
				array1[7] = array1[7]+";"+str(FLAG)+";"+"COMMON_PASS_FAIL"
				array1[6] = "."
				countcpf = countcpf + 1
			else:
				array1[7] = array1[7]+";"+str(FLAG)+";"+"COMMON_FAIL_FAIL"
                                array1[6] = "."
				countcff = countcff +  1
	else:
		if FLAG == "PASS":
			array1[7] = array1[7]+";"+str(FLAG)+";"+"UNIQUE_PASS"
			array1[6] = "."
			countup = countup + 1
		else:
			array1[7] = array1[7]+";"+str(FLAG)+";"+"UNIQUE_FAIL"
                        array1[6] = "."
			countuf = countuf + 1

	Done[array1[0]+"\t"+array1[1]] = "Yes"

	output.write(str("\t".join(array1)))
	output.flush()		

print "# of Variants in File1", count1
print "# of Variants in File2", count2
print "# of Variants in common", countc
print "# of Variants in common and both passed", countcpp
print "# of Variants in common, failed in File 2 but  passed in 1", countcfp
print "# of Variants in common,failed in File 1 but  passed in 2",countcpf
print "# of Variants in ccmmon and failed in both", countcff
print "# of Variants unique to File2 and passed", countup
print "# of Variants unique to File2 and failed", countuf
print "# of SNPs and Indels at same position in File1", counta
print "# of SNPs and Indels at same position in File2", countb
output.close()
   	 
