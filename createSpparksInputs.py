inFile=open("candidateParameterSets\\nextParameters.txt",'r')
numberOfFiles=int(inFile.readline())

for file in range(0,numberOfFiles):
	outFile=open("candidateParameterSets\\candidateSet_"+str(file)+".txt",'w')
	parameterSet=inFile.readline().split(",")
	parameterSet.pop()
	
	outFile.write("seed	12345\n")
	outFile.write("app_style	chemistry\n")
	outFile.write("solve_style	linear\n")
	outFile.write("volume	1\n")
	
	speciesList=['G','GS','M','P','D']
	for species in speciesList:
		outFile.write("add_species 	"+species+"\n")
	
	outFile.write("add_reaction	1 G "+parameterSet[0]+" GS\n")
	outFile.write("add_reaction	2 GS "+parameterSet[1]+" G\n")
	outFile.write("add_reaction	3 GS "+parameterSet[2]+" GS M\n")
	outFile.write("add_reaction	4 M "+parameterSet[3]+" D\n")
	outFile.write("add_reaction	5 M "+parameterSet[4]+" M P\n")
	outFile.write("add_reaction	6 P "+parameterSet[5]+" D\n")
	
	outFile.write("count	G 1\n")
	outFile.write("count	GS 0\n")
	outFile.write("count	M 0\n")
	outFile.write("count	P 0\n")
	outFile.write("count	D 0\n")
	
	outFile.write("run	120\n")
	outFile.close()
	
