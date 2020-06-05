import matplotlib.pyplot as plt


testDataIn=open("testOutData.txt",'r')
trueDataIn=open("trueOutData.txt",'r')

fig=plt.figure()
plotData=[[0 for i in range(5)] for throw in range(6)]

for line in range(5):
	pullIn=testDataIn.readline()
	timePoint, inData=pullIn.split(" ")[0],pullIn.split(" ")[1:]
	plotData[0][line]=float(timePoint)
	inData.pop()
	for index, value in enumerate(inData):
		plotData[index+1][line]=(float(value))
		
for line in range(0,5):
	plt.plot(plotData[0],plotData[line], label="Fit Data "+str(line))


for line in range(5):
	pullIn=trueDataIn.readline()
	timePoint, inData=pullIn.split(" ")[0],pullIn.split(" ")[1:]
	plotData[0][line]=float(timePoint)
	inData.pop()
	for index, value in enumerate(inData):
		plotData[index+1][line]=(float(value))
		
for line in range(0,5):
	plt.plot(plotData[0],plotData[line],label="True Data "+str(line))

plt.ylim(0,225)
plt.legend()
plt.show()