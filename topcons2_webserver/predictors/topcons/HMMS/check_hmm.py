import os


inPath = "topcons.hmg.new"#"/Users/christoph/Documents/phd/project/TopoPredOpt/input/hmm/firstlast_pkrinout_140327.hmg"
outPath = "topcons.hmg.fixed"#"/Users/christoph/Documents/phd/project/TopoPredOpt/input/hmm/fixed_firstlast_pkrinout_140327.hmg"

dict_vertices = {}
hmmOrig = []
hmmList = []
hmmListNew = []

linebreak = ""

with open(inPath) as hmm_file:
	for line in hmm_file:
		hmmOrig.append(line)
		linebreak = line[-1]


for line in hmmOrig:
	if "Vertex " in line and not "	Vertex " in line and not "label" in line and not "type" in line:
		line =  "Vertex " + " -" + line.replace(":","").split()[1] + ":" + linebreak
	elif "	Vertex " in line:# and not "Vertex " in line and not "label" in line and not "type" in line:
		line = "	Vertex " + " -" + line.replace(":","").split()[1] + ": 1.0" + linebreak
	hmmList.append(line)
iCounter = 0

for line in hmmList:
	if "Vertex " in line and not "	Vertex " in line and not "label" in line and not "type" in line:
		dict_vertices[line.replace(":","").split()[1]] = str(iCounter)
		iCounter += 1

for line in hmmList:
	if "Vertex " in line and not "	Vertex " in line and not "label" in line and not "type" in line:
		line =  "Vertex " + dict_vertices[line.replace(":","").split()[1]] + ":" + linebreak
	if "	Vertex " in line :#and not "Vertex " in line and not "label" in line and not "type" in line:
		line = "	Vertex " + dict_vertices[line.replace(":","").split()[1]] + ": 1.0" + linebreak
		print line
	hmmListNew.append(line)


for i in range(0, len(hmmList)):
	if (hmmOrig[i] == hmmListNew[i]) == False:
		print hmmOrig[i]
		print hmmListNew[i]
		print "---"

with open(outPath, "w") as outFile:
	for line in hmmListNew:
		outFile.write(line)
