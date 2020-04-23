import sys, os

if (len(sys.argv))< 5:
    print("Usage python ./rankingPASTAA.py  csv_file, fasta_file, outputRanking, outputFasta")

else:

	ranking = sys.argv[1]
	fastaFile = sys.argv[2]
	output = sys.argv[3]
	outputFasta = sys.argv[4]

	#identify sequences longer than 21
	headers = [] #stores all headers associated with a sequence longer or equal than 21
	with open(fastaFile, 'r') as f, open(outputFasta, 'w') as o:
		for line in f:
			if line[0] == ">":
				if "::" in line:
					line = line.split("::")
					header = line[0] + '\n'
				else:
					header = line
			else:
				line = line.strip()
				if len(line) >= 21:
					o.write(header + line + '\n')
					headers.append(header[1:-1])

	#parse ranking REMs in suitable format for PASTAA
	counter = 1
	rankedCSV= []
	# rank csv file according to their absoult value
	with open(ranking, 'r') as r:
		r.readline() #skip header
		for line in r:
			line = line.strip().split('\t')
			if line[2] in headers:
				rankedCSV.append([line[2], abs(float(line[10]))])
	rankedCSV.sort(key = lambda x: x[1], reverse = True) 
	#print("rankedCSV: " + str(rankedCSV))

	with  open(output, 'w') as o:
		for i in rankedCSV:
			o.write(i[0] + '\t' + str(counter) + '\n')
			counter+=1
