import sys, os

cutoff = 0

if (len(sys.argv))< 5:
    print("Usage python ./rankingPASTAA.py  csv_file, fasta_file, outputRanking, outputFasta")

else:

	ranking = sys.argv[1]
	fastaFile = sys.argv[2]
	output = sys.argv[3]
	outputFasta = sys.argv[4]

	#parse ranking REMs in suitable format for PASTAA
	counter = 1
	rankedCSV= []
	rankedCSV_helper = {}
	# rank csv file according to their absoult value
	with open(ranking, 'r') as r:
		r.readline() #skip header
		for line in r:
			line = line.strip().split('\t')
			rankedCSV_helper[line[2]] =  abs(float(line[10]))

	#identify sequences longer than 24
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
				#if len(line) >= 24 and rankedCSV_helper[header[1:-1]] >= cutoff:
				if len(line) >= 24:
					rankedCSV.append([header[1:-1], rankedCSV_helper[header[1:-1]]])
					o.write(header + line + '\n')
					headers.append(header[1:-1])

	rankedCSV.sort(key = lambda x: x[1], reverse = True) 
#	print("rankedCSV: " + str(rankedCSV))

	with  open(output, 'w') as o:
		for i in rankedCSV:
			o.write(i[0] + '\t' + str(counter) + '\n')
			counter+=1
