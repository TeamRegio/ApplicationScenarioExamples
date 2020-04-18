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
				line = line.split("::")
				header = line[0] + '\n'
			else:
				line = line.strip()
				if len(line) >= 21:
					o.write(header + line + '\n')
					headers.append(header)

	#parse ranking REMs in suitable format for PASTAA
	counter = 1
	with open(ranking, 'r') as r, open(output, 'w') as o: #assumes that csv file is already sorted
		r.readline() #skip header
		for line in r:
			line = line.strip().split('\t')
			key = line[3] + ":" + line[4] + "-" + line[5] + "_" + line[0] #with REM Id to be unique 
			o.write(key + '\t' + str(counter) + '\n')
			counter+=1
