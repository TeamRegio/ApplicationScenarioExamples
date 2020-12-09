import sys, os
from statsmodels.stats.multitest import fdrcorrection

#global variables
PVALUE = 0.01
WPVALUE = 0.001


#if (len(sys.argv) < 6):
#	print("python3 parseOutputFimo.py FIMO_output, epiregio_file,PASTAA_output, output_info, output, fimoTSS.txt")
if (len(sys.argv) < 5):
	print("python3 parseOutputFimo.py FIMO_output, epiregio_file,PASTAA_output, output_info, output")
else:
	FIMO = sys.argv[1]
	epiregio = sys.argv[2]
	PASTAA_file = sys.argv[3]
	output = open(sys.argv[4],'w')
	#write header output
	output.write("TF\tREMID\tregion\tcellTypeScore\tREM_activity(DNase_signal)\tstart\tstop\tstart_relative\tend_relative\tstrand\tmatched_seq\n")
	output2 = open(sys.argv[5], 'w')
	#fimoTSS = sys.argv[6]

	#read PASTAA file and store all significant motifs
	sig_TFs = []
	with open(PASTAA_file, 'r') as result:
		for line in result:
			line = line.strip().split('\t')
			if float(line[1]) <= PVALUE:
				sig_TFs.append(line[0])
	
	print("sig TFs: " + str(len(sig_TFs)))
	

	#read activity and cellType score from epiregio
	REMID_mapping = {} #chr:start-end -> REMID
	cellTypeScore = {} #REMID->score 
	dnaseSignal = {} #REMID-> signal
	with open(epiregio, "r") as input_:
		input_.readline()
		for line in input_:
			line = line.strip().split('\t')
			key = line[3] + "\t" + line[4] + "\t" + line[5]
			id_ = line[2]
			REMID_mapping[id_] = key
			cellTypeScore[id_] = line[10]
			dnaseSignal[id_] = line[11]
	#print("REM_ID_mapping: " + str(REMID_mapping))	
	#print("cellTypeScore: " + str(cellTypeScore))	
	#print("dnaseSignal: " + str(dnaseSignal))	

	#read FIMO 
	hits_per_TF = {} #{TF -> {REMID -> numberHits}}
	#prefill hits_per_TF
	for i in sig_TFs:
		helper = {}
		for j in cellTypeScore.keys():
			helper[j] = 0
		hits_per_TF[i] = helper

	TFs_Fimo = []
	REMs_Fimo = []
	with open(FIMO, 'r') as result:
		result.readline()
		for line in result:
			if not line.strip() or line[0] == "#":
				break
			line = line.strip().split('\t')
			currentTF = line[1]
			#check if TF is a significant one
			if currentTF in sig_TFs:
		#		print(currentTF)
				if currentTF not in TFs_Fimo:
					TFs_Fimo.append(currentTF)
				helper = hits_per_TF[currentTF]
				
				currentREMID = line[2]
				if currentREMID not in REMs_Fimo:
					REMs_Fimo.append(currentREMID)
				helper[currentREMID] = helper[currentREMID] + 1
	 			#write info file FIMO result in such a way that it is suitable for james
				region = REMID_mapping[line[2]]
				region = region.strip().split('\t')
				output.write(line[1] + "\t" + line[2] + '\t' + region[0] + ":" + region[1] + "-" + region[2] + '\t' + cellTypeScore[line[2]] + '\t' + dnaseSignal[line[2]]+ '\t'  + line[3] + '\t' + line[4]  + '\t' + str(int(region[1]) + int(line[3])) + '\t' + str(int(region[1]) + int(line[4]))  +  '\t' + line[5] + '\t' + line[9] + '\n')


	#read FIMO TSS file
	#TSS_long = {} # TF -> binding site counts
	#TSS_short = {} # same just for the short isoform
	#with open(fimoTSS, 'r') as f:
#		f.readline() #skip header
#		for line in f:
#			line = line.strip().split('\t')
#			print(line)
#			TF = line[1]
#			seq_name = line[2].split("_")
#			seq_name = seq_name[2]
#			if seq_name == "shortIsoform":
#				if TF in TSS_short.keys():
##					TSS_short[TF] = TSS_short[TF] + 1
#				else:
#					TSS_short[TF] = 1
#			else:
#				if TF in TSS_long.keys():
#					TSS_long[TF] = TSS_long[TF] + 1
#				else:
#					TSS_long[TF] = 1
	#write output
	output2.write("TF")
	order = []
	for k in cellTypeScore.keys():
		if k in REMs_Fimo:
			output2.write('\t' + k)
			order.append(k)
	#output2.write( '\tTSS_short\tTSS_long\n')
	output2.write( '\n')

	TFs = hits_per_TF.keys()
	#print("len TFs: " + str(len(TFs)))

	#skip rows with zeros
	for TF in TFs:
		if TF in TFs_Fimo: #check if TF occurs in FIMO result
			output2.write(TF)
			REMs = hits_per_TF[TF]
		#	print(REMs)
			#for elem in REMs.keys():
			for elem in order:
				#if elem in REMs_Fimo: #check if REM is considered in FIMO result
					output2.write("\t" +str(REMs[elem]))
			output2.write('\n')
#			if TF in TSS_short.keys():
#				output2.write("\t" + str(TSS_short[TF]))
#			else:
#				output2.write("\t0")
#			if TF in TSS_long.keys():
#				output2.write("\t" + str(TSS_long[TF])+ '\n')
#			else:
#				output2.write("\t0\n")
	print(REMs_Fimo)
	output2.write("activity")
	#for elem in REMs_Fimo:
	for elem in order:
		output2.write('\t' + dnaseSignal[elem])
	#output2.write("\t0.0\t0.0\ncellTypeScore")
	output2.write("\ncellTypeScore")
	#for elem in REMs_Fimo:
	for elem in order:
		output2.write('\t' + cellTypeScore[elem])
	#output2.write('\t0.0\t0.0\n')
	output2.write('\n')

	output.close()
	output2.close()
