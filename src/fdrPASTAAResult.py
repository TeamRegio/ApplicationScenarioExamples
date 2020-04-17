import sys, os
import numpy as np
import statsmodels.stats.multitest


if (len(sys.argv) < 4):
	print("python3 fdrPASTAAResult.py PASTAA_result, outputFile, pvalue")
else:
	PASTAA_result = sys.argv[1]
	output_ = sys.argv[2]
	p = float(sys.argv[3])


	TFs = []
	pvalues = []
	#read result
	with open(PASTAA_result, 'r') as result:
		for line in result:
			line = line.strip().split('\t')
			TFs.append(line[0])	
			pvalues.append(float(line[1]))
	
	#determine fdr

	rec, cor_pvalue = statsmodels.stats.multitest.fdrcorrection(pvals = pvalues, alpha = p,is_sorted = True)
#	print(rec)
#	print(cor_pvalue)
	counter = 0
	with open(output_, 'w') as o:
		o.write("TF\tpvalue(fdr correction)\n")
		for i in rec:
			if i == True:
				o.write(TFs[counter] + '\t' + str(cor_pvalue[counter]) + '\n')
			counter+=1
	
