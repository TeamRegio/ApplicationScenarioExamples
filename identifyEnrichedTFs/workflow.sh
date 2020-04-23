#!/bin/bash

#input 
#note al path have to end with /
setOfMotifs=$1 #set of motifs in transfac format (JASPAR2020_HUMAN_transfac_P0.txt)
pathToRepo=$2 #path to cloned gitHub repo
genome=$3 # genome in fasta format
REMs=$4 #CSV file from EpiRegio
outputDir=$5 #user defined output dir
pvalue=$6 #pvalue threshold for PASTAA


#mkdir -p ${outputDir}
##Step 1: change CSV file format to bed file format  and run bedtools getFasta function
echo "Step1: determine DNA sequence"
awk 'NR!=1{print $4 "\t" $5 "\t" $6 "\t" $3}' ${REMs}  >${outputDir}REMs.bed #remove first line (header) and reorder the columns

bedtools getfasta -name -fi ${genome} -bed ${outputDir}REMs.bed -fo ${outputDir}REMs_.fa #determine fasta seq

##Step 2: parse REM ranking in correct format for PASTAA and delete REMs shorter than the longest moitf (21bp)
echo "Step 2: parse REM ranking in the format PASTAA requires"
python3 ${pathToRepo}src/rankingPASTAA.py  ${REMs} ${outputDir}REMs_.fa ${outputDir}ranking.txt ${outputDir}REMs.fa

rm ${outputDir}REMs_.fa

##Step 3: call PASTAA workflow
echo "run PASTAA workflow"
#transform PFMs to PSEM (energy file)
${pathToRepo}src/PSCM_to_PSEM ${setOfMotifs} >${outputDir}/energy.txt
echo "determined energy file"
#run TRAP
${pathToRepo}src/TRAP_normalized ${outputDir}/energy.txt ${outputDir}REMs.fa >${outputDir}/TRAP_output.txt
echo "determined TRAP"
#run PASTAA
echo "${pathToRepo}src/PASTAA ${outputDir}/TRAP_output.txt ${outputDir}ranking.txt |  sort -k2,2 -g  > ${outputDir}/PASTAA_output.txt"
${pathToRepo}src/PASTAA ${outputDir}/TRAP_output.txt ${outputDir}ranking.txt |  sort -k2,2 -g  > ${outputDir}/PASTAA_output.txt
echo "determined PASTAA"

#Step 4: fdr correction PASTAA result
python3 ${pathToRepo}src/fdrPASTAAResult.py ${outputDir}PASTAA_output.txt ${outputDir}PASTAA_result.txt ${pvalue}
