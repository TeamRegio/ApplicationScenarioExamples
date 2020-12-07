#!/bin/bash

EpiRegioOutput=$1 #e.g. PCAT19_EpiRegio_REMs.txt
outputDir=$2 #e.g. results/
outputDirPASTAA=$3 #e.g. PASTAA_result
genome=$4 #e.g. /MMCI/MS/EpiregDeep/work/TFtoMotifs/hg38.fa 
pathToRepo=$5 #e.g. /MMCI/MS/EpiregDeep2/work/EpiRegioDatabase/EpiRegioDBNina/ApplicationScenarioExamples/

time bash ${pathToRepo}/identifyEnrichedTFs/workflow.sh ${pathToRepo}/identifyEnrichedTFs/JASPAR2020_HUMAN_transfac.txt ${pathToRepo} ${genome}  ${outputDir}/resultMean.txt ${outputDirPASTAA} 0.05

#extract significant TFs (0.01) original pastaa output and run Fimo on their webserver # -> result stored in fimo.tsv
#TODO install fimo proberly

#parse fimo output
python3 ${pathToRepo}/src/parseOutputFIMO.py ${outputDir}/fimo.tsv ${EpiRegioOutput} ${outputDirPASTAA}/PASTAA_output.txt ${outputDir}/info.txt ${outputDir}/heatmapInput.txt

#determine heatmap
Rscript ${pathToRepo}/PASTAA_Fimo_analysis/heatmap.R ${outputDir}/heatmapInput.txt
