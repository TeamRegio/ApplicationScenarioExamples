#!/bin/bash

EpiRegioOutput=$1 #e.g. PCAT19_EpiRegio_REMs.txt
outputDir=$2 #e.g. results/
outputDirPASTAA=$3 #e.g. PASTAA_result
genome=$4 #e.g. /MMCI/MS/EpiregDeep/work/TFtoMotifs/hg38.fa 
pathToRepo=$5 #e.g. /MMCI/MS/EpiregDeep2/work/EpiRegioDatabase/EpiRegioDBNina/ApplicationScenarioExamples/
fdrCuttoff=0.05

#set path to meme suite if necessary
#export PATH=/opt/meme/bin:/opt/meme/libexec/meme-5.2.0:$PATH

mkdir -p ${outputDir} #create output folder
echo "start PASTAA workflow"
time bash ${pathToRepo}/identifyEnrichedTFs/workflow.sh ${pathToRepo}/identifyEnrichedTFs/JASPAR2020_HUMAN_transfac.txt ${pathToRepo} ${genome}  ${EpiRegioOutput} ${outputDirPASTAA} ${fdrCutoff}



#run Fimo (supress output)
echo "run Fimo"
fimo  --oc ${outputDir}/fimo/ --max-stored-scores 100000000 ${pathToRepo}/identifyTFBindingSites/JASPAR2020_HUMAN_meme.txt ${outputDirPASTAA}/REMs.fa  &>/dev/null

#parse fimo output
echo "parse output"
python3 ${pathToRepo}/src/parseOutputFIMO.py ${outputDir}/fimo/fimo.tsv ${EpiRegioOutput} ${outputDirPASTAA}/PASTAA_output.txt ${outputDir}/info.txt ${outputDir}/heatmapInput.txt

#check if any TFs are enriched for a pvalue of 0.01 (default in the parseOutputFimo file)
lines=$(< ${outputDir}/heatmapInput.txt wc -l)
echo "${lines}"

if [ $lines -gt 3 ] 
then
	#determine heatmap
	echo "determine heatmap"
	Rscript ${pathToRepo}/PASTAA_Fimo_analysis/heatmap.R ${outputDir}/heatmapInput.txt ${outputDir}
	echo "Rscript ${pathToRepo}/PASTAA_Fimo_analysis/heatmap.R ${outputDir}/heatmapInput.txt ${outputDir}"
	
else
	echo "There is no TF with a pvalue <= 0.01! -> no heatmap is generated"
fi

