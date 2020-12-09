if (!require("RColorBrewer")) {
 install.packages("RColorBrewer", dependencies = TRUE)
 library(RColorBrewer)
 }

grDevices::dev.set(1)
options(bitmapType='cairo') #important to plot huge heatmaps with ComplexHeatmap package

#install.packages("ComplexHeatmap", repos = "http://cran.us.r-project.org")
library(ComplexHeatmap)
library(circlize)

#input files
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1] #data for the heatmap as matrix TF x REMs
outputDir <-args[2] #outputDir

data <- read.delim(input_file, header = TRUE, sep = "\t", row.names = 1) #REMs x TFs and last column is activity per REM
data = t(data) #for activity

rnames <- rownames(data)
#print(rnames)
data_cellTypeScore = data[, ncol(data)]
print(paste("highest value cellTypeScore", max(data_cellTypeScore)))
print(paste("smallest value cellTypeScore", min(data_cellTypeScore)))

helper = data[,- ncol(data)]
data_activity = helper[, ncol(helper)]
print(paste("highest value activity", max(data_activity)))

helper = helper[,- ncol(helper)]

data_heatmap = helper
print(paste("highest occurrence of a TF", max(data_heatmap)))
#print(data_heatmap)

hclustfunc <- function(x) hclust(x, method="ward.D")
distfunc <- function(x) dist(x,method="manhattan")

cluster_REMs = hclustfunc(distfunc(data_heatmap))
ordering_REMs = rownames(data[cluster_REMs$order,])

cluster_TFs = hclustfunc(distfunc(t(data_heatmap)))

#complex heatmap
#color number TF hits
helper = rev(brewer.pal(9, "BuPu"))
col_fun = colorRamp2(c(0, 1, 2,3,4, 5, 6), c(helper[1], helper[3],helper[5], helper[6], helper[7], helper[8], helper[9])) #TODO: change depending on the highest number of occurence of a TF

#color activity 
helper = rev(brewer.pal(9, "Greys"))
col_fun2 = colorRamp2(c(0.0, 0.5, 1, 2, 2.5, 3, 4), c(helper[1], helper[3], helper[5], helper[6], helper[7], helper[8], helper[9])) #TODO: change depending on the activity

#color cellTypeScore
helper = rev(brewer.pal(11, "RdYlBu"))
col_fun3 = colorRamp2(c(-0.04, -0.03, -0.02, -0.01, 0, 0.01 ,0.02 ,0.03, 0.04, 0.05), c(helper[1], helper[3],helper[4], helper[5], helper[6], helper[7], helper[8], helper[9], helper[10], helper[11])) #TODO: change dependening on the min and max value for the cellTypeScore

pdf(paste(outputDir, "heatmap.pdf", sep = ""))
Heatmap(data_heatmap, name = "?",
	heatmap_legend_param = list(title = "", at = c(0, 1, 2, 3, 4, 5, 6)), #TODO: also change the legend here for the number of TFs
	cluster_rows = as.dendrogram(cluster_REMs), 
	cluster_columns = as.dendrogram(cluster_TFs),
	col = col_fun, 
	column_title = "Transcription factors", 
	column_title_gp = gpar(fontsize = 10),
	column_names_gp = gpar(fontsize = 8), #TODO: maybe change fontsize, dependent on thenumber of TFs shown in the heatmap
	column_title_side = "bottom", #"top"
	row_title = "# TF hits per REM", 
	show_row_names = TRUE,
	row_title_gp = gpar(fontsize = 10),
	row_names_gp = gpar(fontsize = 8), #TODO: maybe change fontsize, dependent on the number of REMs shown in the heatmap
	row_title_side = "right", #"left"
	right_annotation = rowAnnotation(activity = data_activity, CellTypeScore = data_cellTypeScore, col = list(activity = col_fun2, CellTypeScore = col_fun3)),
)
dev.off()
