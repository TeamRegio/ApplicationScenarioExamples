library(gplots)
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

data <- read.delim(input_file, header = TRUE, sep = "\t", row.names = 1) #REMs x TFs and last column is activity per REM
data = t(data) #for activity

rnames <- rownames(data)
#print(rnames)
data_cellTypeScore = data[, ncol(data)]
data_cellTypeScore = head(data_cellTypeScore, -2)
#print(data_cellTypeScore)
#print(max(data_cellTypeScore))
#print(min(data_cellTypeScore))
helper = data[,- ncol(data)]
data_activity = helper[, ncol(helper)]
data_activity = head(data_activity, -2)
#print(data_activity)
#print(max(data_activity))
#print(min(data_activity))

helper = helper[,- ncol(helper)]
#extract TSS 
promoter_long_isoform = helper[nrow(helper),]
helper = helper[-nrow(helper),]
promoter_short_isoform = helper[nrow(helper), ]
helper = helper[-nrow(helper),]
dataTSS = cbind(promoter_short_isoform, promoter_long_isoform)
dataTSS= t(dataTSS)

data_heatmap = helper
#print(data_heatmap)
#print(data_heatmap)


hclustfunc <- function(x) hclust(x, method="ward.D")
distfunc <- function(x) dist(x,method="manhattan")

cluster_REMs = hclustfunc(distfunc(data_heatmap))
ordering_REMs = rownames(data[cluster_REMs$order,])
#dd = as.dendrogram(cluster_REMs)
#ordering_REMs = order.dendrogram(dd)

cluster_TFs = hclustfunc(distfunc(t(data_heatmap)))

#complex heatmap
#color activity
helper = rev(brewer.pal(9, "BuPu"))
#col_fun = colorRamp2(c(0, 1, 2,3,5, 10, 60), c(helper[1], helper[3],helper[5], helper[6], helper[7], helper[8], helper[9]))
col_fun = colorRamp2(c(0, 1, 2,3,5, 10, 16), c(helper[1], helper[3],helper[5], helper[6], helper[7], helper[8], helper[9]))

#color activity 
helper = rev(brewer.pal(9, "Greys"))
col_fun2 = colorRamp2(c(0.0, 0.09, 0.5, 1.0, 2.0,3.0,4.5 ), c(helper[1], helper[3], helper[5], helper[6], helper[7], helper[8], helper[9]))
#color cellTypeScore
helper = rev(brewer.pal(11, "RdYlBu"))
col_fun3 = colorRamp2(c(-0.04, -0.03, -0.02, -0.01, 0, 0.01 ,0.02 ,0.03, 0.04, 0.05 ), c(helper[1], helper[3],helper[4], helper[5], helper[6], helper[7], helper[8], helper[9], helper[10], helper[11]))

h1 = Heatmap(data_heatmap, name = "?",
	heatmap_legend_param = list(title = "", at = c(0, 1, 2, 3, 5, 10)), 
	cluster_rows = as.dendrogram(cluster_REMs), 
	cluster_columns = as.dendrogram(cluster_TFs),
	col = col_fun, 
	column_title = "Transcription factors", 
	column_title_gp = gpar(fontsize = 10),
#	column_names_gp = gpar(fontsize = 8),
	column_names_gp = gpar(fontsize = 1),
	column_title_side = "bottom", #"top"
	row_title = "# TF hits per REM", 
	show_row_names = TRUE,
	row_title_gp = gpar(fontsize = 10),
	#row_names_gp = gpar(fontsize = 8),
	row_names_gp = gpar(fontsize = 6),
	row_title_side = "right", #"left"
	right_annotation = rowAnnotation(activity = data_activity, CellTypeScore = data_cellTypeScore, col = list(activity = col_fun2, CellTypeScore = col_fun3)),
)
h2 = Heatmap(dataTSS, name = "2",
	show_heatmap_legend = F,
	col = col_fun, 
	#show_column_names = F,
	show_column_dend = F,
	#column_title = "Transcription factors", 
	#column_title_gp = gpar(fontsize = 10),
	#column_names_gp = gpar(fontsize = 8),
	#column_names_gp = gpar(fontsize = 2),
	#column_title_side = "bottom", #"top"
	row_title = " ", 
	show_row_names = TRUE,
	row_title_gp = gpar(fontsize = 10),
	#row_names_gp = gpar(fontsize = 8),
	row_names_gp = gpar(fontsize = 6),
	row_title_side = "right", #"left"
)
combined = h2 %v%  h1
draw(combined)
