# By Santiago Carmona 2019
# iSEE shinyapp to explore the CD8 TIL dataset in https://www.biorxiv.org/content/10.1101/800847v1


# Environment -------------------------------------------------------------


library(rsconnect)
library(shiny)
library(ggplot2)
#library(plotly)
#library(reshape)
library(SingleCellExperiment)
library(iSEE)

#options(repos = BiocManager::repositories())


# Load Data ---------------------------------------------------------------

x <- readRDS("input/data.sce.rds")


################################################################################
# Settings for reduced dimension plots
################################################################################

redDimPlotArgs <- new('DataFrame', nrows=5L, rownames=paste0('redDimPlot', seq_len(5)))
redDimPlotArgs[['Type']] <- c(2L, 2L, 1L, 1L, 1L)
redDimPlotArgs[['XAxis']] <- c(1, 1, 1, 1, 1)
redDimPlotArgs[['YAxis']] <- c(2, 2, 2, 2, 2)
redDimPlotArgs[['DataBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['VisualBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SelectByPlot']] <- c("Feature assay plot 1", "---", "---", "---", "---")
#redDimPlotArgs[['SelectByPlot']] <- c("Reduced dimension plot 2", "---", "---", "---", "---")
redDimPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
redDimPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
redDimPlotArgs[['SelectColor']] <- c("#FF0000", "#FF0000", "#FF0000", "red", "red")
redDimPlotArgs[['ShapeBy']] <- c("Column data", "None", "None", "None", "None")
redDimPlotArgs[['ShapeByColData']] <- c("cluster.lab", "Sample", "Sample", "Sample", "Sample")


redDimPlotArgs[['VisualChoices']] <- "Color"
redDimPlotArgs[['PointSize']] <- c(4, 1, 1, 1, 1)
redDimPlotArgs[['PointAlpha']] <- c(0.9, 1, 1, 1, 1)
redDimPlotArgs[['Downsample']] <-  c(TRUE, TRUE, TRUE, TRUE, TRUE)
redDimPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
redDimPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
redDimPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

redDimPlotArgs[['ContourAdd']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['ContourColor']] <- c("#0000FF", "#0000FF", "#0000FF", "blue", "blue")
redDimPlotArgs[['ColorBy']] <- c("Feature name", "Column data", "None", "None", "None")

redDimPlotArgs[['ColorByColData']] <- c("cluster.lab", "cluster.lab", "nGene", "nGene", "nGene")
redDimPlotArgs[['ColorByRowTable']] <- c("Row statistics table 1", "---", "---", "---", "---")
redDimPlotArgs[['ColorByFeatName']] <- c(481L, 481L, 1L, 1L, 1L)
redDimPlotArgs[['ColorByFeatNameAssay']] <- c(2L, 2L, 2L, 2L, 2L)
redDimPlotArgs[['ColorByColTable']] <- c("---", "---", "---", "---", "---")
redDimPlotArgs[['ColorBySampName']] <- c(1L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['ColorBySampNameColor']] <- c("#FF0000", "#FF0000", "#FF0000", "red", "red")
redDimPlotArgs[['FacetByRow']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['FacetByColumn']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['RowFacetByColData']] <- c("Sample", "Sample", "Sample", "Sample", "Sample")
redDimPlotArgs[['ColumnFacetByColData']] <- c("Sample", "Sample", "Sample", "Sample", "Sample")

################################################################################
# Settings for feature assay plots
################################################################################

featAssayPlotArgs <- new('DataFrame', nrows=5L, rownames=paste0('featAssayPlot', seq_len(5)))
featAssayPlotArgs[['Assay']] <- c(2L, 2L, 2L, 2L, 2L)
featAssayPlotArgs[['XAxis']] <- c("Column data", "None", "None", "None", "None")
featAssayPlotArgs[['XAxisColData']] <- c("cluster.lab", "nGene", "nGene", "nGene", "nGene")
featAssayPlotArgs[['XAxisFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['XAxisRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['YAxisFeatName']] <- c(481L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['YAxisRowTable']] <- c("Row statistics table 1", "---", "---", "---", "---")
featAssayPlotArgs[['DataBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['VisualBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SelectByPlot']] <- c("Reduced dimension plot 2", "---", "---", "---", "---")
featAssayPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
featAssayPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
featAssayPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

featAssayPlotArgs[['VisualChoices']] <- "Color"
featAssayPlotArgs[['PointSize']] <- c(1, 1, 1, 1, 1)
featAssayPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
featAssayPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
featAssayPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
featAssayPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

featAssayPlotArgs[['ContourAdd']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['ContourColor']] <- c("#0000FF", "blue", "blue", "blue", "blue")
featAssayPlotArgs[['ColorBy']] <- c("Column data", "None", "None", "None", "None")
featAssayPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
featAssayPlotArgs[['ColorByColData']] <- c("cluster.lab", "nGene", "nGene", "nGene", "nGene")
featAssayPlotArgs[['ShapeBy']] <- c("None", "None", "None", "None", "None")
featAssayPlotArgs[['ShapeByColData']] <- c("Sample", "Sample", "Sample", "Sample", "Sample")
featAssayPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['ColorByFeatNameAssay']] <- c(2L, 2L, 2L, 2L, 2L)
featAssayPlotArgs[['ColorByColTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['ColorBySampName']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['ColorBySampNameColor']] <- c("#FF0000", "red", "red", "red", "red")
featAssayPlotArgs[['FacetByRow']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['FacetByColumn']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['RowFacetByColData']] <- c("Sample", "Sample", "Sample", "Sample", "Sample")
featAssayPlotArgs[['ColumnFacetByColData']] <- c("Sample", "Sample", "Sample", "Sample", "Sample")

################################################################################
# Settings for row statistics tables
################################################################################

rowStatTableArgs <- new('DataFrame', nrows=5L, rownames=paste0('rowStatTable', seq_len(5)))
rowStatTableArgs[['Selected']] <- c(481L, 1L, 1L, 1L, 1L)
rowStatTableArgs[['Search']] <- c("Pdcd1", "", "", "", "")

rowStatTableArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
rowStatTableArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")

################################################################################
# Settings for heat maps
################################################################################

heatMapPlotArgs <- new('DataFrame', nrows=5L, rownames=paste0('heatMapPlot', seq_len(5)))
heatMapPlotArgs[['Assay']] <- c(2L, 2L, 2L, 2L, 2L)
heatMapPlotArgs[['FeatNameBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)

heatMapPlotArgs[['FeatName']] <- list(c(13614L, 11993L, 11317L, 756L, 3558L, 6437L, 3691L, 481L, 14518L, 768L, 4422L, 853L, 
                                   11035L, 11250L, 16329L, 1070L, 9871L, 12994L, 2621L, 13104L, 11030L, 13785L, 11777L)
)

heatMapPlotArgs[['ColDataBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)

heatMapPlotArgs[['ColData']] <- "cluster.lab"
heatMapPlotArgs[['FeatNameSource']] <- c("---", "---", "---", "---", "---")

heatMapPlotArgs[['CenterScale']] <- "Scaled"
heatMapPlotArgs[['Lower']] <- c(NA, -Inf, -Inf, -Inf, -Inf)
heatMapPlotArgs[['Upper']] <- c(NA, Inf, Inf, Inf, Inf)
heatMapPlotArgs[['ColorScale']] <- c("blue-white-orange", "purple-black-yellow", "purple-black-yellow", "purple-black-yellow", 
                                     "purple-black-yellow")

heatMapPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
heatMapPlotArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
heatMapPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
heatMapPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
heatMapPlotArgs[['SelectColor']] <- c("red", "red", "red", "red", "red")



################################################################################
# Initial panel settings
################################################################################

initialPanels <- DataFrame(
  Name=c("Row statistics table 1", "Feature assay plot 1", "Reduced dimension plot 1"
  ),
  Width=c(4L, 4L, 8L),
  Height=c(500L, 500L, 800L)
)
#  #, "Reduced dimension plot 2", "Heat map 1"


################################################################################
# Define cluster colors
################################################################################


cluster_colors <- function(n){
  c("#F8766D","#A58AFF","#00B6EB","#53B400")
}
ecm <- ExperimentColorMap(
  global_discrete =  cluster_colors
)



################################################################################
# Run app
################################################################################

introtour <- read.delim("www/intro_firststeps.txt",sep = ";", header = TRUE)

app <- iSEE(x, initialPanels = initialPanels, redDimArgs = redDimPlotArgs, rowStatArgs = rowStatTableArgs, featAssayArgs = featAssayPlotArgs, heatMapArgs = heatMapPlotArgs, colormap = ecm, appTitle = "CD8 TILs from B16 melanoma tumors (dataset Carmona SJ et al)",runLocal = F,tour = introtour)

#runApp(app) # either execute this line or RStudios' "Run App"
