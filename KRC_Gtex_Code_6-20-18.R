rm(list = ls())

#### 

##-- List packages required in this analysis
cpan.pkg.list <- c('ggplot2', 'ape', 'scales', 'RColorBrewer', 
                   'reshape','VennDiagram')
bioc.pkg.list <- c('ctc',  'limma', 'edgeR', 'DESeq2', 'vsn', 
                   'genefilter', 'pheatmap', 
                   'clusterProfiler', 'pathview',
                   'AnnotationHub')

##-- Set up CPAN repo (required if running IRkernel in Jupyter)
cpan.repos <- 'http://cran.us.r-project.org'

##-- Set up Bioconductor repo
source("https://bioconductor.org/biocLite.R")

##-- Install Bioc packages
biocLite('ctc')
biocLite('limma')
biocLite('edgeR')
biocLite('DESeq2')
biocLite('vsn')
biocLite('genefilter')
biocLite('pheatmap')
biocLite('clusterProfiler')
biocLite('pathview')
biocLite('AnnotationHub')

#Load all libraries
for(pkg in c(cpan.pkg.list, bioc.pkg.list)) {
  print(pkg)
  suppressMessages(library(pkg, character.only = TRUE))
}

library(readr)

##-- Set Parameters
data.base <- 'GTEx' 
fdr <- 0.05
fc <- 1.5
gene.type <- 'coding'
caller <- 'deseq2'
group1 <- 'Whole Blood'
group2 <- 'Spleen'
group3 <- 'EBV'
colors <- c('red', 'blue', 'darkgreen' )

#Set WD and paths
path <- '/Users/Kyle/Box Sync/Gajewski/bioinformatics/GTEx_PRKCD_correlation'

setwd(path)

## check if the working directory has been successfully set
getwd()

## list files in current directory
list.files(path = '.')


##-- Set up working directory
work.dir <- '.'
setwd(work.dir)

##-- Input/Output directories
in.dir <- '../Raw_data'
code.dir <- '../Code'
out.dir <- '../Output'

expr.file <- paste0(data.base, '_Raw_Gene_reads.Rdata')
sample.file <- paste0(data.base, 'ID_noNA_long.csv')
geneinfo.file <- 'gencode.v24.primary_assembly.annotation.gtf.geneinfo'

print(paste0('Database = ', data.base))
print(paste0('gene type = ', gene.type))
print(paste0('DEG fdr cutoff = ', fdr))
print(paste0('DEG fc cutoff = ', fc))
print(paste0('Expression file = ', expr.file))
print(paste0('Sample group file  = ', sample.file))
print(paste0('Gene info file  = ', geneinfo.file))

##-- Input/Output files

# Import annotation for organs and sample IDs,
# Tissue_SampID <- read_delim("~/Box Sync/Gajewski/Proposed Figs/prkcd_coexpression/Annotation_ID_Tissue.txt", 
#                                "\t", escape_double = FALSE, trim_ws = TRUE)
# View(Tissue_SampID)
#
# Raw_Gene_reads <- read_delim("~/Box Sync/Gajewski/Proposed Figs/prkcd_coexpression/Gene_reads", 
#                                  "\t", escape_double = FALSE, trim_ws = TRUE, 
#                                  skip = 2)
#
# save(Raw_Gene_reads, file = "Raw_Gene_reads.Rdata")
# save(Tissue_SampID, file = "Tissue_SampID.Rdata")

#Load raw gene read counts
load("./Raw_data/Raw_Gene_reads.Rdata")

## using the terminal, create a csv file of sample IDs for each organ. 
## grep organ Annotation_ID_Tissue.txt > IDs_organ
## Did this for wholeblood (grep Blood), Spleen (grep Spleen), and LCLs (grep EBV)
## make sure all matricies are in wide format (samples in a single row)
#Load matrix with correct sample IDs for organ type desired. Make sure the samples are in wide format

sample_IDs <- read_csv("./IDonly_noNA_wide.csv", col_names = FALSE)



#Create new DF where only the columns in organ ID that are matched with raw data are present. Repeat for each organ.
genereads_parsed = Raw_Gene_reads[,colnames(Raw_Gene_reads) %in% c('Name','Description',sample_IDs),drop=F]

cds <- genereads_parsed[,-c(1:2)]
# Now you have a DF with only gene read counts and no gene names


#import a file of all the IDs and organ names for each sample merged into 1 long format
sampleID_long_all <- read_csv("./ID_noNA_long.csv", col_names = TRUE)

# Add +1 or +3 to all values, to remove 0 values from matrix, Later log transformations
# cannot handle values of zero
## test <- all_genereads_noID + 1

##-- DESeq2: build DESeqDataSet object, prepare design matrix
dds <- DESeqDataSetFromMatrix(countData = cds,
                              colData = sampleID_long_all, 
                              design = ~ Organ)

#vsd normalization
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

head(assay(vsd), 1)

#rld normalization
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 1)

# Cannot get this to work, error in row names or something
#write.table(data.frame(Name = row.names(assay(vsd)), assay(vsd)),
#            file = 'vsdnorm.txt',
#            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

pdf(paste0('meanvar_vsd.pdf'),width = 7, height = 7)
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()


options(jupyter.plot_mimetypes = "image/svg+xml") 
options(repr.plot.width = 6, repr.plot.height = 5)

##-- DESeq2: sample correlation heatmap
##-- calculate sample distance 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Group, vsd$Libtype, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

##-- DESeq2: plot heatmap
heatmap.colors <- rev(cm.colors(32))[1:16]
pdf(paste0('sm_cor.pdf'),width = 7, height = 7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=heatmap.colors)
dev.off()

##-- DESeq2: Principal component analysis (PCA) plot of the samples
##-- use ggplot2 to customize the PCA plot
pdf(paste0( 'DESeq2_pca.pdf'),width = 7, height = 10)
data.pca <- plotPCA(vsd, intgroup=c('Organ'), returnData=TRUE)
percent.var <- round(100 * attr(data.pca, "percentVar"))
pca.colors <- c(KO = colors[1], WT = colors[2])
p1 <- ggplot(data.pca, aes(PC1, PC2, color = "Organ")) +
  geom_point(size = 5, shape = 17) +
  scale_color_discrete(h = c(0 , 697)) + 
  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
  ylab(paste0("PC2: ",percent.var[2],"% variance"))
plot(p1)
dev.off()

##-- Set up R plot display options in notebook
options(jupyter.plot_mimetypes = "image/svg+xml") 
options(repr.plot.width = 6, repr.plot.height = 5)

##-- DESeq2: remove genes not expressed in any samples
##-- for plottig purposes
notAllZero <- (rowSums(counts(dds))>0)

##-- DESeq2: mean to var plots 
pdf(paste0('meanvar_log2.pdf'),width = 7, height = 7)
meanSdPlot(log2(counts(estimateSizeFactors(dds),
                       normalized=TRUE)[notAllZero,] + 1))
dev.off()


pdf(paste0('meanvar_vsd.pdf'),width = 7, height = 7)
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()


