#!/usr/bin/env Rscript
### Score differential expression in RNA-seq data based on one of three quantification methods
### (salmon, vast-tools, or bam files and samtools idxstats)### using edgeR.
### Outputs a table of gene RPKM in each sample, a table of DE results, and several plots.

cArgs <- commandArgs(TRUE)

## Check dependencies
libMissing <- !require("optparse", quietly=T) && stop("Failed to load R package 'optparse'")


### Capture input
opt.list <- list(
    make_option(c("-s", "--sampleTab"),      action="store", type="character",      metavar="CSVFILE",
                help="Sample table in CSV format with, at least, columns Sample (an arbitrary, unique name),
                and Type (treatment etc.). In case of Salmon counts, also File (locations of the Salmon
                gene quant.sf)."),
    make_option(c("-c", "--contrTab"),       action="store", type="character",      metavar="CSVFILE",
                help="Contrast table in CSV format with, at least, columns 'Experimental' and 'Control',
                matching levels in the 'Type' column of the sample table."),
    make_option(c("-t", "--tx2gene"),        action="store",      metavar="FILE",
    		default="/home/blencowe/blencowe1/ulrich/genomes/mm10/GCvM21_txToGene.txt",
                help="A tab-separated file with, at least, columns transcID, geneID and geneName,
                where transcID and geneID match the IDs used for Salmon quantification."),
    make_option(c("-r", "--rawCounts"),      action="store", type="character",      metavar="FILE",
                help="vast-tools cRPKM_AND_COUNTS table. Mutually exclusive with Salmon output."),
    make_option(c("-i", "--geneInfo"),       action="store", type="character", metavar="FILE",
                help="A tab-separated file with, at least, columns chrom, start, end, strand, geneID,
                geneName, biotype [%default]"),
    make_option(c("-o", "--outDir"),         action="store", default="./DEresults", metavar="DIR",
                help="Output directory. Will be created if it does not exist [%default]"),
    make_option("--noGO",                   action="store_true", default="FALSE",
                help="Do not save files for GO analysis."),
    make_option(c("-e", "--minCPM"),         action="store", default=0.1,           metavar="NUM",
                help="Minimum CPM for DE analysis. Only lines where at least as many samples as the smallest
                group of types at least this CPM will be analyzed [%default]"),
    make_option(c("-l", "--minLFC"),         action="store", default=1,             metavar="NUM",
                help="Minimum log2(fold change) to be considered differential for plotting [%default]"),
    make_option(c("-f", "--maxFDR"),         action="store", default=0.05,          metavar="NUM",
                help="Maximum FDR to be considered differential for plotting [%default]"),
    make_option(c("-x", "--minRPKM"),        action="store", default=5,             metavar="NUM",
                help="Minimum expression RPKM for inclusion in output table [%default]"),
    make_option(c("-C", "--colors"),         action="store",
                default="dodgerblue,grey50,darkmagenta,darkorange1,darkolivegreen2,darkslategray4,burlywood3,cyan3,darkgoldenrod3,firebrick3,navy,seagreen4",
                metavar="NUM",
                help="Colors per sample type for plotting [%default]"),
    make_option("--RlibPath",                action="store", default=NA,             metavar="FILE",
                help="Path to R function library 'include.R' [./R/include.R relative to executable]")
		
)

opt <- parse_args(OptionParser(
    usage = "Perform edgeR differential expression analysis based on Salmon gene-level quantification
             or vast-tools counts table",
    option_list=opt.list),
    args=cArgs)

libMissing <- !require("tximport", quietly=T) && stop("Failed to load R package 'tximport'")
libMissing <- !require("edgeR", quietly=T)    && stop("Failed to load R package 'edgeR'")
libMissing <- !require("gplots", quietly=T)   && stop("Failed to load R package 'gplots'")


### Check and load R library

if (is.na(opt$RlibPath)) {
    argv <- commandArgs(trailingOnly = F)
    scriptPath <- dirname(sub("--file=","",argv[grep("--file",argv)]))
    RlibPath <- file.path(scriptPath, "R", "include.R")
} else {
    RlibPath <- opt$RlibPath
}
if (!file.exists(RlibPath)) {stop("R library not found at ", RlibPath)}
source(RlibPath)


### Body

## Check options
if (is.null(opt$sampleTab))            {stop("--sampleTab must be provided")}
if (is.null(opt$contrTab))             {stop("--contrTab must be provided")}
if (!file.exists(opt$sampleTab))       {stop(opt$sampleTab, " not found")}
if (!file.exists(opt$contrTab))        {stop(opt$contrTab, " not found")}
sampleTab <- read.csv(opt$sampleTab, as.is=T)
contrTab  <- read.csv(opt$contrTab, as.is=T)

if (is.null(opt$rawCounts)) {
    if (is.null(sampleTab$File))       {stop("Provide either paths of Salmon/idxstats output files in the sample Table or vast-tools counts!")}
    if (all(grepl("quant.sf", sampleTab$File))) {
        mode <- "salmon"
        if (is.null(opt$tx2gene))      {stop("--tx2gene not specified")}
        if (!file.exists(opt$tx2gene)) {stop(opt$tx2gene, " not found")}
    } else {
        mode <- "idxstats"
    }
} else {
    if (!is.null(sampleTab$File)) {
        stop("Both Salmon/idxstats output files (in sample table) as well as vast-tools counts provided. Don't know what to do.")
    }
    mode <- "vast"
    if (!file.exists(opt$rawCounts))   {stop(opt$rawCounts, " not found")}
}
if (mode %in% c("vast","idxstats") & !is.null(opt$tx2gene))  {
    cat("--tx2gene provided even though running in", mode, "mode. Ignoring...\n")
    opt$tx2gene <- NULL
}
if (!is.null(opt$geneInfo) && mode == "salmon" && !file.exists(opt$geneInfo)) {
    stop(opt$geneInfo, " not found")
}

pathSlots <- which(names(opt) %in% c("sampleTab","tx2gene","rawCounts","geneInfo","outDir"))
suppressWarnings(opt[pathSlots] <- sapply(opt[pathSlots], normalizePath))

## Check if there are replicates
contrasts <- lapply(1:nrow(contrTab), FUN=function(x) {
    c(contrTab[x, "Experimental"], contrTab[x, "Control"])
})
contrNames <- sapply(contrasts, paste, collapse=":")
ctl <- unique(sapply(contrasts, "[", 2))
if (!all(sapply(contrasts, length) == 2)) {stop("Check format of contrasts string")}
uqTypes <- unique(sampleTab$Type)
uqContr <- unique(unlist(contrasts))
if (!all(uqContr %in% sampleTab$Type)) {
    stop("Sample type(s) ",
         paste(uqContr[!(uqContr %in% sampleTab$Type)], collapse=", "),
         " have no samples associated")
}

reps <- sapply(contrasts, FUN=function(x) {table(sampleTab$Type)[x]})
if (any(colMeans(reps) == 1)) {
    stop("There are no replicates in contrast(s) ",
         paste(sapply(contrasts[colMeans(reps) == 1], paste, collapse="/"), collapse=", "))
}

typeReps <- sapply(uqContr, FUN=function(x) {length(which(sampleTab$Type == x))})
if (any(typeReps == 1)) {
    warning("Sample type(s) ",
            paste(names(typeReps)[typeReps == 1], collapse=", "),
            " has/have only one replicate. This is not recommended.\n")
}

## Load and check input tables
if (mode == "salmon") {
    tx2gene      <- read.delim(opt$tx2gene, header=T, as.is=T)
    files        <- sampleTab$File
    names(files) <- sampleTab$Sample

    if (!all(file.exists(files))) {
        stop("File(s) not found:\n", paste(files[!(file.exists(files))], collapse="\n"))
    } 
    if (!all(c("transcID","geneID","geneName") %in% names(tx2gene))) {
        stop("tx2gene must contain columns 'transcID', 'geneID' and 'geneName'")
    }   
}
if (mode == "idxstats") {
    files        <- sampleTab$File
    names(files) <- sampleTab$Sample

    if (!all(file.exists(files))) {
        stop("File(s) not found:\n", paste(files[!(file.exists(files))], collapse="\n"))
    }
    dat <- lapply(files, read.delim, header=F)
    
    diffFormat <- sapply(dat, FUN=function(x) {
        ncol(x) != 4 || class(x[,1]) != "factor" || class(x[,2]) != "integer" || class(x[,3]) != "integer"
    })
    if (any(diffFormat)) {
        stop("File(s) for sample(s) ", paste(names(diffFormat)[diffFormat], collapse=", "), " have unexpected format")
    }
    
    diffGenes  <- sapply(dat, FUN=function(x) {!all(as.character(x[,1]) == dat[[1]][,1])})
    if (any(diffGenes)) {
        stop("Genes for sample(s) ", paste(names(diffGenes)[diffGenes], collapse=", "), " differ from first sample")
    }

    cts <- sapply(dat, "[[", 3)
    rownames(cts) <- dat[[1]][,1]
}
if (mode == "vast") {
    counts <- read.delim(opt$rawCounts)
    cts <- as.matrix(counts[,seq(4, ncol(counts), 2)])
    rownames(cts) <- counts$ID
    colnames(cts) <- sub("\\.Counts$","",colnames(cts))
    
    sampleFound <- sampleTab$Sample %in% colnames(cts)
    if (!all(sampleFound)) {
        stop("Sample(s) not found in counts table:\n", paste(sampleTab$Sample[!sampleFound], collapse="]n"))
    }
}
if (!is.null(opt$geneInfo) && mode == "salmon") {
    geneInfo  <- read.delim(opt$geneInfo)
    if (is.null(geneInfo$geneID)) {stop("Column geneID not found in --geneInfo")}
} else {
    geneInfo <- NULL
} 

if (!all(uqContr %in% sampleTab$Type)) {
    stop("Sample type(s) ", paste(uqContr[!(uqContr %in% sampleTab$Type)], collapse=", "),
         " not found in sample table")
}

outDir <- sub("/$","", opt$outDir)
if (!dir.exists(outDir)) {dir.create(outDir, recursive=T)}
logName <- file.path(outDir, "DEanalysis.log")


cols <- strsplit(opt$colors, split=",")[[1]]
if (!all(cols %in% colors())) {
    stop("Color(s) not defined: ", paste(cols[!(cols %in% colors())], collapse=", "))
}
if (length(uqTypes) > length(cols)) {cols <- rep(cols, length(uqTypes) %/% length(cols) + 1)}
sampleTab$col <- sapply(sampleTab$Type, FUN=function(x) {cols[which(uqTypes == x)]})


## Save log file
saveLog(logName, opt, sampleTab, contrasts, mode)


## edgeR
if (mode == "salmon")   {res <- edgeRsalmon(tx2gene, files, opt, outDir, ctl, uqContr)}
if (mode == "vast")     {res <- edgeRvast(cts, counts,    opt, sampleTab, outDir, ctl, uqContr)}
if (mode == "idxstats") {res <- edgeRvast(cts, counts=NA, opt, sampleTab, outDir, ctl, uqContr)}

## Output tables
if (!is.null(geneInfo)) {
    geneInfo <- merge(data.frame(geneID=res$out.de$geneID, ind=1:nrow(res$out.de)),
                      geneInfo,
                      by="geneID", all.x=T)
    geneInfo <- geneInfo[order(geneInfo$ind), -2]
    if (length(which(is.na(geneInfo[,2]))) / nrow(res$out.de) > 0.5) {
        warning("< 50% of all genes are represented in --geneInfo, ignoring...")
        geneInfo <- NULL
    }
} 
out.contr <- lapply(contrasts, de.contr, opt=opt, de=res$out.de, geneInfo=geneInfo)

gz <- gzfile(file.path(outDir, "RPKM.samples.tab.gz"), "w")
write.table(res$out.sing, row.names=F, col.names=T, quote=F, sep='\t',
            file=gz)
close(gz)

gz <- gzfile(file.path(outDir, "DE.tab.gz"), "w")
write.table(res$out.de, row.names=F, col.names=T, quote=F, sep='\t',
            file=gz)
close(gz)

for (i in 1:length(contrasts)) {
    write.csv(out.contr[[i]], row.names=F,
                file=file.path(outDir, paste0("DE_", contrasts[[i]][1], ":", contrasts[[i]][2], ".csv")))
}

## MA plots
for (i in 1:length(contrasts)) {
    plot.MA(outDir, contrasts[[i]], res$et[[i]], res$de[[i]], opt,
            rpkmMax = apply(res$rpkm[, sapply(contrasts[[i]], FUN=function(x) {
                grep(paste0(x, ".RP"), colnames(res$rpkm))
            })], MAR=1, max, na.rm=T))
}
                                        


## Plot correlation matrix (of voom-transformed count estimates)
ctsv <- voom(res$dge, normalize.method="quantile")$E
distCts <- dist(t(ctsv))
clust <- as.dendrogram(hclust(distCts))

pdf(file.path(outDir, "SampleCorr_pseudoCts.pdf"),
    wid=4.5 + 0.10 * nrow(sampleTab), hei=4.6 + 0.10 * nrow(sampleTab))
thisCor <- cor(ctsv, use="p")
heatmap.2(thisCor, Rowv=clust, Colv=clust, trace="n",
          mar=rep(10 + 0.11 * length(contrasts), 2),
          col=ifelse(any(as.numeric(thisCor) < 0), "cm.colors", "heat.colors"),
          main="  Sample correlation",
          cexRow=0.2 + 0.8/log10(nrow(sampleTab)), cexCol=0.2 + 0.8/log10(nrow(sampleTab)),
          key.title="", key.xlab="Corr. of log counts")
dev.off()


## Plot MDS
mds <- cmdscale(distCts, k=2)
pdf(file.path(outDir, "MDS_pseudoCts.pdf"), wid=5.5, hei=6)
plot(mds[,c(1,2)], col=sampleTab$col, pch=20,
     main="Multi-dimensional scaling", xlab="DIM 1", ylab="DIM 2")
text(mds[,c(1,2)], labels=rownames(mds), cex=0.8, pos=1, col=sampleTab$col)
dev.off()


## Correlation heatmap of averaged samples
use <- rowSums(res$rpkm > opt$minRPKM) > 0
distMeans <- dist(t(log(0.1 + res$rpkm[use,])))
clustMeans <- as.dendrogram(hclust(distMeans))

pdf(file.path(outDir, "Correlation_RPKM.pdf"))
thisCor <- cor(log(0.1 + res$rpkm[use,]), use="p")
heatmap.2(thisCor, Rowv=clustMeans, Colv=clustMeans, trace="n", mar=c(12,12),
          col=ifelse(any(as.numeric(thisCor) < 0), "cm.colors", "heat.colors"),
          labRow=sub("\\.c*RPKM","",colnames(res$rpkm)),
          labCol=sub("\\.c*RPKM","",colnames(res$rpkm)),
          main=paste0("Correlation of log(RPKM)\n(RPKM > ", opt$minRPKM, ")"),
          key.title="", key.xlab="Correlation coeff.")
dev.off()


## Correlation heatmap of contrasts
if (length(contrasts) > 2) {
    use <- rowSums(res$rpkm > opt$minRPKM) > 0
    distMeans <- dist(t(res$de.all[use, seq(2, ncol(res$de.all) - 2, 3)]))
    clustMeans <- as.dendrogram(hclust(distMeans))
    labels <- sub("\\.log2FC","",colnames(res$de.all[,seq(2, ncol(res$de.all) - 2, 3)]))
    
    pdf(file.path(outDir, "Correlation_log2FC.pdf"))
    thisCor <- cor(res$de.all[use, seq(2, ncol(res$de.all) - 2, 3)], use="p")
    heatmap.2(thisCor,
              col=ifelse(any(as.numeric(thisCor) < 0), "cm.colors", "heat.colors"),
              Rowv=clustMeans, Colv=clustMeans, trace="n", mar=c(14,14),
              labRow=labels, labCol=labels,
              main=paste0("Correlation of log2(FC)\n(RPKM > ", opt$minRPKM, ")"),
              key.title="", key.xlab="Correlation coeff.")
    dev.off()
}

## Scatterplots of contrasts 1-on-1
if (length(contrasts) > 1) {
    pdf(file.path(outDir, "DE_correlation.pairwise.pdf"), wid=5, hei=5.4)
    for (i in 1:length(contrasts)) {
        for (j in 2:length(contrasts)) {
            if (i >= j) next
            scatterVsOther(x=contrasts[[i]], y=contrasts[[j]], opt, res$out.de)
        }
    }
    dev.off()
}

## Save files for GO analysis
if (!opt$noGO) {
    if (!dir.exists(file.path(outDir, "GO", "input"))) {
        dir.create(file.path(outDir, "GO", "input"), recursive=T)
    }
    for (i in 1:length(contrasts)) {
        saveGOfiles(contr=contrasts[[i]], dat=res$out.de, opt=opt, outDir=outDir)
    }
}


## Finish
write.table(paste("Completed", strftime(Sys.time())),
            row.names=F, col.names=F, quote=F, sep='\t',
            file=logName, append=T)




