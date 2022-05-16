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
                matching levels in the 'Type' column of the sample Table."),
    make_option(c("-t", "--tx2gene"),        action="store",      metavar="FILE",
    		default="/home/blencowe/blencowe1/ulrich/genomes/mm10/GCvM21_txToGene.txt",
                help="A tab-separated file with, at least, columns transcID, geneID and geneName,
                where transcID and geneID matche the IDs used for Salmon quantification."),
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
                help="Colors per sample type for plotting [%default]")
)

opt <- parse_args(OptionParser(
    usage = "Perform edgeR differential expression analysis based on Salmon gene-level quantification",
    option_list=opt.list),
    args=cArgs)

libMissing <- !require("tximport", quietly=T) && stop("Failed to load R package 'tximport'")
libMissing <- !require("edgeR", quietly=T)    && stop("Failed to load R package 'edgeR'")
libMissing <- !require("gplots", quietly=T)   && stop("Failed to load R package 'gplots'")


### Function defs

saveLog <- function(logName, opt, sampleTab, contrasts, mode) {
### Save a log file
    modeMsg <- switch(mode,
                      salmon = "*** Differential expression analysis of Salmon counts using edgeR ***\n\n== Options ==",
                      vast   = "*** Differential expression analysis of vast-tools counts using edgeR ***\n\n== Options =="
                      )
    selOpt <- which(names(opt) %in% c("sampleTab","tx2gene","rawCounts","geneInfo",
                                      "minCPM","minLFC","maxFDR","minRPKM","outDir","noGO")
                    )
    optTable <- data.frame(option = names(opt)[selOpt],
                           value  = unlist(opt[selOpt]),
                           stringsAsFactors=F
                           )
    optTable$option <- paste0(optTable$option,
                             sapply(max(nchar(optTable$option)) + 1 - nchar(optTable$option), FUN=function(x) {
                                 paste(rep(" ", x), collapse="")
                             })
                             )
    contrTable <- data.frame(Treatment = sapply(contrasts, "[[", 1),
                             Control   = sapply(contrasts, "[[", 2)
                             )
    suppressWarnings({
        write.table(modeMsg,
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName)
        write.table(paste(optTable$option, optTable$value, sep=": "),
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
        write.table("\n== Samples ==",
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
        write.table(sampleTab[,names(sampleTab) %in% c("Sample","Type","File")],
                    row.names=F, col.names=T, quote=F, sep='\t',
                    file=logName, append=T)
        write.table("\n== Contrasts ==",
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
        write.table(contrTable,
                    row.names=F, col.names=T, quote=F, sep='\t',
                    file=logName, append=T)
        write.table("\n== Time stamp ==",
                    row.names=F, col.names=F, quote=F, sep='\t',
                    file=logName, append=T)
    })
}

mergeToAll <- function(x, allgenes=cts) {
### Merge a reduced table x to a larger one based on the ID in allgenes
    tmp <- merge(data.frame(ID = rownames(allgenes), geneInd = 1:nrow(allgenes)),
                 x,
                 by.x=1, by.y=0, all.x=T)
    tmp <- tmp[order(tmp$geneInd), -2]
    names(tmp)[2:3] <- c("log2FC","log2CPM")
    tmp
}

de.table <- function(de, contrasts) {
### Generate a table with DE for all contrasts
    out <- data.frame(geneID = de[[1]]$ID)
    for (i in 1:length(contrasts)) {
        tmp <- signif(de[[i]][,c(2,4,5)],4)
        names(tmp) <- paste(contrasts[[i]][1], ":", contrasts[[i]][2], ".", names(tmp), sep="")
        out <- cbind(out, tmp)
    }
    out
}
    
de.contr <- function(x, opt, de, geneInfo=NULL) {
### Generate a sorted table of differentially expressed genes for a single contrast
    if (min(grep("\\.RP", names(de))) == 2) {
        de <- data.frame(geneName = de$geneName,
                         geneID   = NA,
                         de[,-1])
        
    }
    out <- de[,c(1,2,
                 grep(paste(x[1], "RP", sep="."), names(de)),
                 grep(paste(x[2], "RP", sep="."), names(de)),
                 0:2 + grep(paste(paste(x, collapse="."), "log2FC", sep="."), names(de))
                 )]
    if (!is.null(geneInfo)) {
        out <- data.frame(out, geneInfo[,!(names(geneInfo) %in% names(out))], stringsAsFactors=F)
    }
    out <- out[which((out[,3] >= opt$minRPKM | out[,4] >= opt$minRPKM) &
                     abs(out[,5]) >= opt$minLFC &
                     out[,7]       < opt$maxFDR)
                    ,]
    out <- out[order(out[,5], decreasing=T),]
    if (!is.null(geneInfo)) {
        out <- out[c(1,2,8:ncol(out),3:7)]
    }
    out
}


scatterVsOther <- function(x,y, opt, dat=res$out.de) {
### Scatter plot of log2FC of two contrasts
    xname <- paste(x, collapse=":")
    yname <- paste(y, collapse=":")
    xi  <- which(colnames(dat) == paste(xname, "log2FC", sep="."))
    yi  <- which(colnames(dat) == paste(yname, "log2FC", sep="."))
    xei <- sapply(x, FUN=function(z) {grep(paste(z, "RP", sep="."), colnames(dat))})
    yei <- sapply(y, FUN=function(z) {grep(paste(z, "RP", sep="."), colnames(dat))})
    xexpr <- apply(dat[,xei], MAR=1, max, na.rm=T) > opt$minRPKM
    yexpr <- apply(dat[,yei], MAR=1, max, na.rm=T) > opt$minRPKM
    use   <- xexpr | yexpr
    change1 <- abs(dat[,xi]) > opt$minLFC & dat[,xi + 2] < opt$maxFDR & xexpr
    change2 <- abs(dat[,yi]) > opt$minLFC & dat[,yi + 2] < opt$maxFDR & yexpr
    change <- which(change1 | change2)

    par(mar=c(5,4,3,1))
    lims <- range(as.numeric(c(dat[use, xi], dat[use, yi])), na.rm=T)
    sigcols <- rep(NA, length(use))
    sigcols[change1] <- "dodgerblue"
    sigcols[change2] <- "brown1"
    sigcols[change1 & change2 & sign(dat[,xi]) == sign(dat[,yi])] <- "darkmagenta"
    sigcols[change1 & change2 & sign(dat[,xi]) != sign(dat[,yi])] <- "darkolivegreen2"
    
    plot(dat[use,xi], dat[use,yi], pch=".", xlim=lims, ylim=lims,
         col=ifelse(is.na(sigcols[use]), "black", NA),
         xlab=sub("\\.log2FC", " (log2FC)", colnames(dat)[xi]),
         ylab=sub("\\.log2FC", " (log2FC)", colnames(dat)[yi]))
    title(main=paste("Expression differences"))
    points(dat[,xi], dat[,yi], pch=20, , cex=0.6, col=sigcols)
    abline(a=0, b=1, lty=2)
    
    legend("topleft",
           legend=c("Significant changes:",
                    paste(c(xname, yname,"both","discordant"), " (",
                          c(length(which(sigcols == "dodgerblue")),  length(which(sigcols == "brown1")),
                            length(which(sigcols == "darkmagenta")), length(which(sigcols == "darkolivegreen2"))),
                        ")", sep="")
                    ),
           pch=20, col=c(NA,"dodgerblue","brown1","darkmagenta","darkolivegreen2"), bty="n")
    
    legend("bottomright",
           legend=paste("r =", round(cor(dat[use,xi], dat[use,yi], use="c"), 2)),
           bty="n")
}

saveGOfiles <- function(contr, dat, opt, outDir) {
### Save tables of up- and downregulated genes and background for GO analysis
    i.expr1 <- which(colnames(dat) == paste(contr[1], "RPKM", sep="."))
    i.expr2 <- which(colnames(dat) == paste(contr[2], "RPKM", sep="."))
    i.de    <- which(colnames(dat) == paste0(contr[1], ":", contr[2], ".log2FC"))

    if (grepl("^ENS", dat$geneID[1])) {dat$geneID <- sub("\\.[0-9]+", "", dat$geneID)}
    expr <- dat[,i.expr1] > opt$minRPKM | dat[,i.expr2] > opt$minRPKM
    up   <- expr & dat[,i.de] >  opt$minLFC & dat[,i.de + 2] < opt$maxFDR
    dn   <- expr & dat[,i.de] < -opt$minLFC & dat[,i.de + 2] < opt$maxFDR

    expr <- dat$geneID[which(expr)]
    up   <- dat$geneID[which(up)][order(dat[,i.de][which(up)], decreasing=T)]
    dn   <- dat$geneID[which(dn)][order(dat[,i.de][which(dn)], decreasing=F)]

    contrName <- paste(contr[1], contr[2], sep=":")
    write.table(expr, row.names=F, col.names=F, quote=F, sep='\t',
            file=file.path(outDir, "GO", "input", paste0(contrName, "_BG.txt")))
    write.table(up,   row.names=F, col.names=F, quote=F, sep='\t',
            file=file.path(outDir, "GO", "input", paste0(contrName, "_UP.txt")))
    write.table(dn,   row.names=F, col.names=F, quote=F, sep='\t',
            file=file.path(outDir, "GO", "input", paste0(contrName, "_DN.txt")))   
}

plot.MA <- function(outDir, contr.i, et.i, de.i, opt, rpkmMax) {
### MA-plots highlighting hits with FDR, LFC and RPKM thresholds
    png(file.path(outDir, paste0("DE_", contr.i[1], ":", contr.i[2], ".png")),
        wid=800, hei=800, point=20)
    xy <- plotSmear(et.i, de.tags=de.i$ID[which(de.i$FDR < opt$maxFDR)], deCol="dodgerblue2",
                    main=paste0("Differential expression ", contr.i[1], "/", contr.i[2]))
    selUp <- (rpkmMax > opt$minRPKM & de.i$log2FC >  opt$minLFC & de.i$FDR < opt$maxFDR)[which(!is.na(de.i$FDR))]
    selDn <- (rpkmMax > opt$minRPKM & de.i$log2FC < -opt$minLFC & de.i$FDR < opt$maxFDR)[which(!is.na(de.i$FDR))]
    points(xy$A[selUp], xy$M[selUp], col="red", pch=".", cex=4)
    points(xy$A[selDn], xy$M[selDn], col="red", pch=".", cex=4)
    abline(h=c(-opt$minLFC, opt$minLFC), col="black", lty=3)

    legend("topright",
           legend=c(length(which(selUp)),
                    length(which(de.i$FDR < opt$maxFDR & de.i$log2FC > 0))),
           text.col=c("red","dodgerblue2"), bty="n")
    legend("bottomright",
           legend=c(length(which(de.i$FDR < opt$maxFDR & de.i$log2FC < 0)),
                    length(which(selDn))),
           text.col=c("dodgerblue2","red"), bty="n")
    legend("topleft", paste(length(which(!is.na(de.i$FDR))), "total"), bty="n")
    legend("bottomleft", c(paste0("FDR < ", opt$maxFDR),
                           paste0("FDR < ", opt$maxFDR, ", |LFC| > ",opt$minLFC, ", maxRPKM > ", opt$minRPKM)),
           text.col=c("dodgerblue2","red"), bty="n")

    dev.off()
}

edgeRsalmon <- function(tx2gene, files, opt, outDir, ctl, uqContr) {
    txi <- tximport(files, type="salmon", dropInfReps=TRUE, tx2gene=tx2gene)
    cts <- txi$counts
    
    idMatch <- rownames(cts) %in% tx2gene$geneID
    if (!all(idMatch)) {
        stop("Only ", length(which(idMatch)), "/", length(idMatch),
             " IDs of Salmon output found in tx2gene table")
    }
    
    normMat <- txi$length
    normMat <- normMat / exp(rowMeans(log(normMat)))
    
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    
    dge <- DGEList(cts,
                   group=relevel(as.factor(sampleTab$Type), ref=ctl[1])
                   )
    dge$offset <- t(t(log(normMat)) + o)

    
    dge <- calcNormFactors(dge)
    dge.all <- dge
    use <- rowSums(cpm(dge) > opt$minCPM) >= min(table(sampleTab$Type[sampleTab$Type %in% uqContr]))
    dge$counts <- dge$counts[use,]
    dge$offset <- dge$offset[use,]
    
    dge <- estimateDisp(dge, robust=T)
    
    png(file.path(outDir, "BCVplot.png"), width=480, height=480)
    plotBCV(dge, main="edgeR dispersion estimation")
    garb <- dev.off()
    
    et <- lapply(contrasts, FUN=function(x) {exactTest(dge, pair=rev(x))})
    
    de <- lapply(et, topTags, n=length(which(use)))
    de <- lapply(de, mergeToAll, allgenes=cts)

    elen <- exp(rowMeans(log(txi$length)))    
    rpkm <- rpkmByGroup(dge.all, gene.length=elen)
    rpkm <- rpkm[,sapply(uqTypes, FUN=function(x) {which(colnames(rpkm) == x)})]
    colnames(rpkm) <- paste0(colnames(rpkm), ".RPKM")
    
    genes <- tx2gene[!duplicated(tx2gene$geneID),]
    genes <- merge(data.frame(geneID = rownames(rpkm), ind=1:nrow(rpkm)),
                   genes,
                   by.x=1, by.y="geneID", all.x=T)
    genes <- genes[order(genes$ind), names(genes) %in% c("geneID", "geneName")]
    if (!all(genes$geneID == rownames(rpkm))) {stop("Having trouble matching gene IDs")}
    
    rpkm.sing <- rpkm(dge.all, gene.length=elen)
    colnames(rpkm.sing) <- paste0(colnames(rpkm.sing), ".RPKM")
    out.sing <- data.frame(genes, round(rpkm.sing, 2))
    
    de.all <- de.table(de, contrasts)
    out.de <- data.frame(genes, round(rpkm, 2), de.all[,-1], check.names=F)

    list(dge      = dge,
         et       = et,
         de       = de,
         rpkm     = rpkm,
         de.all   = de.all,
         out.de   = out.de,
         out.sing = out.sing
         )
}

edgeRvast <- function(cts, counts, opt, sampleTab, outDir, ctl, uqContr) {
    dge <- DGEList(cts[rowSums(cts, na.rm=T) > 0,],
                   group=relevel(as.factor(sampleTab$Type), ref=ctl[1])
                   )
    
    dge <- calcNormFactors(dge)
    dge.all <- dge
    use <- rowSums(cpm(dge) > opt$minCPM) >= min(table(sampleTab$Type[sampleTab$Type %in% uqContr]))
    dge$counts <- dge$counts[use,]
    dge$offset <- dge$offset[use,]
    
    dge <- estimateDisp(dge, robust=T)
    
    png(file.path(outDir, "BCVplot.png"), width=480, height=480)
    plotBCV(dge, main="edgeR dispersion estimation")
    garb <- dev.off()
    
    et <- lapply(contrasts, FUN=function(x) {exactTest(dge, pair=rev(x))})
    
    de <- lapply(et, topTags, n=length(which(use)))
    de <- lapply(de, mergeToAll, allgenes=cts)

    if (length(counts) == 1 && is.na(counts)) {        
        rpkm <- cpmByGroup(dge)
        rpkm <- rpkm[,sapply(uqContr, FUN=function(x) {which(colnames(rpkm) == x)})]
        colnames(rpkm) <- paste0(colnames(rpkm), ".RPM")        

        rpkm.sing <- cpm(dge)
        colnames(rpkm.sing) <- paste0(colnames(rpkm.sing), ".RPM")

        rpkm <- rpkm[match(de[[1]]$ID, rownames(rpkm)),]
        rpkm[is.na(rpkm)] <- 0
        rpkm.sing <- rpkm.sing[match(de[[1]]$ID, rownames(rpkm.sing)),]
        rpkm.sing[is.na(rpkm.sing)] <- 0

        genes <- data.frame(geneName = rownames(cts))        
        out.sing  <- data.frame(genes, round(rpkm.sing, 2))
        
    } else {
        rpkm <- as.matrix(
            sapply(unique(sampleTab$Type), FUN=function(x) {
                rowMeans(as.matrix(counts[,names(counts) %in%
                                           paste0(sampleTab$Sample[sampleTab$Type == x], ".cRPKM")]), na.rm=T)
            })
        )
        colnames(rpkm) <- paste0(colnames(rpkm), ".RPKM")
        
        genes <- data.frame(geneID   = counts$ID,
                            geneName = counts$NAME
                            )
        
        rpkm.sing <- counts[,seq(3, ncol(counts) - 1, 2)]
        out.sing  <- data.frame(genes, round(rpkm.sing, 2))
    }
    
    de.all <- de.table(de, contrasts)
    out.de <- data.frame(genes, round(rpkm, 2), de.all[,-1], check.names=F)

    list(dge      = dge,
         et       = et,
         de       = de,
         rpkm     = rpkm,
         de.all   = de.all,
         out.de   = out.de,
         out.sing = out.sing
         )
}

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




