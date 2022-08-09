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
    } else {
        gz <- gzfile(file.path(outDir, "tximport.pseudocounts.tab.gz"), "w")
        write.table(cts, row.names=T, col.names=T, quote=F, sep='\t',
            file=gz)
        close(gz)

        gz <- gzfile(file.path(outDir, "tximport.lengths.tab.gz"), "w")
        write.table(txi$length, row.names=T, col.names=T, quote=F, sep='\t',
            file=gz)
        close(gz)
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

