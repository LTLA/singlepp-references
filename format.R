library(celldex)
library(SingleR)
library(matrixStats)

formatter <- function(x, prefix, dir) {
    dir.create(dir, showWarnings=FALSE)

    mat <- assay(x, "logcounts")
    rmat <- colRanks(mat, ties.method="min")
    mhandle <- gzfile(file.path(dir, paste0(prefix, "_matrix.csv.gz")))
    write.table(file=mhandle, rmat, sep=",", col.names=FALSE, row.names=FALSE)

    ghandle <- gzfile(file.path(dir, paste0(prefix, "_genes.csv.gz")))
    fdf <- rowData(x)[,c("ensembl", "symbol"),drop=FALSE]
    fdf$ensembl[is.na(fdf$ensembl)] <- ""
    fdf$symbol[is.na(fdf$symbol)] <- ""
    write.table(fdf, file=ghandle, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

    for (lab in colnames(colData(x))) {
        suff <- sub("label\\.", "", lab)
        lhandle <- gzfile(file.path(dir, paste0(prefix, "_labels_", suff, ".csv.gz")))
        write(colData(x)[[lab]], file=lhandle) 

        Y <- getClassicMarkers(x, colData(x)[[lab]], de.n=100)
        curhandle <- gzfile(file.path(dir, paste0(prefix, "_markers_", suff, ".gmt.gz")), open="wb")
        for (i in names(Y)) {
            for (j in names(Y[[i]])) {
                if (i == j) {
                    next
                }
                thing <- c(i, j, match(Y[[i]][[j]], rownames(x)))
                write(paste(thing, collapse="\t"), file=curhandle, append=TRUE)
            }
        }
        close(curhandle)
    }
}

library(AnnotationHub)
ahub <- AnnotationHub()

for (species in c("hs", "mm")) {
    if (species == "hs") {
        ens <- ahub[["AH73881"]]
        datasets <- c("BlueprintEncode", "HumanPrimaryCellAtlas", "MonacoImmune", "DatabaseImmuneCellExpression", "NovershternHematopoietic", "MonacoImmune")
    } else {
        ens <- ahub[["AH73905"]]
        datasets <- c("ImmGen", "MouseRNAseq")
    }

    for (ref in datasets) {
        X <- get(paste0(ref, "Data"))()

        if (!all(grepl("^ENS", rownames(X)))) {
            rowData(X)$ensembl <- mapIds(ens, keys = rownames(X), keytype = "SYMBOL", column = "GENEID")
            rowData(X)$symbol <- rownames(X)
        } else {
            rowData(X)$ensembl <- rownames(X)
            rowData(X)$symbol <- mapIds(ens.hs, keys = rownames(X), keytype = "GENEID", column = "SYMBOL")
        }

        formatter(X, ref, dir=species)
    }
}
