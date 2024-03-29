library(celldex)
library(SingleR)
library(matrixStats)

formatter <- function(x, prefix, dir = "v1.0.0") {
    dir.create(dir, showWarnings=FALSE)

    mat <- assay(x, "logcounts")
    rmat <- colRanks(mat, ties.method="min")
    mhandle <- gzfile(file.path(dir, paste0(prefix, "_matrix.csv.gz")))
    write.table(file=mhandle, rmat, sep=",", col.names=FALSE, row.names=FALSE)

    ghandle <- gzfile(file.path(dir, paste0(prefix, "_genes.csv.gz")))
    fdf <- rowData(x)[,c("ensembl", "symbol", "entrez"),drop=FALSE]
    fdf$ensembl[is.na(fdf$ensembl)] <- ""
    fdf$symbol[is.na(fdf$symbol)] <- ""
    fdf$entrez[is.na(fdf$entrez)] <- ""
    write.table(fdf, file=ghandle, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

    for (lab in colnames(colData(x))) {
        suff <- sub("label\\.", "", lab)
        curlabs <- colData(x)[[lab]]

        ulabs <- sort(unique(curlabs))
        uhandle <- gzfile(file.path(dir, paste0(prefix, "_label_names_", suff, ".csv.gz")))
        write(ulabs, file=uhandle, ncolumns=1)

        lhandle <- gzfile(file.path(dir, paste0(prefix, "_labels_", suff, ".csv.gz")))
        write(match(curlabs, ulabs) - 1L, file=lhandle, ncolumns=1)

        Y <- getClassicMarkers(x, curlabs, de.n=100)
        curhandle <- gzfile(file.path(dir, paste0(prefix, "_markers_", suff, ".gmt.gz")), open="wb")
        for (i in names(Y)) {
            for (j in names(Y[[i]])) {
                if (i == j) {
                    next
                }
                thing <- c(match(c(i, j), ulabs), match(Y[[i]][[j]], rownames(x))) - 1L
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
        ens <- ahub[["AH109336"]] # Ensembl 108.
        datasets <- c("BlueprintEncode", "HumanPrimaryCellAtlas", "MonacoImmune", "DatabaseImmuneCellExpression", "NovershternHematopoietic", "MonacoImmune")
    } else {
        ens <- ahub[["AH109367"]] # Ensembl 108.
        datasets <- c("ImmGen", "MouseRNAseq")
    }

    for (ref in datasets) {
        X <- get(paste0(ref, "Data"))()

        if (!all(grepl("^ENS", rownames(X)))) {
            rowData(X)$ensembl <- mapIds(ens, keys = rownames(X), keytype = "SYMBOL", column = "GENEID")
            rowData(X)$entrez <- mapIds(ens, keys = rownames(X), keytype = "SYMBOL", column = "ENTREZID")
            rowData(X)$symbol <- rownames(X)
        } else {
            rowData(X)$ensembl <- rownames(X)
            rowData(X)$entrez <- mapIds(ens, keys = rownames(X), keytype = "GENID", column = "ENTREZID")
            rowData(X)$symbol <- mapIds(ens.hs, keys = rownames(X), keytype = "GENEID", column = "SYMBOL")
        }

        formatter(X, ref)
    }
}
