library(celldex)
library(SingleR)
library(matrixStats)

formatter <- function(x, prefix) {
    dir <- "output"
    dir.create(dir, showWarnings=FALSE)

    mat <- assay(x, "logcounts")
    rmat <- colRanks(mat, ties.method="first", preserveShape=TRUE)
    mhandle <- gzfile(file.path(dir, paste0(prefix, "_matrix.csv.gz")))
    write.table(file=mhandle, rmat, sep=",", col.names=FALSE, row.names=FALSE)

    ghandle <- gzfile(file.path(dir, paste0(prefix, "_genes.csv.gz")))
    write.table(data.frame(symbol = rownames(x)), file=ghandle, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

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

for (ref in c("BlueprintEncode", "ImmGen", "HumanPrimaryCellAtlas", "MonacoImmune", "MouseRNAseq", "DatabaseImmuneCellExpression", "NovershternHematopoietic")) {
    X <- get(paste0(ref, "Data"))(ensembl=TRUE)
    formatter(X, ref)
}
