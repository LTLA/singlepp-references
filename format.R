library(celldex)
library(SingleR)
library(matrixStats)

formatter <- function(x, dir) {
    mat <- assay(x, "logcounts")
    rmat <- colRanks(mat, ties.method="first", preserveShape=TRUE)
    dir.create(dir, showWarnings=FALSE)
    mhandle <- gzfile(file.path(dir, "matrix.csv.gz"))
    write.table(file=mhandle, rmat, sep=",", col.names=FALSE, row.names=FALSE)

    ghandle <- gzfile(file.path(dir, "genes.csv.gz"))
    write.table(data.frame(symbol = rownames(x)), file=ghandle, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

    for (lab in colnames(colData(x))) {
        suff <- sub("label\\.", "", lab)
        lhandle <- gzfile(file.path(dir, paste0("labels_", suff, ".csv.gz")))
        write(colData(x)[[lab]], file=lhandle) 

        Y <- getClassicMarkers(x, colData(x)[[lab]], de.n=100)
        curhandle <- gzfile(file.path(dir, paste0("markers_", suff, ".gmt.gz")), open="wb")
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

X <- BlueprintEncodeData(ensembl=TRUE)
formatter(X, "BlueprintEncode")