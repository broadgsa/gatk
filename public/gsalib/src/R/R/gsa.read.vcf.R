gsa.read.vcf <- function(vcffile, skip=0, nrows=-1, expandGenotypeFields = FALSE) {
    headers = readLines(vcffile, n=100);
    headerline = headers[grep("#CHROM", headers)];
    header = unlist(strsplit(gsub("#", "", headerline), "\t"))
    
    d = read.table(vcffile, header=FALSE, skip=skip, nrows=nrows, stringsAsFactors=FALSE);
    colnames(d) = header;

    if (expandGenotypeFields) {
        columns = ncol(d);

        offset = columns + 1;
        for (sampleIndex in 10:columns) {
            gt = unlist(lapply(strsplit(d[,sampleIndex], ":"), function(x) x[1]));
            d[,offset] = gt;
            colnames(d)[offset] = sprintf("%s.GT", colnames(d)[sampleIndex]);

            offset = offset + 1;
        }
    }

    return(d);
}
