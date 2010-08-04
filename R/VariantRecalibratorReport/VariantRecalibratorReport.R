library(ellipse);
library(hexbin);
library(rgl);

args = commandArgs(TRUE);

plotRoot = args[1];
if (is.na(plotRoot)) { plotRoot = "test"; }

clusterFile = args[2];
if (is.na(clusterFile)) { clusterFile = "/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v7/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized"; }

vcfTable = args[3];
if (is.na(vcfTable)) { vcfTable = "/home/radon01/kiran/scr1/projects/DataProcessingPaper/scratch/MarkBustedWEx.table"; }

lociFile = args[4];
if (is.na(lociFile)) { lociFile = "/home/radon01/kiran/scr1/projects/DataProcessingPaper/scratch/MarkBustedWEx.loci"; }

maxVariants = args[5];
if (is.na(maxVariants)) { maxVariants = -1; }

getAnnIndex <- function(d, ann) {
    index = -1;
    for (i in c(1:length(names(d)))) {
        if (names(d)[i] == ann) {
            index = i;
        }
    }

    index;
}

getClusterAnnIndex <- function(c, ann) {
    index = -1;

    for (i in c(1:length(c[[1]]$anns))) {
        if (c[[1]]$anns[i] == ann) {
            index = i;
        }
    }

    index;
}

plotAnn <- function(d.known, d.novel, d.loci, ann, plotfile) {
    index = getAnnIndex(d.known, ann);

    k = hist(d.known[,index], breaks=100, plot=FALSE);
    n = hist(d.novel[,index], breaks=100, plot=FALSE);

    pdf(plotfile);

    plot(k$mids, k$density, type="b", col="blue", ylim=c(0, max(k$density)), lwd=2, xlab=ann, ylab="Density", bty="n");
    points(n$mids, n$density, type="b", col="red", lwd=2);

    if (!is.na(d.loci)) {
        legend("topright", c("Known", "Novel", "Suspicious loci"), col=c("blue", "red", "yellow3"), pch=c(21, 21, 18));
    } else {
        legend("topright", c("Known", "Novel"), col=c("blue", "red"), pch=21);
    }

    if (!is.na(d.loci)) {
        for (i in c(1:nrow(d.loci))) {
            points(d.loci[i, index], 0, col="yellow3", pch=18, cex=2.0);
        }
    }

    dev.off();
}

read.clusters <- function(filename) {
    con = file(filename, "r", blocking = FALSE)
    lines = readLines(con)
    close(con);

    anns = c();

    annIndex = 1;
    clusterIndex = 1;
    clusters = c();

    conversions = c();

    for (line in lines) {
        if (length(grep("ANNOTATION", line)) > 0) {
            linePieces = unlist(strsplit(line, ","));

            anns = c(anns, linePieces[2]);
            conversions[[annIndex]] = list(ann = linePieces[2], offset = as.numeric(linePieces[3]), multiplier = as.numeric(linePieces[4]));

            annIndex = annIndex + 1;
        } else if (length(grep("CLUSTER", line)) > 0) {
            linePieces = unlist(strsplit(line, ","));
            
            mixtureWeight = linePieces[2];
            mu = linePieces[3:(3+length(anns)-1)];
            cov = linePieces[(3+length(anns)):length(linePieces)];

            clusters[[clusterIndex]] = list(
                anns = anns,
                conversions = conversions,
                mixtureWeight = as.numeric(mixtureWeight),
                means = as.numeric(mu),
                cov = matrix(cov, nrow=length(anns), ncol=length(anns))
            );
            clusterIndex = clusterIndex + 1;
        }
    }

    clusters;
}

plotClusters <- function(d.known, d.novel, d.loci, c, ann1, ann2, filename, maxVariants = -1) {
    index1 = getAnnIndex(d.known, ann1);
    index2 = getAnnIndex(d.known, ann2);

    cindex1 = getClusterAnnIndex(c, ann1);
    cindex2 = getClusterAnnIndex(c, ann2);

    mult1 = c[[1]]$conversions[[cindex1]]$multiplier;
    off1 = c[[1]]$conversions[[cindex1]]$offset;

    mult2 = c[[1]]$conversions[[cindex2]]$multiplier;
    off2 = c[[1]]$conversions[[cindex2]]$offset;

    xlims = c(min((d.novel[,index1] - off1)/mult1), 1.2*max((d.novel[,index1] - off1)/mult1)); 
    ylims = c(min((d.novel[,index2] - off2)/mult2), max((d.novel[,index2] - off2)/mult2)); 

    clusterColors = c("#A62103", "#F27405", "#F29F05", "#F2B705", "#F2CB05");

    pdf(filename);

    par(mar=c(5, 6, 2, 5));
    plot(0, 0, type="n", xaxt="n", yaxt="n", xlim=xlims, ylim=ylims, xlab=ann1, ylab=ann2, bty="n");

    mv.known = ifelse(maxVariants == -1, nrow(d.known), maxVariants);
    mv.novel = ifelse(maxVariants == -1, nrow(d.novel), maxVariants);
    points(((d.known[,index1] - off1)/mult1)[1:mv.known], ((d.known[,index2] - off2)/mult2)[1:mv.known], pch=19, cex=0.4, col="blue");
    points(((d.novel[,index1] - off1)/mult1)[1:mv.novel], ((d.novel[,index2] - off2)/mult2)[1:mv.novel], pch=19, cex=0.4, col="red");

    for (clusterIndex in c(1:length(c))) {
        mu = c(c[[clusterIndex]]$means[cindex1], c[[clusterIndex]]$means[cindex2]);
        cov = matrix(as.numeric(
                matrix(
                    c(
                        c[[clusterIndex]]$cov[cindex1,cindex1],
                        c[[clusterIndex]]$cov[cindex2,cindex1],
                        c[[clusterIndex]]$cov[cindex1,cindex2],
                        c[[clusterIndex]]$cov[cindex2,cindex2]
                    ),
                    nrow=2, ncol=2
                )
              ), nrow=2, ncol=2
        );

        weight = c[[clusterIndex]]$mixtureWeight;

        if (weight <= 0.20) { color = clusterColors[1]; }
        else if (weight > 0.20 && weight <= 0.40) { color = clusterColors[2]; }
        else if (weight > 0.40 && weight <= 0.60) { color = clusterColors[3]; }
        else if (weight > 0.60 && weight <= 0.80) { color = clusterColors[4]; }
        else if (weight > 0.80) { color = clusterColors[5]; }

        lineweight = ifelse(weight > 0.50, 2, 1);

        points(mu[1], mu[2], pch=21, col=color, cex=0.5);
        points(ellipse(t(cov), centre=mu), type="l", lwd=lineweight, col=color);
    }

    axis(1,
        labels=as.integer(seq(from=min(d.novel[,index1]), to=max(d.novel[,index1]), by=(abs(min(d.novel[,index1])) + abs(max(d.novel[,index1])))/5)),
        at=seq(from=min((d.novel[,index1] - off1)/mult1), to=max((d.novel[,index1] - off1)/mult1), by=(abs(min((d.novel[,index1] - off1)/mult1)) + abs(max((d.novel[,index1] - off1)/mult1)))/5)
    );
    axis(2,
        labels=as.integer(seq(from=min(d.novel[,index2]), to=max(d.novel[,index2]), by=(abs(min(d.novel[,index2])) + abs(max(d.novel[,index2])))/5)),
        at=seq(from=min((d.novel[,index2] - off2)/mult2), to=max((d.novel[,index2] - off2)/mult2), by=(abs(min((d.novel[,index2] - off2)/mult2)) + abs(max((d.novel[,index2] - off2)/mult2)))/5)
    );

    if (!is.na(d.loci)) {
        legend("bottomleft", c("Known", "Novel", "Suspicious loci"), col=c("blue", "red", "yellow3"), pch=19);
    } else {
        legend("bottomleft", c("Known", "Novel"), col=c("blue", "red"), pch=19);
    }

    pieces = 100;
    scale = ((abs(ylims[1]) + abs(ylims[2]))/pieces);
    width = ((abs(xlims[1]) + abs(xlims[2]))/12);
    offset = ylims[1];

    for (i in c(1:pieces)) {
        color = clusterColors[1];
        if (i <= 20) { color = clusterColors[1]; }
        else if (i > 20 && i <= 40) { color = clusterColors[2]; }
        else if (i > 40 && i <= 60) { color = clusterColors[3]; }
        else if (i > 60 && i <= 80) { color = clusterColors[4]; }
        else if (i > 80) { color = clusterColors[5]; }

        polygon(x=xlims[2] + c(0, 0, width, width), y=offset + scale*c((i-1), (i), (i), (i-1)), col=color, border=color);
    }

    axis(4, labels=c(0.0, 0.20, 0.40, 0.60, 0.80, 1.0), at=c(ylims[1], ylims[1] + 20*scale, ylims[1] + 40*scale, ylims[1] + 60*scale, ylims[1] + 80*scale, ylims[2]));
    mtext("Mixture coefficient", 4, line=3);

    if (!is.na(d.loci)) {
        points((d.loci[,index1] - off1)/mult1, (d.loci[,index2] - off2)/mult2, pch=19, cex=0.8, col="yellow3");
    }

    dev.off();
}

d = read.table(vcfTable, header=TRUE);
c = read.clusters(clusterFile);
l = c();
if (!is.na(lociFile)) {
    t = read.table(lociFile, header=TRUE);
    l = t$POS;
}

d.known = d[which(d$DB == 1),];
d.novel = d[which(d$DB == 0),];
d.loci = NA;
if (length(l) > 0) {
    d.loci = d[which(d$POS %in% l),];
}

for (ann1 in c[[1]]$anns) {
    plotAnn(d.known, d.novel, d.loci, ann1, paste(plotRoot, ".anndist.", ann1, ".pdf", sep=""));

    for (ann2 in c[[1]]$anns) {
        if (ann1 != ann2) {
            plotClusters(d.known, d.novel, d.loci, c, ann1, ann2, paste(plotRoot, ".cluster.", ann1, "_vs_", ann2, ".pdf", sep=""), maxVariants=maxVariants);
        }
    }
}
