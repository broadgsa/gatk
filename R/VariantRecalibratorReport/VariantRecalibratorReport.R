library(ellipse);
library(hexbin);

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

plotAnn <- function(d.known, d.novel, d.loci, ann) {
    index = getAnnIndex(d.known, ann);

    k = hist(d.known[,index], breaks=100, plot=FALSE);
    n = hist(d.novel[,index], breaks=100, plot=FALSE);

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

clusterLimits <- function( vals, defaultMin, defaultMax ) {
    x = c(max(defaultMin, min(vals, -2)), min(defaultMax, max(vals, 2)))
    print(x)
    x
}

getClusterColor <- function(clusterIndex, nClusters) {
    clusterColors(nClusters)[clusterIndex]
}

clusterColors <- function(nClusters) {
    rainbow(nClusters)
}


makeAxis <- function( num, vals, off1, mult1, xmin, xmax ) {
    #labels=as.integer(seq(from=min(vals), to=max(vals), by=(abs(min(vals)) + abs(max(vals)))/5))
    #at=seq(from=min((vals - off1)/mult1), to=max((vals - off1)/mult1), by=(abs(min((vals - off1)/mult1)) + abs(max((vals - off1)/mult1)))/5)

    #from = xmin * mult1 + off1
    #to = xmax * mult1 + off1
    #print(list(off1=off1, mult1=mult1, xmin=xmin, xmax=xmax))
    at = as.integer(seq(from=xmin, to=xmax, by=(abs(xmin) + abs(xmax))/5))
    labels = as.integer(at * mult1 + off1)
    #print(list(from=from, to=to, by=(abs(from) + abs(to))/5))
    #print(list(labels=labels, at=at))

    axis(num, labels=labels, at=at);

#    axis(num,
#        labels=as.integer(seq(from=min(vals), to=max(vals), by=(abs(min(vals)) + abs(max(vals)))/5)),
#        at=seq(from=min((vals - off1)/mult1), to=max((vals - off1)/mult1), by=(abs(min((vals - off1)/mult1)) + abs(max((vals - off1)/mult1)))/5)
#    );
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

    xvalsForLims = clusterLimits(d.known[,index1], -4, 4)
    yvalsForLims = clusterLimits(d.known[,index2], -4, 4)
    xlims = c(min(xvalsForLims), 1.2*max(xvalsForLims)); 
    ylims = c(min(yvalsForLims), max(yvalsForLims)); 

    # par(mar=c(5, 6, 2, 5));
    plot(0, 0, type="n", xaxt="n", yaxt="n", xlim=xlims, ylim=ylims, xlab=ann1, ylab=ann2, bty="n");

    mv.known = if (maxVariants == -1 | maxVariants >= nrow(d.known)) { seq(1, nrow(d.known)) } else { as.integer(runif(maxVariants, 1, nrow(d.known)+1))}
    mv.novel = if (maxVariants == -1 | maxVariants >= nrow(d.novel)) { 1:nrow(d.novel) } else { as.integer(runif(maxVariants, 1, nrow(d.novel)+1)) }
    
    print(dim(mv.known))
    print(maxVariants)
    
    points(((d.known[,index1] - off1)/mult1)[mv.known], ((d.known[,index2] - off2)/mult2)[mv.known], pch=19, cex=0.3, col="#0000FF33");
    points(((d.novel[,index1] - off1)/mult1)[mv.novel], ((d.novel[,index2] - off2)/mult2)[mv.novel], pch=19, cex=0.3, col="#FF000033");

    nClusters = length(c)
    for (clusterIndex in c(1:nClusters)) {
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
        color = getClusterColor(clusterIndex, nClusters);
        lineweight = ifelse(weight > 0.50, 4, 3);

        points(mu[1], mu[2], pch=21, col=color, cex=0.5);
        points(ellipse(t(cov), centre=mu), type="l", lwd=lineweight, col=color);
    }

    makeAxis(1, d.novel[,index1], off1, mult1, xvalsForLims[1], xvalsForLims[2]) 
    makeAxis(2, d.novel[,index2], off2, mult2, yvalsForLims[1], yvalsForLims[2]) 

    # add points legend on the lower left
    if (!is.na(d.loci)) {
        legend("bottomleft", c("Known", "Novel", "Suspicious loci"), col=c("blue", "red", "yellow3"), pch=19);
    } else {
        legend("bottomleft", c("Known", "Novel"), col=c("blue", "red"), pch=19);
    }
    
    # add upper right legend with cluster id and weights
    weights = round(sapply(c, function(x) x$mixtureWeight),2)
    clusterNames = paste("C", paste(1:nClusters), sep="")
    clusterLegendNames = paste(clusterNames, weights, sep="-W=") 
    legend("topright", clusterLegendNames, fill=clusterColors(nClusters))

    if (!is.na(d.loci)) {
        points((d.loci[,index1] - off1)/mult1, (d.loci[,index2] - off2)/mult2, pch=19, cex=0.8, col="yellow3");
    }
}

args = commandArgs(TRUE);

plotRoot = args[1];
if (is.na(plotRoot)) { plotRoot = "test"; }

clusterFile = args[2];
if (is.na(clusterFile)) { clusterFile = "/Volumes/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v8/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized"; }

vcfTable = args[3];
if (is.na(vcfTable)) { vcfTable = "/Volumes/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesProcessingPaper/wgs.v8/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.optimized.table"; }

lociFile = args[4];
if (is.na(lociFile) | lociFile == "NA" ) { lociFile = NA; }

maxVariants = args[5];
if (is.na(maxVariants)) { maxVariants = 5000; }
maxVariants = as.integer(maxVariants)

greedy = args[6]
if (is.na(greedy)) { greedy = -1; }
greedy = as.integer(greedy)

l = c();
if (!is.na(lociFile)) {
    t = read.table(lociFile, header=TRUE);
    l = t$POS;
}

print("Greedy reading")
d = read.table(vcfTable, header=TRUE, nrows = greedy);
c = read.clusters(clusterFile);

d.known = d[which(d$DB == 1 | d$ID != "."),];
d.novel = d[which(d$DB == 0 | d$ID == "."),];
d.loci = NA;
if (length(l) > 0) {
    d.loci = d[which(d$POS %in% l),];
}

pdf(paste(plotRoot, ".clusterReport.pdf", sep=""));

for (ann1 in c[[1]]$anns) {
    print(ann1)
    plotAnn(d.known, d.novel, d.loci, ann1);

    for (ann2 in c[[1]]$anns) {
        if (ann1 != ann2) {
	    print(paste("-- v ", ann2))
            plotClusters(d.known, d.novel, d.loci, c, ann1, ann2, maxVariants=maxVariants);
        }
    }
}

dev.off();
