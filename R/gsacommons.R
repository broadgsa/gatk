plot.begin <- function(plotRoot, name, width=10, height=10) {
    if (!is.na(plotRoot)) {
        filename = paste(plotRoot, ".", name, ".pdf", sep="");
        print(sprintf("Plotting '%s' to '%s'", name, filename));

        pdf(filename, width=width, height=height);
    } else {
        print(sprintf("Plotting '%s' to X11", name));

        x11(width=width, height=height);
    }
}

plot.end <- function(plotRoot) {
    if (!is.na(plotRoot)) {
        dev.off();
    }
}

plot.venn <- function(a, b, c=0, a_and_b, a_and_c=0, b_and_c=0,
                     col=c("#FF6342", "#63C6DE", "#ADDE63"),
                     pos=c(0.20, 0.20, 0.80, 0.82),
                     debug=0
                    ) {
    library(png);
    library(graphics);

    # Set up properties
    for (i in 1:length(col)) {
        rgbcol = col2rgb(col[i]);
        col[i] = sprintf("%02X%02X%02X", rgbcol[1], rgbcol[2], rgbcol[3]);
    }

    chco = paste(col[1], col[2], col[3], sep=",");
    chd = paste(a, b, c, a_and_b, a_and_c, b_and_c, sep=",");

    props = c(
        'cht=v',
        'chs=525x525',
        'chds=0,10000000000',
        paste('chco=', chco, sep=""),
        paste('chd=t:', chd, sep="")
    );
    proplist = paste(props[1], props[2], props[3], props[4], props[5], sep='&');

    # Get the venn diagram (as a temporary file)
    filename = tempfile("venn");
    cmd = paste("wget -O ", filename, " 'http://chart.apis.google.com/chart?", proplist, "' > /dev/null 2>&1", sep="");

    if (debug == 1) {
        print(cmd);
    }
    system(cmd);

    # Render the temp png file into a plotting frame
    a = readPNG(filename);
    
    plot(0, 0, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="");
    rasterImage(a, pos[1], pos[2], pos[3], pos[4]);

    # Clean up!
    unlink(filename);
}

.attemptToLoadFile <- function(filename) {
    file = NA;

    if (file.exists(filename) & file.info(filename)$size > 500) {
        file = read.csv(filename, header=TRUE, comment.char="#");
    }

    file;
}

read.eval <- function(evalRoot) {
    fileAlleleCountStats = paste(evalRoot, ".AlleleCountStats.csv", sep="");
    fileCompOverlap = paste(evalRoot, ".Comp_Overlap.csv", sep="");
    fileCountVariants = paste(evalRoot, ".Count_Variants.csv", sep="");
    fileGenotypeConcordance = paste(evalRoot, ".Genotype_Concordance.csv", sep="");
    fileMetricsByAc = paste(evalRoot, ".MetricsByAc.csv", sep="");
    fileMetricsBySample = paste(evalRoot, ".MetricsBySample.csv", sep="");
    fileQuality_Metrics_by_allele_count = paste(evalRoot, ".Quality_Metrics_by_allele_count.csv", sep="");
    fileQualityScoreHistogram = paste(evalRoot, ".QualityScoreHistogram.csv", sep="");
    fileSampleStatistics = paste(evalRoot, ".Sample_Statistics.csv", sep="");
    fileSampleSummaryStatistics = paste(evalRoot, ".Sample_Summary_Statistics.csv", sep="");
    fileTi_slash_Tv_Variant_Evaluator = paste(evalRoot, ".Ti_slash_Tv_Variant_Evaluator.csv", sep="");
    fileTiTvStats = paste(evalRoot, ".TiTvStats.csv", sep="");
    fileVariant_Quality_Score = paste(evalRoot, ".Variant_Quality_Score.csv", sep="");

    eval = list(
        AlleleCountStats = NA,
        CompOverlap = NA,
        CountVariants = NA,
        GenotypeConcordance = NA,
        MetricsByAc = NA,
        MetricsBySample = NA,
        Quality_Metrics_by_allele_count = NA,
        QualityScoreHistogram = NA,
        SampleStatistics = NA,
        SampleSummaryStatistics = NA,
        TiTv = NA,
        TiTvStats = NA,
        Variant_Quality_Score = NA,

        CallsetNames = c(),
        CallsetOnlyNames = c(),
        CallsetFilteredNames = c()
    );

    eval$AlleleCountStats                = .attemptToLoadFile(fileAlleleCountStats);
    eval$CompOverlap                     = .attemptToLoadFile(fileCompOverlap);
    eval$CountVariants                   = .attemptToLoadFile(fileCountVariants);
    eval$GenotypeConcordance             = .attemptToLoadFile(fileGenotypeConcordance);
    eval$MetricsByAc                     = .attemptToLoadFile(fileMetricsByAc);
    eval$MetricsBySample                 = .attemptToLoadFile(fileMetricsBySample);
    eval$Quality_Metrics_by_allele_count = .attemptToLoadFile(fileQuality_Metrics_by_allele_count);
    eval$QualityScoreHistogram           = .attemptToLoadFile(fileQualityScoreHistogram);
    eval$SampleStatistics                = .attemptToLoadFile(fileSampleStatistics);
    eval$SampleSummaryStatistics         = .attemptToLoadFile(fileSampleSummaryStatistics);
    eval$TiTv                            = .attemptToLoadFile(fileTi_slash_Tv_Variant_Evaluator);
    eval$TiTvStats                       = .attemptToLoadFile(fileTiTvStats);
    eval$Variant_Quality_Score           = .attemptToLoadFile(fileVariant_Quality_Score);

    uniqueJexlExpressions = unique(eval$TiTv$jexl_expression);
    eval$CallsetOnlyNames = as.vector(uniqueJexlExpressions[grep("Filtered|Intersection|none", uniqueJexlExpressions, invert=TRUE, ignore.case=TRUE)]);
    eval$CallsetNames = as.vector(gsub("-only", "", eval$CallsetOnlyNames));
    eval$CallsetFilteredNames = as.vector(c(paste(eval$CallsetNames[1], "-filteredIn", eval$CallsetNames[2], sep=""), paste(eval$CallsetNames[2], "-filteredIn", eval$CallsetNames[1], sep="")));

    eval;
}

eval.getMetrics <- function(eval, jexl_expression) {
    callset.counts = eval$CountVariants[which(eval$CountVariants$evaluation_name == "eval" & eval$CountVariants$comparison_name == "dbsnp" & eval$CountVariants$jexl_expression == jexl_expression),]; 
    callset.counts.titv = eval$TiTv[which(eval$TiTv$evaluation_name == "eval" & eval$TiTv$comparison_name == "dbsnp" & eval$TiTv$jexl_expression == jexl_expression),]; 

    callset.calledCounts = callset.counts[which(callset.counts$filter_name == "called" & callset.counts$novelty_name == "all"),]$nVariantLoci;
    callset.calledCounts.titv = callset.counts.titv[which(callset.counts.titv$filter_name == "called" & callset.counts.titv$novelty_name == "all"),]$ti.tv_ratio;

    callset.knownCounts = callset.counts[which(callset.counts$filter_name == "called" & callset.counts$novelty_name == "known"),]$nVariantLoci;
    callset.knownCounts.titv = callset.counts.titv[which(callset.counts.titv$filter_name == "called" & callset.counts.titv$novelty_name == "known"),]$ti.tv_ratio;

    callset.novelCounts = callset.counts[which(callset.counts$filter_name == "called" & callset.counts$novelty_name == "novel"),]$nVariantLoci;
    callset.novelCounts.titv = callset.counts.titv[which(callset.counts.titv$filter_name == "called" & callset.counts.titv$novelty_name == "novel"),]$ti.tv_ratio;

    callset.filteredCounts = callset.counts[which(callset.counts$filter_name == "filtered" & callset.counts$novelty_name == "all"),]$nVariantLoci;
    callset.filteredCounts.titv = callset.counts.titv[which(callset.counts.titv$filter_name == "filtered" & callset.counts.titv$novelty_name == "all"),]$ti.tv_ratio;
    
    metrics = list(
        all = callset.calledCounts,
        all.titv = callset.calledCounts.titv,

        known = callset.knownCounts,
        known.titv = callset.knownCounts.titv,

        novel = callset.novelCounts,
        novel.titv = callset.novelCounts.titv,

        filtered = callset.filteredCounts,
        filtered.titv = callset.filteredCounts.titv
    );
}

.plot.callsetConcordance.getLabelText <- function(name, metrics, union) {
    text = sprintf("%s (%0.1f%% of union)\nAll: %d, Ti/Tv: %0.1f\nKnown: %d, Ti/Tv: %0.1f\nNovel: %d, Ti/Tv: %0.1f\n",
                   name,             100.0*metrics$all/union$all,
                   metrics$all,      metrics$all.titv,
                   metrics$known,    metrics$known.titv,
                   metrics$novel,    metrics$novel.titv
    );
}

plot.callsetConcordance <- function(eval, col=c("#FF6342", "#63C6DE", "#ADDE63")) {
    aonly = eval.getMetrics(eval, eval$CallsetOnlyNames[1]);
    intersection = eval.getMetrics(eval, "Intersection");
    bonly = eval.getMetrics(eval, eval$CallsetOnlyNames[2]);

    union = list(
        all = intersection$all + aonly$all + bonly$all,
        known = intersection$known + aonly$known + bonly$known,
        novel = intersection$novel + aonly$novel + bonly$novel,
        filtered = intersection$filtered + aonly$filtered + bonly$filtered
    );

    #par.def = par(no.readonly = TRUE);
    #par(mar=c(0.1, 0.1, 0.1, 0.1));

    plot.venn(aonly$all + intersection$all, bonly$all + intersection$all, 0, intersection$all, 0, 0, pos=c(0.32, 0.32, 0.68, 0.70), col=col);

    text(0, 0.5, cex=1.2, pos=4, .plot.callsetConcordance.getLabelText(eval$CallsetNames[1], aonly, union));
    text(0.5, 0.75, cex=1.2, adj=c(0.5, 0.33), .plot.callsetConcordance.getLabelText("Intersection", intersection, union));
    text(1, 0.5, cex=1.2, pos=2, .plot.callsetConcordance.getLabelText(eval$CallsetNames[2], bonly, union));

    #par(par.def);
}

plot.callsetConcordanceByAC <- function(eval, normalize=TRUE, novelty_name="all", col=c("#FF6342", "#5C92A4", "#55BBFF")) {
    aonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[1], novelty_name);
    intersection = eval.getMetricsByAc(eval, "Intersection", novelty_name);
    bonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[2], novelty_name);

    title = paste("Callset concordance per allele count (", novelty_name, " variants)", sep="");

    if (length(intersection$AC) > 0 && length(aonly$AC) == 0) {
        aonly = intersection;
        aonly$n = 0;
    }

    if (length(intersection$AC) > 0 && length(bonly$AC) == 0) {
        bonly = intersection;
        bonly$n = 0;
    }

    #par.def = par(no.readonly = TRUE);
    #par(mar=c(5, 5, 3, 5));

    if (normalize == TRUE) {
        norm = aonly$n + intersection$n + bonly$n;
        matnorm = rbind(aonly$n/norm, intersection$n/norm, bonly$n/norm);

        barplot(matnorm, col=col, xlab="Allele count", ylab="Fraction", main=title, names.arg=intersection$AC, xlim=c(0, max(intersection$AC)), ylim=c(0, 1.3), border=NA);
    } else {
        mat = rbind(aonly$n, intersection$n, bonly$n);

        barplot(mat, col=col, xlab="Allele count", ylab="counts", main=title, names.arg=intersection$AC, xlim=c(0, max(intersection$AC)), ylim=c(0, 1), border=NA);
    }

    legend("topright", c(eval$CallsetOnlyNames[1], "Intersection", eval$CallsetOnlyNames[2]), fill=col);

    #par(par.def);
}

plot.alleleCountSpectrum <- function(eval, novelty_name="all", col=c("#FF6342", "#5C92A4", "#55BBFF")) {
    aonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[1], novelty_name);
    intersection = eval.getMetricsByAc(eval, "Intersection", novelty_name);
    bonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[2], novelty_name);

    title = paste("Allele count spectrum (", novelty_name, " variants)", sep="");

    if (length(intersection$AC) > 0 && length(aonly$AC) == 0) {
        aonly = intersection;
        aonly$n = 0;
    }

    if (length(intersection$AC) > 0 && length(bonly$AC) == 0) {
        bonly = intersection;
        bonly$n = 0;
    }

    suppressWarnings(plot(0, 0, type="n", xlim=c(1, length(intersection$AC)), ylim=c(1, max(aonly$n + intersection$n, bonly$n + intersection$n)), xlab="Allele count", ylab="Number of variants", main=title, log="xy", bty="n"));
    suppressWarnings(points(intersection$AC, aonly$n + intersection$n, type="l", lwd=2, col=col[1]));
    suppressWarnings(points(intersection$AC, intersection$n, type="l", lwd=2, col=col[2]));
    suppressWarnings(points(intersection$AC, bonly$n + intersection$n, type="l", lwd=2, col=col[3]));

    legend("topright", c(eval$CallsetOnlyNames[1], "Intersection", eval$CallsetOnlyNames[2]), fill=col);
}

eval.getMetricsByAc <- function(eval, jexl_expression, novelty_name="all") {
    piece = eval$MetricsByAc[which(eval$MetricsByAc$evaluation_name == "eval" & eval$MetricsByAc$comparison_name == "dbsnp" & as.character(eval$MetricsByAc$jexl_expression) == as.character(jexl_expression) & eval$MetricsByAc$filter_name == "called" & eval$MetricsByAc$novelty_name == novelty_name),]; 
}

plot.titvSpectrum <- function(eval, novelty_name="all", col=c("#FF6342", "#5C92A4", "#55BBFF")) {
    aonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[1], novelty_name);
    intersection = eval.getMetricsByAc(eval, "Intersection", novelty_name);
    bonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[2], novelty_name);

    title = paste("Ti/Tv spectrum (", novelty_name, " variants)", sep="");

    if (length(intersection$AC) > 0 && length(aonly$AC) == 0) {
        aonly = intersection;
        aonly$n = 0;
        aonly$nTi = 0;
        aonly$nTv = 0;
    }

    if (length(intersection$AC) > 0 && length(bonly$AC) == 0) {
        bonly = intersection;
        bonly$n = 0;
        bonly$nTi = 0;
        bonly$nTv = 0;
    }

    titv.aonly = (aonly$nTi + intersection$nTi)/(aonly$nTv + intersection$nTv);
    titv.aonly.finite = titv.aonly[which(is.finite(titv.aonly))];

    titv.intersection.finite = intersection$Ti.Tv[which(is.finite(intersection$Ti.Tv))];

    titv.bonly = (bonly$nTi + intersection$nTi)/(bonly$nTv + intersection$nTv);
    titv.bonly.finite = titv.bonly[which(is.finite(titv.bonly))];

    titv.min = min(titv.aonly.finite, titv.intersection.finite, titv.bonly.finite);
    titv.max = max(titv.aonly.finite, titv.intersection.finite, titv.bonly.finite);

    plot(0, 0, type="n", xlim=c(1, length(intersection$AC)), ylim=c(titv.min, titv.max), xlab="Allele count", ylab="Transition/transversion (Ti/Tv) ratio", main=title, bty="n");
    points(intersection$AC, (aonly$nTi + intersection$nTi)/(aonly$nTv + intersection$nTv), type="l", lwd=2, col=col[1]);
    points(intersection$AC, intersection$Ti.Tv, type="l", lwd=2, col=col[2]);
    points(intersection$AC, (bonly$nTi + intersection$nTi)/(bonly$nTv + intersection$nTv), type="l", lwd=2, col=col[3]);

    abline(h=2.3, lty=2);
    mtext("2.3", side=4, at=2.3, cex=0.9);

    abline(h=3.3, lty=2);
    mtext("3.3", side=4, at=3.3, cex=0.9);

    legend("topleft", c(eval$CallsetOnlyNames[1], "Intersection", eval$CallsetOnlyNames[2]), fill=col);
}

plot.variantsPerSample <- function(eval) {
    if (!is.na(eval$MetricsBySample)) {
        metrics.all = eval$MetricsBySample[which(eval$MetricsBySample$evaluation_name == "eval" & eval$MetricsBySample$comparison_name == "dbsnp" & as.character(eval$MetricsBySample$jexl_expression) == "none" & eval$MetricsBySample$filter_name == "called" & eval$MetricsBySample$novelty_name == "all"),]; 
        metrics.known = eval$MetricsBySample[which(eval$MetricsBySample$evaluation_name == "eval" & eval$MetricsBySample$comparison_name == "dbsnp" & as.character(eval$MetricsBySample$jexl_expression) == "none" & eval$MetricsBySample$filter_name == "called" & eval$MetricsBySample$novelty_name == "known"),]; 
        metrics.novel = eval$MetricsBySample[which(eval$MetricsBySample$evaluation_name == "eval" & eval$MetricsBySample$comparison_name == "dbsnp" & as.character(eval$MetricsBySample$jexl_expression) == "none" & eval$MetricsBySample$filter_name == "called" & eval$MetricsBySample$novelty_name == "novel"),]; 

        title = "Calls per sample";
        indices = order(metrics.all$nVariants, decreasing=TRUE);

        plot(0, 0, type="n", xaxt="n", xlim=c(1, length(metrics.all$sample)), ylim=c(0, max(metrics.all$nVariants)), xlab="", ylab="Number of variants", main=title, bty="n");
        points(c(1:length(metrics.all$sample)), (metrics.all$nVariants)[indices], pch=21, col="black");
        points(c(1:length(metrics.known$sample)), (metrics.known$nVariants)[indices], pch=21, col="blue");
        points(c(1:length(metrics.novel$sample)), (metrics.novel$nVariants)[indices], pch=21, col="red");

        legend("topright", c("All", "Known", "Novel"), pch=21, col=c("black", "blue", "red"));

        axis(1, at=c(1:length(metrics.all$sample)), labels=(metrics.all$sample)[indices], las=2, cex.axis=0.4);
    }
}
