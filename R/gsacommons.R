suppressPackageStartupMessages(library(gplots));

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
    if (a >= b) {
        rasterImage(a, pos[1], pos[2], pos[3], pos[4]);
    } else {
        rasterImage(a, 0.37+pos[1], 0.37+pos[2], 0.37+pos[3], 0.37+pos[4], angle=180);
    }

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
    fileSimpleMetricsBySample = paste(evalRoot, ".SimpleMetricsBySample.csv", sep="");
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
        SimpleMetricsBySample = NA,
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
    eval$SimpleMetricsBySample           = .attemptToLoadFile(fileSimpleMetricsBySample);
    eval$TiTv                            = .attemptToLoadFile(fileTi_slash_Tv_Variant_Evaluator);
    eval$TiTvStats                       = .attemptToLoadFile(fileTiTvStats);
    eval$Variant_Quality_Score           = .attemptToLoadFile(fileVariant_Quality_Score);

    uniqueJexlExpressions = unique(eval$TiTv$jexl_expression);
    eval$CallsetOnlyNames = as.vector(uniqueJexlExpressions[grep("FilteredIn|Intersection|none", uniqueJexlExpressions, invert=TRUE, ignore.case=TRUE)]);
    eval$CallsetNames = as.vector(gsub("-only", "", eval$CallsetOnlyNames));
    eval$CallsetFilteredNames = as.vector(c(
        paste(gsub("^(\\w)", "In\\U\\1", eval$CallsetNames[1], perl=TRUE), "-Filtered", gsub("^(\\w)", "In\\U\\1", eval$CallsetNames[2], perl=TRUE), sep=""),
        paste(gsub("^(\\w)", "In\\U\\1", eval$CallsetNames[2], perl=TRUE), "-Filtered", gsub("^(\\w)", "In\\U\\1", eval$CallsetNames[1], perl=TRUE), sep=""))
    );

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

    callset.allFilteredCounts = callset.counts[which(callset.counts$filter_name == "filtered" & callset.counts$novelty_name == "all"),]$nVariantLoci;
    callset.allFilteredCounts.titv = callset.counts.titv[which(callset.counts.titv$filter_name == "filtered" & callset.counts.titv$novelty_name == "all"),]$ti.tv_ratio;

    callset.knownFilteredCounts = callset.counts[which(callset.counts$filter_name == "filtered" & callset.counts$novelty_name == "known"),]$nVariantLoci;
    callset.knownFilteredCounts.titv = callset.counts.titv[which(callset.counts.titv$filter_name == "filtered" & callset.counts.titv$novelty_name == "known"),]$ti.tv_ratio;

    callset.novelFilteredCounts = callset.counts[which(callset.counts$filter_name == "filtered" & callset.counts$novelty_name == "novel"),]$nVariantLoci;
    callset.novelFilteredCounts.titv = callset.counts.titv[which(callset.counts.titv$filter_name == "filtered" & callset.counts.titv$novelty_name == "novel"),]$ti.tv_ratio;
    
    metrics = list(
        all = callset.calledCounts,
        all.titv = callset.calledCounts.titv,

        known = callset.knownCounts,
        known.titv = callset.knownCounts.titv,

        novel = callset.novelCounts,
        novel.titv = callset.novelCounts.titv,

        filtered.all = callset.allFilteredCounts,
        filtered.all.titv = callset.allFilteredCounts.titv,

        filtered.known = callset.knownFilteredCounts,
        filtered.known.titv = callset.knownFilteredCounts.titv,

        filtered.novel = callset.novelFilteredCounts,
        filtered.novel.titv = callset.novelFilteredCounts.titv
    );
}

.plot.callsetConcordance.getLabelText <- function(name, othername, metrics, filtered.metrics=NA, union) {
    if (is.na(filtered.metrics)) {
        text = sprintf("%s (%0.01f%% of union)\nCalled:\nAll: %d, Ti/Tv: %0.2f\nKnown: %d, Ti/Tv: %0.2f\nNovel: %d, Ti/Tv: %0.2f",
                   name,             100*metrics$all/union$all.withfiltered,
                   metrics$all,      metrics$all.titv,
                   metrics$known,    metrics$known.titv,
                   metrics$novel,    metrics$novel.titv
        );
    } else {
        text = sprintf("%s (%0.01f%% of union)\nCalled in %s, filtered in %s:\nAll: %d, Ti/Tv: %0.2f\nKnown: %d, Ti/Tv: %0.2f\nNovel: %d, Ti/Tv: %0.2f\n\nCalled in %s, absent in %s:\nAll: %d, Ti/Tv: %0.2f\nKnown: %d, Ti/Tv: %0.2f\nNovel: %d, Ti/Tv: %0.2f",
                   name,             100*(metrics$all + filtered.metrics$all)/union$all.withfiltered,

                   name,                      othername,
                   filtered.metrics$all,      filtered.metrics$all.titv,
                   filtered.metrics$known,    filtered.metrics$known.titv,
                   filtered.metrics$novel,    filtered.metrics$novel.titv,

                   name,             othername,
                   metrics$all,      metrics$all.titv,
                   metrics$known,    metrics$known.titv,
                   metrics$novel,    metrics$novel.titv
        );
    }
}

plot.titlePage <- function(title, author) {
    textplot(sprintf("Automated Data Processing Report\n\n%s\n%s\n%s\n", title, author, Sys.Date()));
}

.plot.variantTable.getRowText <- function(eval, jexl_expression) {
    allVariants = eval$CountVariants[which(eval$CountVariants$jexl_expression == jexl_expression & eval$CountVariants$filter_name == "called" & eval$CountVariants$novelty_name == "all"),]$nVariantLoci;
    knownVariants = eval$CountVariants[which(eval$CountVariants$jexl_expression == jexl_expression & eval$CountVariants$filter_name == "called" & eval$CountVariants$novelty_name == "known"),]$nVariantLoci;
    novelVariants = eval$CountVariants[which(eval$CountVariants$jexl_expression == jexl_expression & eval$CountVariants$filter_name == "called" & eval$CountVariants$novelty_name == "novel"),]$nVariantLoci;

    allTiTv = eval$TiTv[which(eval$TiTv$jexl_expression == jexl_expression & eval$TiTv$filter_name == "called" & eval$TiTv$novelty_name == "all"),]$ti.tv_ratio;
    knownTiTv = eval$TiTv[which(eval$TiTv$jexl_expression == jexl_expression & eval$TiTv$filter_name == "called" & eval$TiTv$novelty_name == "known"),]$ti.tv_ratio;
    novelTiTv = eval$TiTv[which(eval$TiTv$jexl_expression == jexl_expression & eval$TiTv$filter_name == "called" & eval$TiTv$novelty_name == "novel"),]$ti.tv_ratio;

    cbind(allVariants, knownVariants, sprintf("%0.2f", knownTiTv), novelVariants, sprintf("%0.2f", novelTiTv));
}

plot.variantTable <- function(eval, title) {
    aonly.row = .plot.variantTable.getRowText(eval, eval$CallsetOnlyNames[1]);
    aonly.filtered.row = .plot.variantTable.getRowText(eval, eval$CallsetFilteredNames[1]);
    intersection.row = .plot.variantTable.getRowText(eval, "Intersection");
    bonly.row = .plot.variantTable.getRowText(eval, eval$CallsetOnlyNames[2]);
    bonly.filtered.row = .plot.variantTable.getRowText(eval, eval$CallsetFilteredNames[2]);

    variantsummary = as.data.frame(rbind(bonly.row, bonly.filtered.row, intersection.row, aonly.filtered.row, aonly.row));

    rownames(variantsummary) = c(
        sprintf("Called in %s, absent in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
        sprintf("Called in %s, filtered in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
        "Intersection",
        sprintf("Called in %s, filtered in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2]),
        sprintf("Called in %s, absent in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2])
    );
    colnames(variantsummary) = c("counts (all)", "counts (known)", "ti/tv (known)", "counts (novel)", "ti/tv (novel)");

    textplot(variantsummary);
}

plot.callsetConcordance <- function(eval, col=c("#FF6342", "#63C6DE", "#ADDE63")) {
    aonly = eval.getMetrics(eval, eval$CallsetOnlyNames[1]);
    aonly.filtered = eval.getMetrics(eval, eval$CallsetFilteredNames[1]);
    intersection = eval.getMetrics(eval, "Intersection");
    bonly = eval.getMetrics(eval, eval$CallsetOnlyNames[2]);
    bonly.filtered = eval.getMetrics(eval, eval$CallsetFilteredNames[2]);

    union = list(
        all = intersection$all + aonly$all + bonly$all,
        all.withfiltered = intersection$all + aonly$all + bonly$all + aonly.filtered$all + bonly.filtered$all
    );

    plot.venn(aonly$all + intersection$all + aonly.filtered$all, bonly$all + intersection$all + bonly.filtered$all, 0, intersection$all, 0, 0, pos=c(0.32, 0.32, 0.68, 0.70), col=col);

    text(0, 0.45, cex=1.2, pos=4, .plot.callsetConcordance.getLabelText(eval$CallsetNames[1], eval$CallsetNames[2], aonly, aonly.filtered, union));
    text(0.5, 0.75, cex=1.2, adj=c(0.5, 0.33), .plot.callsetConcordance.getLabelText("Intersection", NA, intersection, NA, union));
    text(1, 0.45, cex=1.2, pos=2, .plot.callsetConcordance.getLabelText(eval$CallsetNames[2], eval$CallsetNames[1], bonly, bonly.filtered, union));
}

plot.callsetConcordanceByAC <- function(eval, normalize=TRUE, novelty_name="all", col=c("#FF6342", "#FF9675", "#5C92A4", "#88EEFF", "#55BBFF")) {
    aonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[1], novelty_name);
    aonly.filtered = eval.getMetricsByAc(eval, eval$CallsetFilteredNames[1]);
    intersection = eval.getMetricsByAc(eval, "Intersection", novelty_name);
    bonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[2], novelty_name);
    bonly.filtered = eval.getMetricsByAc(eval, eval$CallsetFilteredNames[2]);

    title = paste("Callset concordance per allele count (", novelty_name, " variants)", sep="");

    if (length(intersection$AC) > 0 && length(aonly$AC) == 0) {
        aonly = intersection;
        aonly$n = 0;
        aonly.filtered$n = 0;
    }

    if (length(intersection$AC) > 0 && length(bonly$AC) == 0) {
        bonly = intersection;
        bonly$n = 0;
        bonly.filtered$n = 0;
    }

    #par.def = par(no.readonly = TRUE);
    #par(mar=c(5, 5, 3, 5));

    if (normalize == TRUE) {
        norm = aonly$n + aonly.filtered$n + intersection$n + bonly$n + bonly.filtered$n;
        matnorm = rbind(aonly$n/norm, aonly.filtered$n/norm, intersection$n/norm, bonly.filtered$n/norm, bonly$n/norm);

        barplot(matnorm, col=col, xlab="Allele count", ylab="", main=title, names.arg=intersection$AC, xlim=c(1, 1.2*max(intersection$AC)), ylim=c(0, 1.3), border=NA, yaxt="n", cex=1.3, cex.axis=1.3, cex.lab=1.3);
        axis(2, at=seq(from=0, to=1, by=0.2), seq(from=0, to=1, by=0.2), cex=1.3, cex.axis=1.3);
        mtext("Fraction", side=2, at=0.5, padj=-3.0, cex=1.3);
    } else {
        mat = rbind(aonly$n, aonly.filtered$n, intersection$n, bonly.filtered$n, bonly$n);

        barplot(mat, col=col, xlab="Allele count", ylab="counts", main=title, names.arg=intersection$AC, xlim=c(1, max(intersection$AC)), ylim=c(0, 1), border=NA, cex=1.3, cex.axis=1.3, cex.lab=1.3);
    }

    legend(
        "topright",
        c(
            sprintf("Called in %s, absent in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
            sprintf("Called in %s, filtered in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
            "Intersection",
            sprintf("Called in %s, filtered in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2]),
            sprintf("Called in %s, absent in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2])
        ),
        fill=rev(col),
        cex=1.3
    );

    #par(par.def);
}

plot.alleleCountSpectrum <- function(eval, novelty_name="all", col=c("#FF6342", "#FF9675", "#5C92A4", "#88EEFF", "#55BBFF")) {
    aonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[1], novelty_name);
    aonly.filtered = eval.getMetricsByAc(eval, eval$CallsetFilteredNames[1]);
    intersection = eval.getMetricsByAc(eval, "Intersection", novelty_name);
    intersection.all = eval.getMetrics(eval, "Intersection");
    bonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[2], novelty_name);
    bonly.filtered = eval.getMetricsByAc(eval, eval$CallsetFilteredNames[2]);

    title = paste("Allele count spectrum (", novelty_name, " variants)", sep="");

    if (length(intersection$AC) > 0 && length(aonly$AC) == 0) {
        aonly = intersection;
        aonly$n = 0;
        aonly.filtered$n = 0;
    }

    if (length(intersection$AC) > 0 && length(bonly$AC) == 0) {
        bonly = intersection;
        bonly$n = 0;
        bonly.filtered$n = 0;
    }

    suppressWarnings(plot(0, 0, type="n", xlim=c(1, length(intersection$AC)), ylim=c(1, max(aonly$n + aonly.filtered$n + intersection$n, bonly$n + bonly.filtered$n + intersection$n)), xlab="Allele count", ylab="Number of variants", main=title, log="xy", bty="n", cex=1.3, cex.lab=1.3, cex.axis=1.3));
    suppressWarnings(points(intersection$AC, aonly$n + aonly.filtered$n + intersection$n, type="l", lwd=2, col=col[1]));
    suppressWarnings(points(intersection$AC, aonly$n + intersection$n, type="l", lwd=2, lty=2, col=col[1]));
    suppressWarnings(points(intersection$AC, intersection$n, type="l", lwd=2, col=col[3]));
    suppressWarnings(points(intersection$AC, bonly$n + intersection$n, type="l", lwd=2, lty=2, col=col[4]));
    suppressWarnings(points(intersection$AC, bonly$n + bonly.filtered$n + intersection$n, type="l", lwd=2, col=col[5]));

    loci = (unique(eval$CountVariants$nProcessedLoci))[1];

    points(c(1:max(intersection$AC)), 0.9*(1/1000)*loci*(1/c(1:max(intersection$AC))), type="l", lwd=2, lty=2, col="black");

    legend(
        "bottomleft",
        c(
            sprintf("Intersection + called in %s, absent or filtered in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
            sprintf("Intersection + called in %s, absent in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
            "Intersection",
            sprintf("Intersection + called in %s, absent in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2]),
            sprintf("Intersection + called in %s, absent or filtered in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2]),
            sprintf("Neutral expectation ( 0.9*(1/1000)*%d*(1/c(1:max(%d))) )", loci, max(intersection$AC))
        ),
        lwd=c(2, 2, 3, 2, 2, 2),
        lty=c(1, 2, 1, 2, 1, 2),
        col=c(rev(col), "black"),
        cex=1.3
    );
}

eval.getMetricsByAc <- function(eval, jexl_expression, novelty_name="all") {
    piece = eval$MetricsByAc[which(eval$MetricsByAc$evaluation_name == "eval" & eval$MetricsByAc$comparison_name == "dbsnp" & as.character(eval$MetricsByAc$jexl_expression) == as.character(jexl_expression) & eval$MetricsByAc$filter_name == "called" & eval$MetricsByAc$novelty_name == novelty_name),]; 
}

plot.titvSpectrum <- function(eval, novelty_name="all", col=c("#FF6342", "#FF9675", "#5C92A4", "#88EEFF", "#55BBFF")) {
    aonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[1], novelty_name);
    aonly.filtered = eval.getMetricsByAc(eval, eval$CallsetFilteredNames[1]);
    intersection = eval.getMetricsByAc(eval, "Intersection", novelty_name);
    bonly = eval.getMetricsByAc(eval, eval$CallsetOnlyNames[2], novelty_name);
    bonly.filtered = eval.getMetricsByAc(eval, eval$CallsetFilteredNames[2]);

    title = paste("Ti/Tv spectrum (", novelty_name, " variants)", sep="");

    if (length(intersection$AC) > 0 && length(aonly$AC) == 0) {
        aonly = intersection;
        aonly$n = 0;
        aonly$nTi = 0;
        aonly$nTv = 0;
        aonly.filtered$n = 0;
        aonly.filtered$nTi = 0;
        aonly.filtered$nTv = 0;
    }

    if (length(intersection$AC) > 0 && length(bonly$AC) == 0) {
        bonly = intersection;
        bonly$n = 0;
        bonly$nTi = 0;
        bonly$nTv = 0;
        bonly.filtered$n = 0;
        bonly.filtered$nTi = 0;
        bonly.filtered$nTv = 0;
    }

    titv.aonly.withfiltered = (aonly$nTi + aonly.filtered$nTi + intersection$nTi)/(aonly$nTv + aonly.filtered$nTv + intersection$nTv);
    titv.aonly.withfiltered.finite = titv.aonly.withfiltered[which(is.finite(titv.aonly.withfiltered))];

    titv.aonly = (aonly$nTi + intersection$nTi)/(aonly$nTv + intersection$nTv);
    titv.aonly.finite = titv.aonly[which(is.finite(titv.aonly))];

    titv.intersection.finite = intersection$Ti.Tv[which(is.finite(intersection$Ti.Tv))];

    titv.bonly = (bonly$nTi + intersection$nTi)/(bonly$nTv + intersection$nTv);
    titv.bonly.finite = titv.bonly[which(is.finite(titv.bonly))];

    titv.bonly.withfiltered = (bonly$nTi + bonly.filtered$nTi + intersection$nTi)/(bonly$nTv + bonly.filtered$nTv + intersection$nTv);
    titv.bonly.withfiltered.finite = titv.bonly.withfiltered[which(is.finite(titv.bonly.withfiltered))];

    titv.min = min(titv.aonly.withfiltered.finite, titv.aonly.finite, titv.intersection.finite, titv.bonly.finite, titv.bonly.withfiltered.finite);
    titv.max = max(titv.aonly.withfiltered.finite, titv.aonly.finite, titv.intersection.finite, titv.bonly.finite, titv.bonly.withfiltered.finite);

    plot(0, 0, type="n", xlim=c(1, length(intersection$AC)), ylim=c(0, 4), xlab="Allele count", ylab="Transition/transversion (Ti/Tv) ratio", main=title, bty="n", cex=1.3, cex.lab=1.3, cex.axis=1.3);
    points(intersection$AC, (aonly.filtered$nTi + intersection$nTi)/(aonly.filtered$nTv + intersection$nTv), type="l", lwd=2, col=col[1]);
    points(intersection$AC, (aonly$nTi + intersection$nTi)/(aonly$nTv + intersection$nTv), type="l", lwd=2, lty=2, col=col[2]);
    points(intersection$AC, intersection$Ti.Tv, type="l", lwd=2, col=col[3]);
    points(intersection$AC, (bonly$nTi + intersection$nTi)/(bonly$nTv + intersection$nTv), type="l", lwd=2, lty=2, col=col[4]);
    points(intersection$AC, (bonly.filtered$nTi + intersection$nTi)/(bonly.filtered$nTv + intersection$nTv), type="l", lwd=2, col=col[5]);

    abline(h=2.3, lty=2);
    mtext("2.3", side=4, at=2.3, cex=0.9);

    abline(h=3.3, lty=2);
    mtext("3.3", side=4, at=3.3, cex=0.9);

    #legend("topleft", c(eval$CallsetOnlyNames[1], "Intersection", eval$CallsetOnlyNames[2]), fill=col);

    legend(
        "topleft",
        c(
            sprintf("Intersection + called in %s, absent or filtered in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
            sprintf("Intersection + called in %s, absent in %s", eval$CallsetOnlyNames[2], eval$CallsetOnlyNames[1]),
            "Intersection",
            sprintf("Intersection + called in %s, absent in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2]),
            sprintf("Intersection + called in %s, absent or filtered in %s", eval$CallsetOnlyNames[1], eval$CallsetOnlyNames[2])
        ),
        lwd=c(2, 2, 3, 2, 2),
        lty=c(1, 2, 1, 2, 1),
        col=rev(col),
        cex=1.3
    );
}

plot.variantsPerSample2 <- function(eval) {
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

plot.variantsPerSample <- function(eval, novelty_name="all") {
    if (!is.na(eval$SimpleMetricsBySample)) {
        metrics = eval$SimpleMetricsBySample[which(eval$SimpleMetricsBySample$evaluation_name == "eval" & eval$SimpleMetricsBySample$comparison_name == "dbsnp" & as.character(eval$SimpleMetricsBySample$jexl_expression) == "none" & eval$SimpleMetricsBySample$filter_name == "called" & eval$SimpleMetricsBySample$novelty_name == novelty_name),]; 

        title = paste("Calls per sample (", novelty_name, ")", sep="");
        indices = order(metrics$CountVariants, decreasing=TRUE);

        par.def = par(no.readonly = TRUE);
        par(mar=c(5, 4, 4, 4));

        plot(0, 0, type="n", xaxt="n", xlim=c(1, length(metrics$row)), ylim=c(0, max(metrics$CountVariants)), xlab="", ylab="Number of variants", main=title, bty="n");
        points(c(1:length(metrics$row)), (metrics$CountVariants)[indices], pch=21, col="black");

        axis(1, at=c(1:length(metrics$row)), labels=(metrics$row)[indices], las=2, cex.axis=0.4);

        par(new=TRUE);
        plot(0, 0, type="n", xaxt="n", yaxt="n", xlim=c(1, length(metrics$row)), ylim=c(min(metrics$TiTvRatio), 1.2*max(metrics$TiTvRatio)), xlab="", ylab="", main=title, bty="n");
        points(c(1:length(metrics$row)), (metrics$TiTvRatio)[indices], pch=19, col="black");

        titvaxis = c(min(metrics$TiTvRatio), max(metrics$TiTvRatio));
        axis(4, at=titvaxis, labels=titvaxis, las=2);

        par(par.def);
    }
}
