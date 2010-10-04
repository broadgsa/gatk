.attemptToLoadFile <- function(filename) {
    file = NA;

    if (file.exists(filename) & file.info(filename)$size > 500) {
        file = read.csv(filename, header=TRUE, comment.char="#");
    }

    file;
}

read.eval <-
function(evalRoot) {
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

