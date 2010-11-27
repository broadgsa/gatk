\name{gsa.read.eval}
\alias{gsa.read.eval}
\title{
Read a VariantEval file
}
\description{
Read a VariantEval file that's output in R format.
}
\usage{
gsa.read.eval(evalRoot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{evalRoot}{
%%     ~~Describe \code{evalRoot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
    eval$CallsetFilteredNames = as.vector(c()); 
    eval;
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
