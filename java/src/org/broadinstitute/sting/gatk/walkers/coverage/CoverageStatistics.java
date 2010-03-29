package org.broadinstitute.sting.gatk.walkers.coverage;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodRefSeq;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * A parallelizable walker designed to quickly aggregate relevant coverage statistics across samples in the input
 * file. Assesses the mean and median granular coverages of each sample, and generates part of a cumulative
 * distribution of % bases and % targets covered for certain depths. The granularity of DOC can be set by command
 * line arguments.
 *
 *
 * @Author chartl
 * @Date Feb 22, 2010
 */
// todo [DONE] -- add per target (e.g. regional) aggregation
// todo [DONE] -- add ability to print out the calculated bins and quit (for pre-analysis bin size selection)
// todo [DONE] -- refactor the location of the ALL_SAMPLE metrics [keep out of the per-sample HashMaps]
// todo [DONE] -- per locus output through -o
// todo [DONE] -- support for using read groups instead of samples
// todo [DONE] -- coverage including deletions
// todo [DONE] -- base counts
// todo -- cache the map from sample names to means in the print functions, rather than regenerating each time
// todo -- support for granular histograms for total depth; maybe n*[start,stop], bins*sqrt(n)
// todo -- alter logarithmic scaling to spread out bins more
// todo -- allow for user to set linear binning (default is logarithmic)
// todo -- formatting --> do something special for end bins in getQuantile(int[] foo), this gets mushed into the end+-1 bins for now
@By(DataSource.REFERENCE)
public class CoverageStatistics extends LocusWalker<Map<String,int[]>, CoverageAggregator> implements TreeReducible<CoverageAggregator> {
    @Argument(fullName = "start", doc = "Starting (left endpoint) for granular binning", required = false)
    int start = 1;
    @Argument(fullName = "stop", doc = "Ending (right endpoint) for granular binning", required = false)
    int stop = 500;
    @Argument(fullName = "nBins", doc = "Number of bins to use for granular binning", required = false)
    int nBins = 499;
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth. Defaults to 50.", required = false)
    byte minMappingQuality = 50;
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth. Defaults to 20.", required = false)
    byte minBaseQuality = 20;
    @Argument(fullName = "printBaseCounts", shortName = "baseCounts", doc = "Will add base counts to per-locus output.", required = false)
    boolean printBaseCounts = false;
    @Argument(fullName = "omitLocusTable", shortName = "omitLocusTable", doc = "Will not calculate the per-sample per-depth counts of loci, which should result in speedup", required = false)
    boolean omitLocusTable = false;
    @Argument(fullName = "omitIntervalStatistics", shortName = "omitIntervals", doc = "Will omit the per-interval statistics section, which should result in speedup", required = false)
    boolean omitIntervals = false;
    @Argument(fullName = "omitDepthOutputAtEachBase", shortName = "omitBaseOutput", doc = "Will omit the output of the depth of coverage at each base, which should result in speedup", required = false)
    boolean omitDepthOutput = false;
    @Argument(fullName = "printBinEndpointsAndExit", doc = "Prints the bin values and exits immediately. Use to calibrate what bins you want before running on data.", required = false)
    boolean printBinEndpointsAndExit = false;
    @Argument(fullName = "omitPerSampleStats", shortName = "omitSampleSummary", doc = "Omits the summary files per-sample. These statistics are still calculated, so this argument will not improve runtime.", required = false)
    boolean omitSampleSummary = false;
    @Argument(fullName = "useReadGroups", shortName = "rg", doc = "Split depth of coverage output by read group rather than by sample", required = false)
    boolean useReadGroup = false;
    @Argument(fullName = "useBothSampleAndReadGroup", shortName = "both", doc = "Split output by both read group and by sample", required = false)
    boolean useBoth = false;
    @Argument(fullName = "includeDeletions", shortName = "dels", doc = "Include information on deletions", required = false)
    boolean includeDeletions = false;
    @Argument(fullName = "ignoreDeletionSites", doc = "Ignore sites consisting only of deletions", required = false)
    boolean ignoreDeletionSites = false;
    @Argument(fullName = "calculateCoverageOverGenes", shortName = "geneList", doc = "Calculate the coverage statistics over this list of genes. Currently accepts RefSeq.", required = false)
    File refSeqGeneList = null;

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public boolean includeReadsWithDeletionAtLoci() { return includeDeletions && ! ignoreDeletionSites; }

    public void initialize() {

        if ( printBinEndpointsAndExit ) {
            int[] endpoints = DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins);
            System.out.print("[ ");
            for ( int e : endpoints ) {
                System.out.print(e+" ");
            }
            System.out.println("]");
            System.exit(0);
        }

        if ( getToolkit().getArguments().outFileName == null ) {
            throw new StingException("This walker requires that you specify an output file (-o)");
        }

        if ( ! omitDepthOutput ) { // print header
            out.printf("%s\t%s\t%s","Locus","Total_Depth","Average_Depth");
            // get all the samples
            HashSet<String> allSamples = getSamplesFromToolKit(useReadGroup);
            if ( useBoth ) {
                allSamples.addAll(getSamplesFromToolKit(!useReadGroup));
            }

            for ( String s : allSamples) {
                out.printf("\t%s_%s","Depth_for",s);
                if ( printBaseCounts ) {
                    out.printf("\t%s_%s",s,"base_counts");
                }
            }

            out.printf("%n");

        } else {
            out.printf("Per-Locus Depth of Coverage output was omitted");
        }
    }

    private HashSet<String> getSamplesFromToolKit( boolean getReadGroupsInstead ) {
        HashSet<String> partitions = new HashSet<String>(); // since the DOCS object uses a HashMap, this will be in the same order

        if ( getReadGroupsInstead ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partitions.add(rg.getSample()+"_rg_"+rg.getReadGroupId());
            }
        } else {
            for ( Set<String> sampleSet : getToolkit().getSamplesByReaders()) {
                for (String s : sampleSet) {
                    partitions.add(s);
                }
            }
        }

        return partitions;
    }

    public boolean isReduceByInterval() {
        return ( ! omitIntervals );
    }

    public CoverageAggregator reduceInit() {

        CoverageAggregator.AggregationType agType;
        if ( useBoth ) {
            agType = CoverageAggregator.AggregationType.BOTH;
        } else {
            agType = useReadGroup ? CoverageAggregator.AggregationType.READ : CoverageAggregator.AggregationType.SAMPLE;
        }

        CoverageAggregator aggro = new CoverageAggregator(agType,start,stop,nBins);

        if ( ! useReadGroup || useBoth ) {
            aggro.addSamples(getSamplesFromToolKit(false));
        }

        if ( useReadGroup || useBoth ) {
            aggro.addReadGroups(getSamplesFromToolKit(true));
        }

        aggro.initialize(includeDeletions,omitLocusTable);

        return aggro;
    }

    public Map<String,int[]> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if ( ! omitDepthOutput ) {
            out.printf("%s",ref.getLocus()); // yes: print locus in map, and the rest of the info in reduce (for eventual cumulatives)
            //System.out.printf("\t[log]\t%s",ref.getLocus());
        }

        Map<String,int[]> countsBySample = CoverageUtils.getBaseCountsBySample(context,minMappingQuality,minBaseQuality,
                useReadGroup ? CoverageUtils.PartitionType.BY_READ_GROUP : CoverageUtils.PartitionType.BY_SAMPLE);

        if ( useBoth ) {
            Map<String,int[]> countsByOther = CoverageUtils.getBaseCountsBySample(context,minMappingQuality,minBaseQuality,
                !useReadGroup ? CoverageUtils.PartitionType.BY_READ_GROUP : CoverageUtils.PartitionType.BY_SAMPLE);
            for ( String s : countsByOther.keySet()) {
                countsBySample.put(s,countsByOther.get(s));
            }
        }

        return countsBySample;
    }

    public CoverageAggregator reduce(Map<String,int[]> thisMap, CoverageAggregator prevReduce) {
        if ( ! omitDepthOutput ) {
            printDepths(out,thisMap, prevReduce.getAllSamples());
            // this is an additional iteration through thisMap, plus dealing with IO, so should be much slower without
            // turning on omit
        }

        prevReduce.update(thisMap); // note that in "useBoth" cases, this method alters the thisMap object

        return prevReduce;
    }

    public CoverageAggregator treeReduce(CoverageAggregator left, CoverageAggregator right) {
        left.merge(right);
        return left;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // INTERVAL ON TRAVERSAL DONE
    ////////////////////////////////////////////////////////////////////////////////////

    public void onTraversalDone( List<Pair<GenomeLoc,CoverageAggregator>> statsByInterval ) {
        if ( refSeqGeneList != null && ( useBoth || ! useReadGroup ) ) {
            printGeneStats(statsByInterval);
        }

        if ( useBoth || ! useReadGroup ) {
            File intervalStatisticsFile = deriveFromStream("sample_interval_statistics");
            File intervalSummaryFile = deriveFromStream("sample_interval_summary");
            printIntervalStats(statsByInterval,intervalSummaryFile, intervalStatisticsFile, true);
        }

        if ( useBoth || useReadGroup ) {
            File intervalStatisticsFile = deriveFromStream("read_group_interval_statistics");
            File intervalSummaryFile = deriveFromStream("read_group_interval_summary");
            printIntervalStats(statsByInterval,intervalSummaryFile, intervalStatisticsFile, false);
        }

        onTraversalDone(mergeAll(statsByInterval));

    }

    public CoverageAggregator mergeAll(List<Pair<GenomeLoc,CoverageAggregator>> stats) {
        CoverageAggregator first = stats.remove(0).second;
        for ( Pair<GenomeLoc,CoverageAggregator> iStat : stats ) {
            treeReduce(first,iStat.second);
        }

        return first;
    }

    private DepthOfCoverageStats printIntervalStats(List<Pair<GenomeLoc,CoverageAggregator>> statsByInterval, File summaryFile, File statsFile, boolean isSample) {
        PrintStream summaryOut;
        PrintStream statsOut;

        try {
            summaryOut = summaryFile == null ? out : new PrintStream(summaryFile);
            statsOut = statsFile == null ? out : new PrintStream(statsFile);
        } catch ( IOException e ) {
            throw new StingException("Unable to open interval file on reduce", e);
        }


        Pair<GenomeLoc,CoverageAggregator> firstPair = statsByInterval.get(0);
        CoverageAggregator firstAggregator = firstPair.second;
        DepthOfCoverageStats firstStats = isSample ? firstAggregator.getCoverageBySample() : firstAggregator.getCoverageByReadGroup();

        StringBuilder summaryHeader = new StringBuilder();
        summaryHeader.append("Target");
        summaryHeader.append("\ttotal_coverage");
        summaryHeader.append("\taverage_coverage");

        for ( String s : firstStats.getAllSamples() ) {
            summaryHeader.append("\t");
            summaryHeader.append(s);
            summaryHeader.append("_total_cvg");
            summaryHeader.append("\t");
            summaryHeader.append(s);
            summaryHeader.append("_mean_cvg");
            summaryHeader.append("\t");
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q1");
            summaryHeader.append("\t");
            summaryHeader.append(s);
            summaryHeader.append("_granular_median");
            summaryHeader.append("\t");
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q3");
        }

        summaryOut.printf("%s%n",summaryHeader);

        int[][] nTargetsByAvgCvgBySample = new int[firstStats.getHistograms().size()][firstStats.getEndpoints().length+1];

        for ( Pair<GenomeLoc,CoverageAggregator> targetAggregator : statsByInterval ) {

            Pair<GenomeLoc,DepthOfCoverageStats> targetStats = new Pair<GenomeLoc,DepthOfCoverageStats>(
                    targetAggregator.first, isSample ? targetAggregator.second.getCoverageBySample() :
                    targetAggregator.second.getCoverageByReadGroup());
            printTargetSummary(summaryOut,targetStats);
            updateTargetTable(nTargetsByAvgCvgBySample,targetStats.second);
        }

        printIntervalTable(statsOut,nTargetsByAvgCvgBySample,firstStats.getEndpoints());

        if ( ! getToolkit().getArguments().outFileName.contains("stdout")) {
            summaryOut.close();
            statsOut.close();
        }

        return firstStats;
    }

    private void printGeneStats(List<Pair<GenomeLoc,CoverageAggregator>> statsByTarget) {
        LocationAwareSeekableRODIterator refseqIterator = initializeRefSeq();
        List<Pair<String,DepthOfCoverageStats>> statsByGene = new ArrayList<Pair<String,DepthOfCoverageStats>>();// maintains order
        Map<String,DepthOfCoverageStats> geneNamesToStats = new HashMap<String,DepthOfCoverageStats>(); // alows indirect updating of objects in list

        for ( Pair<GenomeLoc,CoverageAggregator> targetStats : statsByTarget ) {
            String gene = getGeneName(targetStats.first,refseqIterator);
            if ( geneNamesToStats.keySet().contains(gene) ) {
                geneNamesToStats.get(gene).merge(targetStats.second.getCoverageBySample());
            } else {
                geneNamesToStats.put(gene,targetStats.second.getCoverageBySample());
                statsByGene.add(new Pair<String,DepthOfCoverageStats>(gene,targetStats.second.getCoverageBySample()));
            }
        }

        PrintStream geneSummaryOut = getCorrectStream(out,deriveFromStream("gene_summary"));

        for ( Pair<String,DepthOfCoverageStats> geneStats : statsByGene ) {
            printTargetSummary(geneSummaryOut,geneStats);
        }

        if ( ! getToolkit().getArguments().outFileName.contains("stdout")) {
            geneSummaryOut.close();
        }

    }

    //blatantly stolen from Andrew Kernytsky
    private String getGeneName(GenomeLoc target, LocationAwareSeekableRODIterator refseqIterator) {
        if (refseqIterator == null) { return "UNKNOWN"; }

        RODRecordList annotationList = refseqIterator.seekForward(target);
        if (annotationList == null) { return "UNKNOWN"; }

        for(ReferenceOrderedDatum rec : annotationList) {
            if ( ((rodRefSeq)rec).overlapsExonP(target) ) {
                return ((rodRefSeq)rec).getGeneName();
            }
        }

        return "UNKNOWN";

    }

    private LocationAwareSeekableRODIterator initializeRefSeq() {
        ReferenceOrderedData<rodRefSeq> refseq = new ReferenceOrderedData<rodRefSeq>("refseq",
                refSeqGeneList, rodRefSeq.class);
        return refseq.iterator();
    }

    private void printTargetSummary(PrintStream output, Pair<?,DepthOfCoverageStats> intervalStats) {
        DepthOfCoverageStats stats = intervalStats.second;
        int[] bins = stats.getEndpoints();
        StringBuilder targetSummary = new StringBuilder();
        targetSummary.append(intervalStats.first.toString());
        targetSummary.append("\t");
        targetSummary.append(stats.getTotalCoverage());
        targetSummary.append("\t");
        targetSummary.append(String.format("%.2f",stats.getTotalMeanCoverage()));

        for ( String s : stats.getAllSamples() ) {
            targetSummary.append("\t");
            targetSummary.append(stats.getTotals().get(s));
            targetSummary.append("\t");
            targetSummary.append(String.format("%.2f", stats.getMeans().get(s)));
            targetSummary.append("\t");
            int median = getQuantile(stats.getHistograms().get(s),0.5);
            int q1 = getQuantile(stats.getHistograms().get(s),0.25);
            int q3 = getQuantile(stats.getHistograms().get(s),0.75);
            targetSummary.append(formatBin(bins,q1));
            targetSummary.append("\t");
            targetSummary.append(formatBin(bins,median));
            targetSummary.append("\t");
            targetSummary.append(formatBin(bins,q3));

        }

        output.printf("%s%n", targetSummary);
    }

    private String formatBin(int[] bins, int quartile) {
        if ( quartile >= bins.length ) {
            return String.format(">%d",bins[bins.length-1]);
        } else if ( quartile < 0 ) {
            return String.format("<%d",bins[0]);
        } else {
            return String.format("%d",bins[quartile]);
        }
    }

    private void printIntervalTable(PrintStream output, int[][] intervalTable, int[] cutoffs) {
        output.printf("\tdepth>=%d",0);
        for ( int col = 0; col < intervalTable[0].length-1; col ++ ) {
            output.printf("\tdepth>=%d",cutoffs[col]);
        }

        output.printf(String.format("%n"));
        for ( int row = 0; row < intervalTable.length; row ++ ) {
            output.printf("At_least_%d_samples",row+1);
            for ( int col = 0; col < intervalTable[0].length; col++ ) {
                output.printf("\t%d",intervalTable[row][col]);
            }
            output.printf(String.format("%n"));
        }
    }

    /*
     * @updateTargetTable
     * The idea is to have counts for how many *targets* have at least K samples with
     * median coverage of at least X.
     * To that end:
     * Iterate over the samples the DOCS object, determine how many there are with
     * median coverage > leftEnds[0]; how many with median coverage > leftEnds[1]
     * and so on. Then this target has at least N, N-1, N-2, ... 1, 0 samples covered
     * to leftEnds[0] and at least M,M-1,M-2,...1,0 samples covered to leftEnds[1]
     * and so on.
     */
    private void updateTargetTable(int[][] table, DepthOfCoverageStats stats) {
        int[] cutoffs = stats.getEndpoints();
        int[] countsOfMediansAboveCutoffs = new int[cutoffs.length+1]; // 0 bin to catch everything
        for ( int i = 0; i < countsOfMediansAboveCutoffs.length; i ++) {
            countsOfMediansAboveCutoffs[i]=0;
        }

        for ( String s : stats.getAllSamples() ) {
            int medianBin = getQuantile(stats.getHistograms().get(s),0.5);
            for ( int i = 0; i <= medianBin; i ++) {
                countsOfMediansAboveCutoffs[i]++;
            }
        }

        for ( int medianBin = 0; medianBin < countsOfMediansAboveCutoffs.length; medianBin++) {
            for ( ; countsOfMediansAboveCutoffs[medianBin] > 0; countsOfMediansAboveCutoffs[medianBin]-- ) {
                table[countsOfMediansAboveCutoffs[medianBin]-1][medianBin]++;
                // the -1 is due to counts being 1-based and offsets being 0-based
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // FINAL ON TRAVERSAL DONE
    ////////////////////////////////////////////////////////////////////////////////////

    public void onTraversalDone(CoverageAggregator coverageProfiles) {
        ///////////////////
        // OPTIONAL OUTPUTS
        //////////////////

        if ( ! omitSampleSummary ) {
            logger.info("Printing summary info");
            if ( ! useReadGroup || useBoth ) {
                File summaryStatisticsFile = deriveFromStream("sample_summary_statistics");
                File perSampleStatisticsFile = deriveFromStream("sample_statistics");
                printSummary(out,summaryStatisticsFile,coverageProfiles.getCoverageBySample());
                printPerSample(out,perSampleStatisticsFile,coverageProfiles.getCoverageBySample());
            }

            if ( useReadGroup || useBoth ) {
                File rgStatsFile = deriveFromStream("read_group_summary");
                File rgSumFile = deriveFromStream("read_group_statistics");
                printSummary(out,rgStatsFile,coverageProfiles.getCoverageByReadGroup());
                printPerSample(out,rgSumFile,coverageProfiles.getCoverageByReadGroup());
            }
        }

        if ( ! omitLocusTable ) {
            logger.info("Printing locus summary");
            if ( ! useReadGroup || useBoth ) {
                File perLocusStatisticsFile = deriveFromStream("sample_locus_statistics");
                printPerLocus(perLocusStatisticsFile,coverageProfiles.getCoverageBySample());
            }

            if ( useReadGroup || useBoth ) {
                File perLocusRGStats = deriveFromStream("read_group_locus_statistics");
                printPerLocus(perLocusRGStats,coverageProfiles.getCoverageByReadGroup());
            }
        }
    }

    public File deriveFromStream(String append) {
        String name = getToolkit().getArguments().outFileName;
        if ( name.contains("stdout") || name.contains("Stdout") || name.contains("STDOUT")) {
            return null;
        } else {
            return new File(name+"."+append);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // HELPER OUTPUT METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    private void printPerSample(PrintStream out, File optionalFile, DepthOfCoverageStats stats) {
        PrintStream output = getCorrectStream(out,optionalFile);
        int[] leftEnds = stats.getEndpoints();
        StringBuilder hBuilder = new StringBuilder();
        hBuilder.append("\t");
        hBuilder.append(String.format("from_0_to_%d)\t",leftEnds[0]));
        for ( int i = 1; i < leftEnds.length; i++ )
            hBuilder.append(String.format("from_%d_to_%d)\t",leftEnds[i-1],leftEnds[i]));
        hBuilder.append(String.format("from_%d_to_inf%n",leftEnds[leftEnds.length-1]));
        output.print(hBuilder.toString());
        Map<String,int[]> histograms = stats.getHistograms();
        for ( String s : histograms.keySet() ) {
            StringBuilder sBuilder = new StringBuilder();
            sBuilder.append(String.format("sample_%s",s));
            for ( int count : histograms.get(s) ) {
                sBuilder.append(String.format("\t%d",count));
            }
            sBuilder.append(String.format("%n"));
            output.print(sBuilder.toString());
        }
    }

    private void printPerLocus(File locusFile, DepthOfCoverageStats stats) {
        PrintStream output = getCorrectStream(out,locusFile);
        if ( output == null ) {
            return;
        }

        int[] endpoints = stats.getEndpoints();
        int samples = stats.getHistograms().size();

        int[][] baseCoverageCumDist = stats.getLocusCounts();

        // rows - # of samples
        // columns - depth of coverage

        StringBuilder header = new StringBuilder();
        header.append(String.format("\tgte_0"));
        for ( int d : endpoints ) {
            header.append(String.format("\tgte_%d",d));
        }
        header.append(String.format("%n"));

        output.print(header);

        for ( int row = 0; row < samples; row ++ ) {
            output.printf("%s_%d\t","NSamples",row+1);
            for ( int depthBin = 0; depthBin < baseCoverageCumDist[0].length; depthBin ++ ) {
                output.printf("%d\t",baseCoverageCumDist[row][depthBin]);
            }
            output.printf("%n");
        }
    }

    private PrintStream getCorrectStream(PrintStream out, File optionalFile) {
        PrintStream output;
        if ( optionalFile == null ) {
            output = out;
        } else {
            try {
                output = new PrintStream(optionalFile);
            } catch ( IOException e ) {
                logger.warn("Error opening the output file "+optionalFile.getAbsolutePath()+". Defaulting to stdout");
                output = out;
            }
        }

        return output;
    }

    private void printSummary(PrintStream out, File optionalFile, DepthOfCoverageStats stats) {
        PrintStream output = getCorrectStream(out,optionalFile);

        output.printf("%s\t%s\t%s\t%s\t%s\t%s%n","sample_id","total","mean","granular_third_quartile","granular_median","granular_first_quartile");

        Map<String,int[]> histograms = stats.getHistograms();
        Map<String,Double> means = stats.getMeans();
        Map<String,Long> totals = stats.getTotals();
        int[] leftEnds = stats.getEndpoints();

        for ( String s : histograms.keySet() ) {
            int[] histogram = histograms.get(s);
            int median = getQuantile(histogram,0.5);
            int q1 = getQuantile(histogram,0.25);
            int q3 = getQuantile(histogram,0.75);
            // if any of these are larger than the higest bin, put the median as in the largest bin
            median =  median == histogram.length-1 ? histogram.length-2 : median;
            q1 = q1 == histogram.length-1 ? histogram.length-2 : q1;
            q3 = q3 == histogram.length-1 ? histogram.length-2 : q3;
            output.printf("%s\t%d\t%.2f\t%d\t%d\t%d%n",s,totals.get(s),means.get(s),leftEnds[q3],leftEnds[median],leftEnds[q1]);
        }

        output.printf("%s\t%d\t%.2f\t%s\t%s\t%s%n","Total",stats.getTotalCoverage(),stats.getTotalMeanCoverage(),"N/A","N/A","N/A");
    }

    private int getQuantile(int[] histogram, double prop) {
        int total = 0;

        for ( int i = 0; i < histogram.length; i ++ ) {
            total += histogram[i];
        }

        int counts = 0;
        int bin = -1;
        while ( counts < prop*total ) {
            counts += histogram[bin+1];
            bin++;
        }

        return bin;
    }

    private void printDepths(PrintStream stream, Map<String,int[]> countsBySample, Set<String> allSamples) {
        // get the depths per sample and build up the output string while tabulating total and average coverage
        // todo -- update me to deal with base counts/indels
        StringBuilder perSampleOutput = new StringBuilder();
        int tDepth = 0;
        for ( String s : allSamples ) {
            perSampleOutput.append("\t");
            long dp = countsBySample.keySet().contains(s) ? sumArray(countsBySample.get(s)) : 0;
            perSampleOutput.append(dp);
            if ( printBaseCounts ) {
                perSampleOutput.append("\t");
                perSampleOutput.append(baseCounts(countsBySample.get(s)));
            }
            tDepth += dp;
        }
        // remember -- genome locus was printed in map()
        stream.printf("\t%d\t%.2f\t%s%n",tDepth,( (double) tDepth/ (double) allSamples.size()), perSampleOutput);
        //System.out.printf("\t%d\t%.2f\t%s%n",tDepth,( (double) tDepth/ (double) allSamples.size()), perSampleOutput);

    }

    private long sumArray(int[] array) {
        long i = 0;
        for ( int j : array ) {
            i += j;
        }
        return i;
    }

    private String baseCounts(int[] counts) {
        if ( counts == null ) {
            counts = new int[6];
        }
        StringBuilder s = new StringBuilder();
        int nbases = 0;
        for ( char b : BaseUtils.EXTENDED_BASES ) {
            nbases++;
            if ( includeDeletions || b != 'D' ) {
                s.append(b);
                s.append(":");
                s.append(counts[BaseUtils.extendedBaseToBaseIndex(b)]);
                if ( nbases < 6 ) {
                    s.append(" ");
                }
            }
        }

        return s.toString();
    }
}

class CoverageAggregator {
    private DepthOfCoverageStats coverageByRead;
    private DepthOfCoverageStats coverageBySample;

    enum AggregationType { READ, SAMPLE, BOTH }

    private AggregationType agType;
    private Set<String> sampleNames;
    private Set<String> allSamples;

    public CoverageAggregator(AggregationType type, int start, int stop, int nBins) {
        if ( type == AggregationType.READ || type == AggregationType.BOTH) {
            coverageByRead = new DepthOfCoverageStats(DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins));
        }

        if ( type == AggregationType.SAMPLE || type == AggregationType.BOTH) {
            coverageBySample = new DepthOfCoverageStats(DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins));
        }

        agType = type;
        allSamples = new HashSet<String>();
    }

    public void merge(CoverageAggregator otherAggregator) {
        if ( agType == AggregationType.SAMPLE || agType == AggregationType.BOTH ) {
            this.coverageBySample.merge(otherAggregator.coverageBySample);
        }

        if ( agType == AggregationType.READ || agType == AggregationType.BOTH) {
            this.coverageByRead.merge(otherAggregator.coverageByRead);
        }
    }

    public DepthOfCoverageStats getCoverageByReadGroup() {
        return coverageByRead;
    }

    public DepthOfCoverageStats getCoverageBySample() {
        return coverageBySample;
    }

    public void addSamples(Set<String> samples) {
        for ( String s : samples ) {
            coverageBySample.addSample(s);
            allSamples.add(s);
        }

        if ( agType == AggregationType.BOTH ) {
            sampleNames = samples;
        }
    }

    public void addReadGroups(Set<String> readGroupNames) {
        for ( String s : readGroupNames ) {
            coverageByRead.addSample(s);
            allSamples.add(s);
        }
    }


    public void initialize(boolean useDels, boolean omitLocusTable) {
        if ( agType == AggregationType.SAMPLE || agType == AggregationType.BOTH ) {
            if ( useDels ) {
                coverageBySample.initializeDeletions();
            }

            if ( ! omitLocusTable ) {
                coverageBySample.initializeLocusCounts();
            }
        }

        if ( agType == AggregationType.READ || agType == AggregationType.BOTH) {
            if ( useDels ) {
                coverageByRead.initializeDeletions();
            }

            if ( ! omitLocusTable ) {
                coverageByRead.initializeLocusCounts();
            }
        }
    }

    public void update(Map<String,int[]> countsByIdentifier) {
        if ( agType != AggregationType.BOTH ) {
            if ( agType == AggregationType.SAMPLE ) {
                coverageBySample.update(countsByIdentifier);
            } else {
                coverageByRead.update(countsByIdentifier);
            }
        } else { // have to separate samples and read groups
            HashMap<String,int[]> countsBySample = new HashMap<String,int[]>(sampleNames.size());
            HashMap<String,int[]> countsByRG = new HashMap<String,int[]>(allSamples.size()-sampleNames.size());
            for ( String s : countsByIdentifier.keySet() ) {
                if ( sampleNames.contains(s) ) {
                    countsBySample.put(s,countsByIdentifier.get(s));
                } else {                                               // cannot use .remove() to save time due to concurrency issues
                    countsByRG.put(s,countsByIdentifier.get(s));
                }
            }

            coverageBySample.update(countsBySample);
            coverageByRead.update(countsByRG);
        }
    }

    public Set<String> getAllSamples() {
        return allSamples;
    }

}