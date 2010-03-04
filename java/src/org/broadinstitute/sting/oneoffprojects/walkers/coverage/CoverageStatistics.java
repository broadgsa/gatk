package org.broadinstitute.sting.oneoffprojects.walkers.coverage;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
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
import org.broadinstitute.sting.gatk.walkers.coverage.DepthOfCoverageWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.pileup.PileupElement;

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
public class CoverageStatistics extends LocusWalker<Map<String,int[]>, DepthOfCoverageStats> implements TreeReducible<DepthOfCoverageStats> {
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
    @Argument(fullName = "omitLocusTable", shortName = "omitLocus", doc = "Will not calculate the per-sample per-depth counts of loci, which should result in speedup", required = false)
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
            //System.out.printf("\t[log]\t%s\t%s\t%s","Locus","Total_Depth","Average_Depth");
            // get all the samples
            HashSet<String> allSamples = getSamplesFromToolKit(useReadGroup);

            for ( String s : allSamples) {
                out.printf("\t%s_%s","Depth_for",s);
                //System.out.printf("\t%s_%s","Depth_for",s);
                if ( printBaseCounts ) {
                    out.printf("\t%s_%s",s,"base_counts");
                    //System.out.printf("\t%s_%s",s,"base_counts");
                }
            }

            out.printf("%n");
            //System.out.printf("%n");

        } else {
            out.printf("Per-Locus Depth of Coverage output was omitted");
        }
    }
    
    private HashSet<String> getSamplesFromToolKit( boolean getReadGroupsInstead ) {
        HashSet<String> partitions = new HashSet<String>(); // since the DOCS object uses a HashMap, this will be in the same order

        if ( getReadGroupsInstead ) {
            for ( Set<String> rgSet : getToolkit().getMergedReadGroupsByReaders() ) {
                for ( String rg : rgSet ) {
                    partitions.add(rg);
                }
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

    public DepthOfCoverageStats reduceInit() {
        DepthOfCoverageStats stats = new DepthOfCoverageStats(DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins));

        for ( String sample : getSamplesFromToolKit(useReadGroup) ) {
            stats.addSample(sample);
        }


        if ( ! omitLocusTable ) {
            stats.initializeLocusCounts();
        }

        if ( includeDeletions ) {
            stats.initializeDeletions();
        }

        return stats;
    }

    public Map<String,int[]> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        
        if ( ! omitDepthOutput ) {
            out.printf("%s",ref.getLocus()); // yes: print locus in map, and the rest of the info in reduce (for eventual cumulatives)
            //System.out.printf("\t[log]\t%s",ref.getLocus());
        }

        Map<String,int[]> countsBySample = CoverageUtils.getBaseCountsBySample(context,minMappingQuality,minBaseQuality,
                useReadGroup ? CoverageUtils.PartitionType.BY_READ_GROUP : CoverageUtils.PartitionType.BY_SAMPLE);

        return countsBySample;
    }

    public DepthOfCoverageStats reduce(Map<String,int[]> thisMap, DepthOfCoverageStats prevReduce) {
        prevReduce.update(thisMap);
        if ( ! omitDepthOutput ) {
            printDepths(out,thisMap, prevReduce.getAllSamples());
            // this is an additional iteration through thisMap, plus dealing with IO, so should be much slower without
            // turning on omit
        }

        return prevReduce;
    }

    public DepthOfCoverageStats treeReduce(DepthOfCoverageStats left, DepthOfCoverageStats right) {
        left.merge(right);
        return left;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // INTERVAL ON TRAVERSAL DONE
    ////////////////////////////////////////////////////////////////////////////////////

    public void onTraversalDone( List<Pair<GenomeLoc,DepthOfCoverageStats>> statsByInterval ) {
        if ( refSeqGeneList != null ) {
            printGeneStats(statsByInterval);
        }
        File intervalStatisticsFile = deriveFromStream("interval_statistics");
        File intervalSummaryFile = deriveFromStream("interval_summary");
        DepthOfCoverageStats mergedStats = printIntervalStatsAndMerge(statsByInterval,intervalSummaryFile, intervalStatisticsFile);
        this.onTraversalDone(mergedStats);
    }

    private DepthOfCoverageStats printIntervalStatsAndMerge(List<Pair<GenomeLoc,DepthOfCoverageStats>> statsByInterval, File summaryFile, File statsFile) {
        PrintStream summaryOut;
        PrintStream statsOut;

        try {
            summaryOut = summaryFile == null ? out : new PrintStream(summaryFile);
            statsOut = statsFile == null ? out : new PrintStream(statsFile);
        } catch ( IOException e ) {
            throw new StingException("Unable to open interval file on reduce", e);
        }


        Pair<GenomeLoc,DepthOfCoverageStats> firstPair = statsByInterval.remove(0);
        DepthOfCoverageStats firstStats = firstPair.second;

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
        for ( int i = 0; i < nTargetsByAvgCvgBySample.length; i ++ ) {
            for ( int b = 0; b < nTargetsByAvgCvgBySample[0].length; b++) {
                nTargetsByAvgCvgBySample[i][b] = 0;
            }
        }

        printTargetSummary(summaryOut,firstPair);
        updateTargetTable(nTargetsByAvgCvgBySample,firstStats);

        for ( Pair<GenomeLoc,DepthOfCoverageStats> targetStats : statsByInterval ) {
            printTargetSummary(summaryOut,targetStats);
            updateTargetTable(nTargetsByAvgCvgBySample,targetStats.second);
            firstStats = this.treeReduce(firstStats,targetStats.second);
        }

        printIntervalTable(statsOut,nTargetsByAvgCvgBySample,firstStats.getEndpoints());

        if ( ! getToolkit().getArguments().outFileName.contains("stdout")) {
            summaryOut.close();
            statsOut.close();
        }

        return firstStats;
    }

    private void printGeneStats(List<Pair<GenomeLoc,DepthOfCoverageStats>> statsByTarget) {
        LocationAwareSeekableRODIterator refseqIterator = initializeRefSeq();
        List<Pair<String,DepthOfCoverageStats>> statsByGene = new ArrayList<Pair<String,DepthOfCoverageStats>>();// maintains order
        Map<String,DepthOfCoverageStats> geneNamesToStats = new HashMap<String,DepthOfCoverageStats>(); // alows indirect updating of objects in list

        for ( Pair<GenomeLoc,DepthOfCoverageStats> targetStats : statsByTarget ) {
            String gene = getGeneName(targetStats.first,refseqIterator);
            if ( geneNamesToStats.keySet().contains(gene) ) {
                geneNamesToStats.get(gene).merge(targetStats.second);
            } else {
                geneNamesToStats.put(gene,targetStats.second);
                statsByGene.add(new Pair<String,DepthOfCoverageStats>(gene,targetStats.second));
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
            targetSummary.append(bins[q1]);
            targetSummary.append("\t");
            targetSummary.append(bins[median]);
            targetSummary.append("\t");
            targetSummary.append(bins[q3]);

        }

        output.printf("%s%n", targetSummary);
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

    public void onTraversalDone(DepthOfCoverageStats coverageProfiles) {
        ///////////////////
        // OPTIONAL OUTPUTS
        //////////////////
        
        if ( ! omitSampleSummary ) {
            logger.info("Printing sample summary");
            File summaryStatisticsFile = deriveFromStream("summary_statistics");
            File perSampleStatisticsFile = deriveFromStream("sample_statistics");
            printSummary(out,summaryStatisticsFile,coverageProfiles);
            printPerSample(out,perSampleStatisticsFile,coverageProfiles);
        }

        if ( ! omitLocusTable ) {
            logger.info("Printing locus summary");
            File perLocusStatisticsFile = deriveFromStream("locus_statistics");
            printPerLocus(perLocusStatisticsFile,coverageProfiles);
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
        hBuilder.append(String.format("[0,%d)\t",leftEnds[0]));
        for ( int i = 1; i < leftEnds.length; i++ )
            hBuilder.append(String.format("[%d,%d)\t",leftEnds[i-1],leftEnds[i]));
        hBuilder.append(String.format("[%d,inf)%n",leftEnds[leftEnds.length-1]));
        output.print(hBuilder.toString());
        Map<String,int[]> histograms = stats.getHistograms();
        for ( String s : histograms.keySet() ) {
            StringBuilder sBuilder = new StringBuilder();
            sBuilder.append(String.format("%s",s));
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
        header.append(String.format("\t>=0"));
        for ( int d : endpoints ) {
            header.append(String.format("\t>=%d",d));
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