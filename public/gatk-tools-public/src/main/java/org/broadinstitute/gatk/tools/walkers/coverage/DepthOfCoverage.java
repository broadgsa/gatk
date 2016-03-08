/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.coverage;

import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.SeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrack;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.codecs.refseq.RefSeqCodec;
import org.broadinstitute.gatk.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library
 *
 * <p>
 * This tool processes a set of bam files to determine coverage at different levels of partitioning and
 * aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by
 * sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles,
 * and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by
 * mapping or base quality score.
 * </p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>One or more bam files (with proper headers) to be analyzed for coverage statistics</li>
 *     <li>(Optional) A REFSEQ file to aggregate coverage to the gene level (for information about creating the REFSEQ Rod, please consult the online documentation)</li>
 * </ul>

 * <h3>Output</h3>
 * <p>
 * Tables pertaining to different coverage summaries. Suffix on the table files declares the contents:
 * </p>
 * <ul>
 *     <li>no suffix: per locus coverage</li>
 *     <li>_summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases</li>
 *     <li>_statistics: coverage histograms (# locus with X coverage), aggregated over all bases</li>
 *     <li>_interval_summary: total, mean, median, quartiles, and threshold proportions, aggregated per interval</li>
 *     <li>_interval_statistics: 2x2 table of # of intervals covered to >= X depth in >=Y samples</li>
 *     <li>_gene_summary: total, mean, median, quartiles, and threshold proportions, aggregated per gene</li>
 *     <li>_gene_statistics: 2x2 table of # of genes covered to >= X depth in >= Y samples</li>
 *     <li>_cumulative_coverage_counts: coverage histograms (# locus with >= X coverage), aggregated over all bases</li>
 *     <li>_cumulative_coverage_proportions: proprotions of loci with >= X coverage, aggregated over all bases</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T DepthOfCoverage \
 *   -R reference.fasta \
 *   -o file_name_base \
 *   -I input_bams.list
 *   [-geneList refSeq.sorted.txt] \
 *   [-pt readgroup] \
 *   [-ct 4 -ct 6 -ct 10] \
 *   [-L my_capture_genes.interval_list]
 * </pre>
 */
// todo -- cache the map from sample names to means in the print functions, rather than regenerating each time
// todo -- support for granular histograms for total depth; maybe n*[start,stop], bins*sqrt(n)
// todo -- alter logarithmic scaling to spread out bins more
// todo -- allow for user to set linear binning (default is logarithmic)
// todo -- formatting --> do something special for end bins in getQuantile(int[] foo), this gets mushed into the end+-1 bins for now
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class}, gotoDev = HelpConstants.MC)
@By(DataSource.REFERENCE)
@PartitionBy(PartitionType.NONE)
@Downsample(by= DownsampleType.NONE, toCoverage=Integer.MAX_VALUE)
public class DepthOfCoverage extends LocusWalker<Map<DoCOutputType.Partition,Map<String,int[]>>, CoveragePartitioner> implements TreeReducible<CoveragePartitioner> {
    private final static Logger logger = Logger.getLogger(DepthOfCoverage.class);

    /**
     * Warning message for when the incompatible arguments --calculateCoverageOverGenes and --omitIntervalStatistics are used together.
     */
    private static final String incompatibleArgsMsg = "The arguments --calculateCoverageOverGenes and --omitIntervalStatistics are incompatible. Using them together will result in an empty gene summary output file.";

    @Output
    @Multiplex(value=DoCOutputMultiplexer.class,arguments={"partitionTypes","refSeqGeneList","omitDepthOutput","omitIntervals","omitSampleSummary","omitLocusTable"})
    Map<DoCOutputType,PrintStream> out;
    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    /**
     * Reads with mapping quality values higher than this threshold will be skipped. The default value is the largest number that can be represented as an integer by the program.
     */
    @Argument(fullName = "maxMappingQuality", doc = "Maximum mapping quality of reads to count towards depth", required = false, minValue = 0, maxValue = Integer.MAX_VALUE)
    int maxMappingQuality = Integer.MAX_VALUE;
    /**
     * Bases with quality scores lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth", required = false, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte minBaseQuality = -1;
    /**
     * Bases with quality scores higher than this threshold will be skipped. The default value is the largest number that can be represented as a byte.
     */
    @Argument(fullName = "maxBaseQuality", doc = "Maximum quality of bases to count towards depth", required = false, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte maxBaseQuality = Byte.MAX_VALUE;

    @Argument(fullName = "countType", doc = "How should overlapping reads from the same fragment be handled?", required = false)
    CoverageUtils.CountPileupType countType = CoverageUtils.CountPileupType.COUNT_READS;

    /**
     * Instead of reporting depth, the program will report the base pileup at each locus
     */
    @Argument(fullName = "printBaseCounts", shortName = "baseCounts", doc = "Add base counts to per-locus output", required = false)
    boolean printBaseCounts = false;

    /**
     * Disabling the tabulation of locus statistics (# loci covered by sample by coverage) should speed up processing.
     */
    @Argument(fullName = "omitLocusTable", shortName = "omitLocusTable", doc = "Do not calculate per-sample per-depth counts of loci", required = false)
    boolean omitLocusTable = false;

    /**
     * Disabling the tabulation of interval statistics (mean, median, quartiles AND # intervals by sample by coverage) should speed up processing. This option is required in order to use -nt parallelism.
     */
    @Argument(fullName = "omitIntervalStatistics", shortName = "omitIntervals", doc = "Do not calculate per-interval statistics", required = false)
    boolean omitIntervals = false;
    /**
     * Disabling the tabulation of total coverage at every base should speed up processing.
     */
    @Argument(fullName = "omitDepthOutputAtEachBase", shortName = "omitBaseOutput", doc = "Do not output depth of coverage at each base", required = false)
    boolean omitDepthOutput = false;

    /**
     * Specify a RefSeq file for use in aggregating coverage statistics over genes.
     *
     * This argument is incompatible with --calculateCoverageOverGenes and --omitIntervalStatistics. A warning will be logged and no output file will be produced for the gene list if these arguments are enabled together.
     *
     */
    @Argument(fullName = "calculateCoverageOverGenes", shortName = "geneList", doc = "Calculate coverage statistics over this list of genes", required = false)
    File refSeqGeneList = null;

    /**
     * Output file format (e.g. csv, table, rtable); defaults to r-readable table.
     */
    @Argument(fullName = "outputFormat", doc = "The format of the output file", required = false)
    String outputFormat = "rtable";


    // ---------------------------------------------------------------------------
    //
    // Advanced arguments
    //
    // ---------------------------------------------------------------------------

    /**
     * Normally, sites where the reference is N (or another non-canonical base) are skipped. If this option is enabled, these sites will be included in DoC calculations if there is coverage from neighboring reads.
     */
    @Advanced
    @Argument(fullName = "includeRefNSites", doc = "Include sites where the reference is N", required = false)
    boolean includeRefNBases = false;
    /**
     * Use this option to calibrate what bins you want before performing full calculations on your data.
     */
    @Advanced
    @Argument(fullName = "printBinEndpointsAndExit", doc = "Print the bin values and exit immediately", required = false)
    boolean printBinEndpointsAndExit = false;
    /**
     * Sets the low-coverage cutoff for granular binning. All loci with depth < START are counted in the first bin.
     */
    @Advanced
    @Argument(fullName = "start", doc = "Starting (left endpoint) for granular binning", required = false, minValue = 0)
    int start = 1;
    /**
     * Sets the high-coverage cutoff for granular binning. All loci with depth > STOP are counted in the last bin.
     */
    @Advanced
    @Argument(fullName = "stop", doc = "Ending (right endpoint) for granular binning", required = false, minValue = 1)
    int stop = 500;
    /**
     * Sets the number of bins for granular binning
     */
    @Advanced
    @Argument(fullName = "nBins", doc = "Number of bins to use for granular binning", required = false, minValue = 0, minRecommendedValue = 1)
    int nBins = 499;

    /**
     * This option simply disables writing separate files for per-sample summary statistics (total, mean, median, quartile coverage per sample). These statistics are still calculated internally, so enabling this option will not improve runtime.
     */
    @Argument(fullName = "omitPerSampleStats", shortName = "omitSampleSummary", doc = "Do not output the summary files per-sample", required = false)
    boolean omitSampleSummary = false;
    /**
     * By default, coverage is partitioning by sample, but it can be any combination of sample, readgroup and/or library.
     */
    @Argument(fullName = "partitionType", shortName = "pt", doc = "Partition type for depth of coverage", required = false)
    Set<DoCOutputType.Partition> partitionTypes = EnumSet.of(DoCOutputType.Partition.sample);

    /**
     * Consider a spanning deletion as contributing to coverage. Also enables deletion counts in per-base output.
     */
    @Advanced
    @Argument(fullName = "includeDeletions", shortName = "dels", doc = "Include information on deletions", required = false)
    boolean includeDeletions = false;

    @Advanced
    @Argument(fullName = "ignoreDeletionSites", doc = "Ignore sites consisting only of deletions", required = false)
    boolean ignoreDeletionSites = false;
    
    /**
     * For summary file outputs, report the percentage of bases covered to an amount equal to or greater than this number  (e.g. % bases >= CT for each sample). Defaults to 15; can take multiple arguments.
     */
    @Advanced
    @Argument(fullName = "summaryCoverageThreshold", shortName = "ct", doc = "Coverage threshold (in percent) for summarizing statistics", required = false)
    int[] coverageThresholds = {15};

    String[] OUTPUT_FORMATS = {"table","rtable","csv"};
    String separator = "\t";
    Map<DoCOutputType.Partition,List<String>> orderCheck = new HashMap<DoCOutputType.Partition,List<String>>();

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public boolean includeReadsWithDeletionAtLoci() { return includeDeletions && ! ignoreDeletionSites; }

    public static String incompatibleArgsMsg() { return incompatibleArgsMsg; }

    public void initialize() {

        if ( omitIntervals && refSeqGeneList != null ){
            logger.warn(incompatibleArgsMsg);
        }

        if ( printBinEndpointsAndExit ) {
            int[] endpoints = DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins);
            System.out.print("[ ");
            for ( int e : endpoints ) {
                System.out.print(e+" ");
            }
            System.out.println("]");
            System.exit(0);
        }

        // Check the output format
        boolean goodOutputFormat = false;
        for ( String f : OUTPUT_FORMATS ) {
            goodOutputFormat = goodOutputFormat || f.equals(outputFormat);
        }

        if ( ! goodOutputFormat ) {
            throw new IllegalArgumentException("Improper output format. Can be one of table,rtable,csv. Was "+outputFormat);
        }

        if ( outputFormat.equals("csv") ) {
            separator = ",";
        }

        if ( ! omitDepthOutput ) { // print header
            PrintStream out = getCorrectStream(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
            out.printf("%s\t%s","Locus","Total_Depth");
            for (DoCOutputType.Partition type : partitionTypes ) {
                out.printf("\t%s_%s","Average_Depth",type.toString());
            }
            
            // get all the samples
            HashSet<String> allSamples = getSamplesFromToolKit(partitionTypes);
            ArrayList<String> allSampleList = new ArrayList<String>(allSamples.size());
            for ( String s : allSamples ) {
                allSampleList.add(s);
            }
            Collections.sort(allSampleList);

            for ( String s : allSampleList) {
                out.printf("\t%s_%s","Depth_for",s);
                if ( printBaseCounts ) {
                    out.printf("\t%s_%s",s,"base_counts");
                }
            }

            out.printf("%n");

        } else {
            logger.info("Per-Locus Depth of Coverage output was omitted");
        }

        for (DoCOutputType.Partition type : partitionTypes ) {
            orderCheck.put(type,new ArrayList<String>());
            for ( String id : getSamplesFromToolKit(type) ) {
                orderCheck.get(type).add(id);
            }
            Collections.sort(orderCheck.get(type));
        }
    }

    private HashSet<String> getSamplesFromToolKit( Collection<DoCOutputType.Partition> types ) {
        HashSet<String> partitions = new HashSet<String>(); // since the DOCS object uses a HashMap, this will be in the same order
        for (DoCOutputType.Partition t : types ) {
            partitions.addAll(getSamplesFromToolKit(t));
        }

        return partitions;
    }

    private HashSet<String> getSamplesFromToolKit(DoCOutputType.Partition type) {
        HashSet<String> partition = new HashSet<String>();
        if ( type == DoCOutputType.Partition.sample ) {
            partition.addAll(ReadUtils.getSAMFileSamples(getToolkit().getSAMFileHeader()));
        } else if ( type == DoCOutputType.Partition.readgroup ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partition.add(rg.getSample()+"_rg_"+rg.getReadGroupId());
            }
        } else if ( type == DoCOutputType.Partition.library ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partition.add(rg.getLibrary());
            }
        } else if ( type == DoCOutputType.Partition.center ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partition.add(rg.getSequencingCenter());
            }
        } else if ( type == DoCOutputType.Partition.platform ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partition.add(rg.getPlatform());
            }
        } else if ( type == DoCOutputType.Partition.sample_by_center ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partition.add(String.format("%s_cn_%s",rg.getSample(),rg.getSequencingCenter()));
            }
        } else if ( type == DoCOutputType.Partition.sample_by_platform ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partition.add(String.format("%s_pl_%s",rg.getSample(),rg.getPlatform()));
            }
        } else if ( type == DoCOutputType.Partition.sample_by_platform_by_center ) {
            for ( SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
                partition.add(String.format("%s_pl_%s_cn_%s",rg.getSample(),rg.getPlatform(),rg.getSequencingCenter()));
            }
        } else {
            throw new ReviewedGATKException("Invalid aggregation type sent to getSamplesFromToolKit");
        }

        return partition;
    }

    public boolean isReduceByInterval() {
        return ( ! omitIntervals );
    }

    public CoveragePartitioner reduceInit() {
        CoveragePartitioner aggro = new CoveragePartitioner(partitionTypes,start,stop,nBins);
        for (DoCOutputType.Partition t : partitionTypes ) {
            aggro.addIdentifiers(t,getSamplesFromToolKit(t));
        }
        aggro.initialize(includeDeletions,omitLocusTable);
        checkOrder(aggro);
        return aggro;
    }

    public Map<DoCOutputType.Partition,Map<String,int[]>> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (includeRefNBases || BaseUtils.isRegularBase(ref.getBase())) {
            if ( ! omitDepthOutput ) {
                getCorrectStream(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary).printf("%s",ref.getLocus()); // yes: print locus in map, and the rest of the info in reduce (for eventual cumulatives)
                //System.out.printf("\t[log]\t%s",ref.getLocus());
            }

            return CoverageUtils.getBaseCountsByPartition(context,minMappingQuality,maxMappingQuality,minBaseQuality,maxBaseQuality,countType,partitionTypes);
        } else {
            return null;
        }
    }

    public CoveragePartitioner reduce(Map<DoCOutputType.Partition,Map<String,int[]>> thisMap, CoveragePartitioner prevReduce) {
        if ( thisMap != null ) { // skip sites we didn't want to include in the calculation (ref Ns)
            if ( ! omitDepthOutput ) {
                //checkOrder(prevReduce); // tests prevReduce.getIdentifiersByType().get(t) against the initialized header order
                printDepths(getCorrectStream(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary),thisMap,prevReduce.getIdentifiersByType());
                // this is an additional iteration through thisMap, plus dealing with IO, so should be much slower without
                // turning on omit
            }

            prevReduce.update(thisMap); // note that in "useBoth" cases, this method alters the thisMap object
        }

        return prevReduce;
    }

    public CoveragePartitioner treeReduce(CoveragePartitioner left, CoveragePartitioner right) {
        left.merge(right);
        return left;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // INTERVAL ON TRAVERSAL DONE
    ////////////////////////////////////////////////////////////////////////////////////

    public void onTraversalDone( List<Pair<GenomeLoc, CoveragePartitioner>> statsByInterval ) {
        if ( refSeqGeneList != null && partitionTypes.contains(DoCOutputType.Partition.sample) ) {
            printGeneStats(statsByInterval);
        }

        if ( statsByInterval.size() > 0 ) {
            for(DoCOutputType.Partition partition: partitionTypes) {
                if ( checkType(statsByInterval.get(0).getSecond().getCoverageByAggregationType(partition) ,partition) ) {
                    printIntervalStats(statsByInterval,
                            getCorrectStream(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary),
                            getCorrectStream(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.statistics),
                            partition);
                } else {
                    throw new ReviewedGATKException("Partition type "+partition.toString()+" had no entries. Please check that your .bam header has all appropriate partition types.");
                }
            }
        } else {
            throw new UserException.CommandLineException("Cannot reduce by interval without a list of intervals. Please provide an interval list using the -L argument.");
        }

        onTraversalDone(mergeAll(statsByInterval));

    }

    public CoveragePartitioner mergeAll(List<Pair<GenomeLoc, CoveragePartitioner>> stats) {
        CoveragePartitioner first = stats.remove(0).second;
        for ( Pair<GenomeLoc, CoveragePartitioner> iStat : stats ) {
            treeReduce(first,iStat.second);
        }

        return first;
    }

    private DepthOfCoverageStats printIntervalStats(List<Pair<GenomeLoc, CoveragePartitioner>> statsByInterval, PrintStream summaryOut, PrintStream statsOut, DoCOutputType.Partition type) {
        Pair<GenomeLoc, CoveragePartitioner> firstPair = statsByInterval.get(0);
        CoveragePartitioner firstAggregator = firstPair.second;
        DepthOfCoverageStats firstStats = firstAggregator.getCoverageByAggregationType(type);

        StringBuilder summaryHeader = new StringBuilder();
        summaryHeader.append("Target");
        summaryHeader.append(separator);
        summaryHeader.append("total_coverage");
        summaryHeader.append(separator);
        summaryHeader.append("average_coverage");

        for ( String s : firstStats.getAllSamples() ) {
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_total_cvg");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_mean_cvg");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q1");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_median");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q3");
            for ( int thresh : coverageThresholds ) {
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_%_above_");
                summaryHeader.append(thresh);
            }
        }

        summaryOut.printf("%s%n",summaryHeader);

        int[][] nTargetsByAvgCvgBySample = new int[firstStats.getHistograms().size()][firstStats.getEndpoints().length+1];

        for ( Pair<GenomeLoc, CoveragePartitioner> targetAggregator : statsByInterval ) {

            Pair<GenomeLoc,DepthOfCoverageStats> targetStats = new Pair<GenomeLoc,DepthOfCoverageStats>(
                    targetAggregator.first, targetAggregator.second.getCoverageByAggregationType(type));
            printTargetSummary(summaryOut,targetStats);
            updateTargetTable(nTargetsByAvgCvgBySample,targetStats.second);
        }

        printIntervalTable(statsOut,nTargetsByAvgCvgBySample,firstStats.getEndpoints());

        return firstStats;
    }

    private void printGeneStats(List<Pair<GenomeLoc, CoveragePartitioner>> statsByTarget) {
        logger.debug("statsByTarget size is "+Integer.toString(statsByTarget.size()));
        logger.debug("Initializing refseq...");
        LocationAwareSeekableRODIterator refseqIterator = initializeRefSeq();
        logger.debug("Refseq init done.");
        List<Pair<String,DepthOfCoverageStats>> statsByGene = new ArrayList<Pair<String,DepthOfCoverageStats>>();// maintains order
        Map<String,DepthOfCoverageStats> geneNamesToStats = new HashMap<String,DepthOfCoverageStats>(); // allows indirect updating of objects in list

        for ( Pair<GenomeLoc, CoveragePartitioner> targetStats : statsByTarget ) {
            List<String> genes = getGeneNames(targetStats.first,refseqIterator);
            for (String gene : genes) {
            	if ( geneNamesToStats.keySet().contains(gene) ) {
            		logger.debug("Merging "+geneNamesToStats.get(gene).toString()+" and "+targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample).toString());
            		geneNamesToStats.get(gene).merge(targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample));
            	} else {
            		DepthOfCoverageStats merger = new DepthOfCoverageStats(targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample));
            		geneNamesToStats.put(gene,merger);
            		statsByGene.add(new Pair<String,DepthOfCoverageStats>(gene,merger));
            	}
            }
        }

        PrintStream geneSummaryOut = getCorrectStream(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
        StringBuilder summaryHeader = new StringBuilder();
        summaryHeader.append("Gene");
        summaryHeader.append(separator);
        summaryHeader.append("total_coverage");
        summaryHeader.append(separator);
        summaryHeader.append("average_coverage");

        for ( String s : statsByTarget.get(0).second.getCoverageByAggregationType(DoCOutputType.Partition.sample).getAllSamples() ) {
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_total_cvg");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_mean_cvg");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q1");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_median");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q3");
            for ( int thresh : coverageThresholds ) {
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_%_above_");
                summaryHeader.append(thresh);
            }
        }

        geneSummaryOut.printf("%s%n",summaryHeader);

        for ( Pair<String,DepthOfCoverageStats> geneStats : statsByGene ) {
            printTargetSummary(geneSummaryOut,geneStats);
        }
    }

    //blatantly stolen from Andrew Kernytsky
    // edited by Pawel Sztromwasser to support overlap with multiple exons/genes
    private List<String> getGeneNames(GenomeLoc target, LocationAwareSeekableRODIterator refseqIterator) {
        logger.debug("Examining "+target.toString());
        
        List<String> unknown = Arrays.asList("UNKNOWN");
        
        if (refseqIterator == null) { return unknown; }

        RODRecordList annotationList = refseqIterator.seekForward(target);
        logger.debug("Annotation list is " + (annotationList == null ? "null" : annotationList.getName()));
        if (annotationList == null) { return unknown; }

        List<String> geneNames = new ArrayList<String>();
        for(GATKFeature rec : annotationList) {
            if ( ((RefSeqFeature)rec.getUnderlyingObject()).overlapsExonP(target) ) {
                logger.debug("We do overlap "+ rec.getUnderlyingObject().toString());
                geneNames.add(((RefSeqFeature)rec.getUnderlyingObject()).getGeneName());
            } else {
            	logger.debug("No overlap");
            }
        }

        if (geneNames.isEmpty()) { geneNames = unknown; }
        
        return geneNames;
        
    }

    private LocationAwareSeekableRODIterator initializeRefSeq() {
        RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                                                      getToolkit().getGenomeLocParser(),
                                                      getToolkit().getArguments().unsafe,
                                                      getToolkit().getArguments().disableAutoIndexCreationAndLockingWhenReadingRods,
                                                      null);
        RMDTrack refseq = builder.createInstanceOfTrack(RefSeqCodec.class,refSeqGeneList);
        return new SeekableRODIterator(refseq.getHeader(),refseq.getSequenceDictionary(),getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                getToolkit().getGenomeLocParser(),refseq.getIterator());
    }

    private void printTargetSummary(PrintStream output, Pair<?,DepthOfCoverageStats> intervalStats) {
        DepthOfCoverageStats stats = intervalStats.second;
        int[] bins = stats.getEndpoints();

        StringBuilder targetSummary = new StringBuilder();
        targetSummary.append(intervalStats.first.toString());
        targetSummary.append(separator);
        targetSummary.append(stats.getTotalCoverage());
        targetSummary.append(separator);
        targetSummary.append(String.format("%.2f",stats.getTotalMeanCoverage()));

        for ( String s : stats.getAllSamples() ) {
            targetSummary.append(separator);
            targetSummary.append(stats.getTotals().get(s));
            targetSummary.append(separator);
            targetSummary.append(String.format("%.2f", stats.getMeans().get(s)));
            targetSummary.append(separator);
            int median = getQuantile(stats.getHistograms().get(s),0.5);
            int q1 = getQuantile(stats.getHistograms().get(s),0.25);
            int q3 = getQuantile(stats.getHistograms().get(s),0.75);
            targetSummary.append(formatBin(bins,q1));
            targetSummary.append(separator);
            targetSummary.append(formatBin(bins,median));
            targetSummary.append(separator);
            targetSummary.append(formatBin(bins,q3));
            for ( int thresh : coverageThresholds ) {
                targetSummary.append(String.format("%s%.1f",separator,getPctBasesAbove(stats.getHistograms().get(s),stats.value2bin(thresh))));
            }

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
        String colHeader = outputFormat.equals("rtable") ? "" : "Number_of_sources";
        output.printf(colHeader + separator+"depth>=%d",0);
        for ( int col = 0; col < intervalTable[0].length-1; col ++ ) {
            output.printf(separator+"depth>=%d",cutoffs[col]);
        }

        output.printf(String.format("%n"));
        for ( int row = 0; row < intervalTable.length; row ++ ) {
            output.printf("At_least_%d_samples",row+1);
            for ( int col = 0; col < intervalTable[0].length; col++ ) {
                output.printf(separator+"%d",intervalTable[row][col]);
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

    public void onTraversalDone(CoveragePartitioner coverageProfiles) {
        ///////////////////
        // OPTIONAL OUTPUTS
        //////////////////

        if ( ! omitSampleSummary ) {
            logger.info("Printing summary info");
            for (DoCOutputType.Partition type : partitionTypes ) {
                outputSummaryFiles(coverageProfiles,type);
            }
        }

        if ( ! omitLocusTable ) {
            logger.info("Printing locus summary");
            for (DoCOutputType.Partition type : partitionTypes ) {
                outputLocusFiles(coverageProfiles,type);
            }
        }
    }

    private void outputLocusFiles(CoveragePartitioner coverageProfiles, DoCOutputType.Partition type ) {
        printPerLocus(getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_counts),
                getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_proportions),
                coverageProfiles.getCoverageByAggregationType(type),type);
    }

    private void outputSummaryFiles(CoveragePartitioner coverageProfiles, DoCOutputType.Partition type ) {
        printPerSample(getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.statistics),coverageProfiles.getCoverageByAggregationType(type));
        printSummary(getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.summary),coverageProfiles.getCoverageByAggregationType(type));
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // HELPER OUTPUT METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    private void printPerSample(PrintStream output,DepthOfCoverageStats stats) {
        int[] leftEnds = stats.getEndpoints();

        StringBuilder hBuilder = new StringBuilder();
        if ( ! outputFormat.equals("rTable")) {
            hBuilder.append("Source_of_reads");
        }
        hBuilder.append(separator);
        hBuilder.append(String.format("from_0_to_%d)%s",leftEnds[0],separator));
        for ( int i = 1; i < leftEnds.length; i++ )
            hBuilder.append(String.format("from_%d_to_%d)%s",leftEnds[i-1],leftEnds[i],separator));
        hBuilder.append(String.format("from_%d_to_inf%n",leftEnds[leftEnds.length-1]));
        output.print(hBuilder.toString());
        Map<String,long[]> histograms = stats.getHistograms();

        for ( Map.Entry<String, long[]> p : histograms.entrySet() ) {
            StringBuilder sBuilder = new StringBuilder();
            sBuilder.append(String.format("sample_%s",p.getKey()));
            for ( long count : p.getValue() ) {
                sBuilder.append(String.format("%s%d",separator,count));
            }
            sBuilder.append(String.format("%n"));
            output.print(sBuilder.toString());
        }
    }

    private void printPerLocus(PrintStream output, PrintStream coverageOut, DepthOfCoverageStats stats, DoCOutputType.Partition partitionType) {
        int[] endpoints = stats.getEndpoints();
        int samples = stats.getHistograms().size();

        long[][] baseCoverageCumDist = stats.getLocusCounts();

        // rows - # of samples
        // columns - depth of coverage

        boolean printSampleColumnHeader = outputFormat.equals("csv") || outputFormat.equals("table");

        StringBuilder header = new StringBuilder();
        if ( printSampleColumnHeader ) {
            // mhanna 22 Aug 2010 - Deliberately force this header replacement to make sure integration tests pass.
            // TODO: Update integration tests and get rid of this.
            header.append(partitionType == DoCOutputType.Partition.readgroup ? "read_group" : partitionType.toString());
        }
        header.append(String.format("%sgte_0",separator));
        for ( int d : endpoints ) {
            header.append(String.format("%sgte_%d",separator,d));
        }
        header.append(String.format("%n"));

        output.print(header);
        coverageOut.print(header);

        for ( int row = 0; row < samples; row ++ ) {
            output.printf("%s_%d","NSamples",row+1);
            for ( int depthBin = 0; depthBin < baseCoverageCumDist[0].length; depthBin ++ ) {
                output.printf("%s%d",separator,baseCoverageCumDist[row][depthBin]);
            }
            output.printf("%n");
        }

        for ( String sample : stats.getAllSamples() ) {
            coverageOut.printf("%s",sample);
            double[] coverageDistribution = stats.getCoverageProportions(sample);
            for ( int bin = 0; bin < coverageDistribution.length; bin ++ ) {
                coverageOut.printf("%s%.2f",separator,coverageDistribution[bin]);
            }
            coverageOut.printf("%n");
        }
    }

    private PrintStream getCorrectStream(DoCOutputType.Partition partition, DoCOutputType.Aggregation aggregation, DoCOutputType.FileType fileType) {
        DoCOutputType outputType = new DoCOutputType(partition,aggregation,fileType);
        if(!out.containsKey(outputType))
            throw new UserException.CommandLineException(String.format("Unable to find appropriate stream for partition = %s, aggregation = %s, file type = %s",partition,aggregation,fileType));
        return out.get(outputType);
    }

    private void printSummary(PrintStream output, DepthOfCoverageStats stats) {
        if ( ! outputFormat.equals("csv") ) {
            output.printf("%s\t%s\t%s\t%s\t%s\t%s","sample_id","total","mean","granular_third_quartile","granular_median","granular_first_quartile");
        } else {
            output.printf("%s,%s,%s,%s,%s,%s","sample_id","total","mean","granular_third_quartile","granular_median","granular_first_quartile");
        }

        for ( int thresh : coverageThresholds ) {
            output.printf("%s%s%d",separator,"%_bases_above_",thresh);
        }

        output.printf("%n");

        Map<String,long[]> histograms = stats.getHistograms();
        Map<String,Double> means = stats.getMeans();
        Map<String,Long> totals = stats.getTotals();
        int[] leftEnds = stats.getEndpoints();

        for ( Map.Entry<String, long[]> p : histograms.entrySet() ) {
            String s = p.getKey();
            long[] histogram = p.getValue();
            int median = getQuantile(histogram,0.5);
            int q1 = getQuantile(histogram,0.25);
            int q3 = getQuantile(histogram,0.75);
            // if any of these are larger than the higest bin, put the median as in the largest bin
            median =  median == histogram.length-1 ? histogram.length-2 : median;
            q1 = q1 == histogram.length-1 ? histogram.length-2 : q1;
            q3 = q3 == histogram.length-1 ? histogram.length-2 : q3;
            if ( ! outputFormat.equals("csv") ) {
                output.printf("%s\t%d\t%.2f\t%d\t%d\t%d",s,totals.get(s),means.get(s),leftEnds[q3],leftEnds[median],leftEnds[q1]);
            } else {
                output.printf("%s,%d,%.2f,%d,%d,%d",s,totals.get(s),means.get(s),leftEnds[q3],leftEnds[median],leftEnds[q1]);
            }

            for ( int thresh : coverageThresholds ) {
                output.printf("%s%.1f",separator,getPctBasesAbove(histogram,stats.value2bin(thresh)));
            }

            output.printf("%n");
        }

        if ( ! outputFormat.equals("csv") ) {
            output.printf("%s\t%d\t%.2f\t%s\t%s\t%s%n","Total",stats.getTotalCoverage(),stats.getTotalMeanCoverage(),"N/A","N/A","N/A");
        } else {
            output.printf("%s,%d,%.2f,%s,%s,%s%n","Total",stats.getTotalCoverage(),stats.getTotalMeanCoverage(),"N/A","N/A","N/A");
        }
    }

    private int getQuantile(long[] histogram, double prop) {
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

        return bin == -1 ? 0 : bin;
    }

    private double getPctBasesAbove(long[] histogram, int bin) {
        long below = 0l;
        long above = 0l;
        for ( int index = 0; index < histogram.length; index++) {
            if ( index < bin ) {
                below+=histogram[index];
            } else {
                above+=histogram[index];
            }
        }

        return 100*( (double) above )/( above + below );
    }

    private void printDepths(PrintStream stream, Map<DoCOutputType.Partition,Map<String,int[]>> countsBySampleByType, Map<DoCOutputType.Partition,List<String>> identifiersByType) {
        // get the depths per sample and build up the output string while tabulating total and average coverage
        StringBuilder perSampleOutput = new StringBuilder();
        int tDepth = 0;
        boolean depthCounted = false;
        for (DoCOutputType.Partition type : partitionTypes ) {
            Map<String,int[]> countsByID = countsBySampleByType.get(type);
            for ( String s : identifiersByType.get(type) ) {
                perSampleOutput.append(separator);
                long dp = (countsByID != null && countsByID.keySet().contains(s)) ? sumArray(countsByID.get(s)) : 0 ;
                perSampleOutput.append(dp);
                if ( printBaseCounts ) {
                    perSampleOutput.append(separator);
                    perSampleOutput.append(baseCounts(countsByID != null ? countsByID.get(s) : null ));
                }
                if ( ! depthCounted ) {
                    tDepth += dp;
                }
            }
            depthCounted = true; // only sum the total depth once
        }

        // remember -- genome locus was printed in map()
        stream.printf("%s%d",separator,tDepth);
        for (DoCOutputType.Partition type : partitionTypes ) {
            stream.printf("%s%.2f",separator, ( (double) tDepth / identifiersByType.get(type).size() ) );
        }
        stream.printf("%s%n",perSampleOutput);
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
        for ( byte b : BaseUtils.EXTENDED_BASES ) {
            nbases++;
            if ( includeDeletions || b != BaseUtils.Base.D.base ) {
                s.append((char)b);
                s.append(":");
                s.append(counts[BaseUtils.extendedBaseToBaseIndex(b)]);
                if ( nbases < 6 ) {
                    s.append(" ");
                }
            }
        }

        return s.toString();
    }

    private void checkOrder(CoveragePartitioner ag) {
        // make sure the ordering stored at initialize() is propagated along reduce
        for (DoCOutputType.Partition t : partitionTypes ) {
            List<String> order = orderCheck.get(t);
            List<String> namesInAg = ag.getIdentifiersByType().get(t);

            // todo -- chris check me
            Set<String> namesInDOCS = ag.getCoverageByAggregationType(t).getAllSamples();
            int index = 0;
            for ( String s : namesInAg ) {
                if ( ! s.equalsIgnoreCase(order.get(index)) ) {
                    throw new ReviewedGATKException("IDs are out of order for type "+t+"! Aggregator has different ordering");
                }
                index++;
            }
        }
    }

    public boolean checkType(DepthOfCoverageStats stats, DoCOutputType.Partition type ) {
        if ( stats.getHistograms().size() < 1 ) {
            logger.warn("The histogram per partition type "+type.toString()+" was empty\n"+
                    "Do your read groups have this type? (Check your .bam header).");
            return false;
        } else {
            return true;
        }
    }

}

class DoCOutputMultiplexer implements Multiplexer<DoCOutputType> {
    private final Set<DoCOutputType.Partition> partitions;
    private final File refSeqGeneList;
    private final boolean omitDepthOutput;
    private final boolean omitIntervals;
    private final boolean omitSampleSummary;
    private final boolean omitLocusTable;

    /**
     * Create a new multiplexer type using the values of all variable fields.
     * @param partitions
     * @param refSeqGeneList
     * @param omitDepthOutput
     * @param omitIntervals
     * @param omitSampleSummary
     * @param omitLocusTable
     */
    public DoCOutputMultiplexer(final Set<DoCOutputType.Partition> partitions,
                                final File refSeqGeneList,
                                final boolean omitDepthOutput,
                                final boolean omitIntervals,
                                final boolean omitSampleSummary,
                                final boolean omitLocusTable) {
        this.partitions = partitions;
        this.refSeqGeneList = refSeqGeneList;
        this.omitDepthOutput = omitDepthOutput;
        this.omitIntervals = omitIntervals;
        this.omitSampleSummary = omitSampleSummary;
        this.omitLocusTable = omitLocusTable;
    }

    public Collection<DoCOutputType> multiplex() {
        List<DoCOutputType> outputs = new ArrayList<DoCOutputType>();
        if(!omitDepthOutput) outputs.add(new DoCOutputType(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary));

        if(!omitIntervals) {
            for(DoCOutputType.Partition partition: partitions) {
                outputs.add(new DoCOutputType(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary));
                outputs.add(new DoCOutputType(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.statistics));
            }
        }

        if(refSeqGeneList != null && partitions.contains(DoCOutputType.Partition.sample)) {
            DoCOutputType geneSummaryOut = new DoCOutputType(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
            outputs.add(geneSummaryOut);
        }

        if(!omitSampleSummary) {
            for(DoCOutputType.Partition partition: partitions) {
                outputs.add(new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.summary));
                outputs.add(new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.statistics));
            }
        }

        if(!omitLocusTable) {
            for(DoCOutputType.Partition partition: partitions) {
                outputs.add(new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_counts));
                outputs.add(new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_proportions));
            }
        }

        return outputs;
    }

    public String transformArgument(final DoCOutputType outputType, final String argument) {
        return outputType.getFileName(argument);
    }

}

class CoveragePartitioner {
    private Collection<DoCOutputType.Partition> types;
    private Map<DoCOutputType.Partition,DepthOfCoverageStats> coverageProfiles;
    private Map<DoCOutputType.Partition,List<String>> identifiersByType;
    private Set<String> allIdentifiers;
    public CoveragePartitioner(Collection<DoCOutputType.Partition> typesToUse, int start, int stop, int nBins) {
        coverageProfiles = new TreeMap<DoCOutputType.Partition,DepthOfCoverageStats>();
        identifiersByType = new HashMap<DoCOutputType.Partition,List<String>>();
        types = typesToUse;
        for ( DoCOutputType.Partition type : types ) {
            coverageProfiles.put(type,new DepthOfCoverageStats(DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins)));
            identifiersByType.put(type,new ArrayList<String>());
        }
        allIdentifiers = new HashSet<String>();
    }

    public void merge(CoveragePartitioner otherAggregator) {
        for ( DoCOutputType.Partition type : types ) {
            this.coverageProfiles.get(type).merge(otherAggregator.coverageProfiles.get(type));
        }
    }

    public DepthOfCoverageStats getCoverageByAggregationType(DoCOutputType.Partition t) {
        return coverageProfiles.get(t);
    }

    public void addIdentifiers(DoCOutputType.Partition t, Set<String> ids) {
        for ( String s : ids ) {
            coverageProfiles.get(t).addSample(s);
            identifiersByType.get(t).add(s);
            allIdentifiers.add(s);
        }
        Collections.sort(identifiersByType.get(t));
    }

    public void initialize(boolean useDels, boolean omitLocusTable) {
        for ( DoCOutputType.Partition t : types ) {
            if ( useDels ) {
                coverageProfiles.get(t).initializeDeletions();
            }
            if ( ! omitLocusTable ) {
                coverageProfiles.get(t).initializeLocusCounts();
            }
        }
    }

    public void update(Map<DoCOutputType.Partition,Map<String,int[]>> countsByIdentifierByType) {
        for ( DoCOutputType.Partition t : types ) {
            coverageProfiles.get(t).update(countsByIdentifierByType.get(t));
        }
    }

    public Set<String> getAllIdentifiers() {
        return allIdentifiers;
    }

    public Map<DoCOutputType.Partition,List<String>> getIdentifiersByType() {
        return identifiersByType;
    }
}
