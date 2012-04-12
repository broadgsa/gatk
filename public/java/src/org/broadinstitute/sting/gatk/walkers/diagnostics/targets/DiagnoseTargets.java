/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import net.sf.picard.util.PeekableIterator;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

import java.util.*;

/**
 * Analyzes coverage distribution and validates read mates for a given interval and sample.
 * <p/>
 * <p>
 * Used to diagnose regions with bad coverage, mapping, or read mating. Analyzes each sample independently in addition
 * to interval wide analysis.
 * </p>
 * <p/>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * <ul>
 * <li>A reference file</li>
 * <li>one or more input BAMs</li>
 * <li>One or more intervals</li>
 * </ul>
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * A modified VCF detailing each interval by sample
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *              -T DiagnoseTargets \
 *              -R reference.fasta \
 *              -o output.vcf \
 *              -I sample1.bam \
 *              -I sample2.bam \
 *              -I sample3.bam \
 *              -L intervals.interval_list
 *  </pre>
 *
 * @author Mauricio Carneiro, Roger Zurawicki
 * @since 5/8/12
 */
@By(value = DataSource.READS)
@PartitionBy(PartitionType.INTERVAL)
public class DiagnoseTargets extends LocusWalker<Long, Long> {

    @Output(doc = "File to which variants should be written", required = true)
    private VariantContextWriter vcfWriter = null;

    @Argument(fullName = "minimum_base_quality", shortName = "BQ", doc = "The minimum Base Quality that is considered for calls", required = false)
    private int minimumBaseQuality = 20;

    @Argument(fullName = "minimum_mapping_quality", shortName = "MQ", doc = "The minimum read mapping quality considered for calls", required = false)
    private int minimumMappingQuality = 20;

    @Argument(fullName = "minimum_coverage", shortName = "min", doc = "The minimum allowable coverage, used for calling LOW_COVERAGE", required = false)
    private int minimumCoverage = 5;

    @Argument(fullName = "maximum_coverage", shortName = "max", doc = "The maximum allowable coverage, used for calling EXCESSIVE_COVERAGE", required = false)
    private int maximumCoverage = 700;

    @Argument(fullName = "minimum_median_depth", shortName = "med", doc = "The minimum allowable median coverage, used for calling LOW_MEDIAN_DEPTH", required = false)
    private int minMedianDepth = 20;

    @Argument(fullName = "maximum_insert_size", shortName = "ins", doc = "The maximum allowed distance between a read and its mate", required = false)
    private int maxInsertSize = 50;

    @Argument(fullName = "voting_status_threshold", shortName = "stV", doc = "The needed percentage of samples containing a call for the interval to adopt the call ", required = false)
    private double votePercentage = 0.50;

    @Argument(fullName = "low_median_depth_status_threshold", shortName = "stMED", doc = "The percentage of the loci needed for calling LOW_MEDIAN_DEPTH", required = false)
    private double lowMedianDepthPercentage = 0.20;

    @Argument(fullName = "bad_mate_status_threshold", shortName = "stBM", doc = "The percentage of the loci needed for calling BAD_MATE", required = false)
    private double badMateStatusThreshold = 0.50;

    @Argument(fullName = "coverage_status_threshold", shortName = "stC", doc = "The percentage of the loci needed for calling LOW_COVERAGE and COVERAGE_GAPS", required = false)
    private double coverageStatusThreshold = 0.20;

    @Argument(fullName = "excessive_coverage_status_threshold", shortName = "stXC", doc = "The percentage of the loci needed for calling EXCESSIVE_COVERAGE", required = false)
    private double excessiveCoverageThreshold = 0.20;

    @Argument(fullName = "quality_status_threshold", shortName = "stQ", doc = "The percentage of the loci needed for calling POOR_QUALITY", required = false)
    private double qualityStatusThreshold = 0.50;

    @Argument(fullName = "print_debug_log", shortName = "dl", doc = "Used only for debugging the walker. Prints extra info to screen", required = false)
    private boolean debug = false;

    private HashMap<GenomeLoc, IntervalStatistics> intervalMap = null;                                                  // interval => statistics
    private PeekableIterator<GenomeLoc> intervalListIterator;                                                           // an iterator to go over all the intervals provided as we traverse the genome
    private Set<String> samples = null;                                                                                 // all the samples being processed
    private final Allele SYMBOLIC_ALLELE = Allele.create("<DT>", false);                                                // avoid creating the symbolic allele multiple times
    private ThresHolder thresholds = null;

    @Override
    public void initialize() {
        super.initialize();

        if (getToolkit().getIntervals() == null)
            throw new UserException("This tool currently only works if you provide one or more interval");

        thresholds = new ThresHolder(minimumBaseQuality, minimumMappingQuality, minimumCoverage, maximumCoverage, minMedianDepth, maxInsertSize, votePercentage, lowMedianDepthPercentage, badMateStatusThreshold, coverageStatusThreshold, excessiveCoverageThreshold, qualityStatusThreshold);

        intervalMap = new HashMap<GenomeLoc, IntervalStatistics>();
        intervalListIterator = new PeekableIterator<GenomeLoc>(getToolkit().getIntervals().iterator());

        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());                                       // get all of the unique sample names for the VCF Header
        vcfWriter.writeHeader(new VCFHeader(getHeaderInfo(), samples));                                                 // initialize the VCF header
    }

    @Override
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc refLocus = ref.getLocus();

        removePastIntervals(refLocus, ref.getBase());                                                                   // process and remove any intervals in the map that are don't overlap the current locus anymore
        addNewOverlappingIntervals(refLocus);                                                                           // add all new intervals that may overlap this reference locus    

        for (IntervalStatistics intervalStatistics : intervalMap.values())
            intervalStatistics.addLocus(context, ref, thresholds);                                                                       // Add current locus to stats

        return 1L;
    }

    @Override
    public Long reduceInit() {
        return 0L;
    }

    /**
     * Not sure what we are going to do here
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return a long
     */
    @Override
    public Long reduce(Long value, Long sum) {
        return sum + value;
    }

    /**
     * Process all remaining intervals
     *
     * @param result number of loci processed by the walker
     */
    @Override
    public void onTraversalDone(Long result) {
        for (GenomeLoc interval : intervalMap.keySet())
            outputStatsToVCF(intervalMap.get(interval), Allele.create("A", true));
    }

    private GenomeLoc getIntervalMapSpan() {
        GenomeLoc loc = null;
        for (GenomeLoc interval : intervalMap.keySet()) {
            if (loc == null) {
                loc = interval;
            } else
                loc = interval.union(loc);
        }

        return loc;
    }

    /**
     * Removes all intervals that are behind the current reference locus from the intervalMap
     *
     * @param refLocus the current reference locus
     * @param refBase  the reference allele
     */
    private void removePastIntervals(GenomeLoc refLocus, byte refBase) {
        // if all intervals are safe
        if (getIntervalMapSpan() != null && getIntervalMapSpan().isBefore(refLocus)) {
            for (GenomeLoc interval : intervalMap.keySet()) {
                outputStatsToVCF(intervalMap.get(interval), Allele.create(refBase, true));
                intervalMap.remove(interval);
            }
        }

        GenomeLoc interval = intervalListIterator.peek();                                                               // clean up all intervals that we might have skipped because there was no data
        while (interval != null && interval.isBefore(refLocus)) {
            interval = intervalListIterator.next();
            outputStatsToVCF(createIntervalStatistic(interval), Allele.create(refBase, true));
            interval = intervalListIterator.peek();
        }
    }

    /**
     * Adds all intervals that overlap the current reference locus to the intervalMap
     *
     * @param refLocus the current reference locus
     */
    private void addNewOverlappingIntervals(GenomeLoc refLocus) {
        GenomeLoc interval = intervalListIterator.peek();
        while (interval != null && !interval.isPast(refLocus)) {
            intervalMap.put(interval, createIntervalStatistic(interval));
            intervalListIterator.next();                                                                                // discard the interval (we've already added it to the map)
            interval = intervalListIterator.peek();
        }
    }

    /**
     * Takes the interval, finds it in the stash, prints it to the VCF
     *
     * @param stats     The statistics of the interval
     * @param refAllele the reference allele
     */
    private void outputStatsToVCF(IntervalStatistics stats, Allele refAllele) {
        GenomeLoc interval = stats.getInterval();

        List<Allele> alleles = new ArrayList<Allele>();
        Map<String, Object> attributes = new HashMap<String, Object>();
        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();

        alleles.add(refAllele);
        alleles.add(SYMBOLIC_ALLELE);
        VariantContextBuilder vcb = new VariantContextBuilder("DiagnoseTargets", interval.getContig(), interval.getStart(), interval.getStart(), alleles);

        vcb = vcb.log10PError(VariantContext.NO_LOG10_PERROR);                                                          // QUAL field makes no sense in our VCF
        vcb.filters(statusesToStrings(stats.callableStatuses(thresholds)));

        attributes.put(VCFConstants.END_KEY, interval.getStop());
        attributes.put(VCFConstants.DEPTH_KEY, stats.averageCoverage());

        vcb = vcb.attributes(attributes);

        for (String sample : samples) {
            Map<String, Object> infos = new HashMap<String, Object>();
            SampleStatistics sampleStat = stats.getSample(sample);
            infos.put(VCFConstants.DEPTH_KEY, sampleStat.averageCoverage());
            infos.put("Q1", sampleStat.getQuantileDepth(0.25));
            infos.put("MED", sampleStat.getQuantileDepth(0.50));
            infos.put("Q3", sampleStat.getQuantileDepth(0.75));

            Set<String> filters = new HashSet<String>();
            filters.addAll(statusesToStrings(stats.getSample(sample).getCallableStatuses(thresholds)));


            genotypes.add(new Genotype(sample, null, VariantContext.NO_LOG10_PERROR, filters, infos, false));
        }
        vcb = vcb.genotypes(genotypes);

        if (debug) {
            System.out.printf("Output -- Interval: %s, Coverage: %.2f%n", stats.getInterval(), stats.averageCoverage());
        }

        vcfWriter.add(vcb.make());

    }

    /**
     * Gets the header lines for the VCF writer
     *
     * @return A set of VCF header lines
     */
    private static Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

        // INFO fields for overall data
        headerLines.add(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1, VCFHeaderLineType.Integer, "Stop position of the interval"));
        headerLines.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a lci divided by interval size."));
        headerLines.add(new VCFInfoHeaderLine("Diagnose Targets", 0, VCFHeaderLineType.Flag, "DiagnoseTargets mode"));

        // FORMAT fields for each genotype
        // todo -- find the appropriate VCF constants
        headerLines.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a lci divided by interval size."));
        headerLines.add(new VCFFormatHeaderLine("Q1", 1, VCFHeaderLineType.Float, "Lower Quartile of depth distribution."));
        headerLines.add(new VCFFormatHeaderLine("MED", 1, VCFHeaderLineType.Float, "Median of depth distribution."));
        headerLines.add(new VCFFormatHeaderLine("Q3", 1, VCFHeaderLineType.Float, "Upper Quartile of depth Distribution."));


        // FILTER fields
        for (CallableStatus stat : CallableStatus.values())
            headerLines.add(new VCFHeaderLine(stat.name(), stat.description));

        return headerLines;
    }

    /**
     * Function that process a set of statuses into strings
     *
     * @param statuses the set of statuses to be converted
     * @return a matching set of strings
     */
    private Set<String> statusesToStrings(Set<CallableStatus> statuses) {
        Set<String> output = new HashSet<String>(statuses.size());

        for (CallableStatus status : statuses)
            output.add(status.name());

        return output;
    }

    private IntervalStatistics createIntervalStatistic(GenomeLoc interval) {
        return new IntervalStatistics(samples, interval /*, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality*/);
    }
}
