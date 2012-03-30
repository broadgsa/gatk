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

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocComparator;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.*;

/**
 * Short one line description of the walker.
 * <p/>
 * <p>
 * [Long description of the walker]
 * </p>
 * <p/>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * [Description of the Input]
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * [Description of the Output]
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T [walker name]
 *  </pre>
 *
 * @author Mauricio Carneiro
 * @since 2/1/12
 */
@By(value = DataSource.READS)
@PartitionBy(PartitionType.INTERVAL)
public class DiagnoseTargets extends LocusWalker<Long, Long> implements AnnotatorCompatibleWalker {
    @Input(fullName = "interval_track", shortName = "int", doc = "", required = true)
    private IntervalBinding<Feature> intervalTrack = null;

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName = "expand_interval", shortName = "exp", doc = "", required = false)
    private int expandInterval = 50;

    @Argument(fullName = "minimum_base_quality", shortName = "mbq", doc = "", required = false)
    private int minimumBaseQuality = 20;

    @Argument(fullName = "minimum_mapping_quality", shortName = "mmq", doc = "", required = false)
    private int minimumMappingQuality = 20;

    @Argument(fullName = "minimum_coverage", shortName = "mincov", doc = "", required = false)
    private int minimumCoverage = 5;

    @Argument(fullName = "maximum_coverage", shortName = "maxcov", doc = "", required = false)
    private int maximumCoverage = 700;

    private TreeSet<GenomeLoc> intervalList = null;                                                                     // The list of intervals of interest (plus expanded intervals if user wants them)
    private HashMap<GenomeLoc, IntervalStatistics> intervalMap = null;                                                  // interval => statistics
    private Iterator<GenomeLoc> intervalListIterator;                                                                   // An iterator to go over all the intervals provided as we traverse the genome
    private GenomeLoc currentInterval = null;                                                                           // The "current" interval loaded
    private IntervalStatistics currentIntervalStatistics = null;                                                        // The "current" interval being filled with statistics
    private Set<String> samples = null;                                                                                 // All the samples being processed
    private GenomeLocParser parser;                                                                                     // just an object to allow us to create genome locs (for the expanded intervals)

    @Override
    public void initialize() {
        super.initialize();

        if (intervalTrack == null)
            throw new UserException("This tool currently only works if you provide an interval track");

        parser = new GenomeLocParser(getToolkit().getMasterSequenceDictionary());                                       // Important to initialize the parser before creating the intervals below

        List<GenomeLoc> originalList = intervalTrack.getIntervals(getToolkit());                                        // The original list of targets provided by the user that will be expanded or not depending on the options provided
        intervalList = new TreeSet<GenomeLoc>(new GenomeLocComparator());
        intervalMap = new HashMap<GenomeLoc, IntervalStatistics>();
        for (GenomeLoc interval : originalList)
            intervalList.add(interval);
        //addAndExpandIntervalToMap(interval);

        intervalListIterator = intervalList.iterator();

        // get all of the unique sample names
        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());

        // initialize the header
        Set<VCFHeaderLine> headerInfo = getHeaderInfo();

        vcfWriter.writeHeader(new VCFHeader(headerInfo, samples));
    }

    /**
     * Gets the header lines for the VCF writer
     *
     * @return A set of VCF header lines
     */
    private Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

        // INFO fields for overall data
        headerLines.add(new VCFInfoHeaderLine("END", 1, VCFHeaderLineType.Integer, "Stop position of the interval"));
        headerLines.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total depth in the site. Sum of the depth of all pools"));
        headerLines.add(new VCFInfoHeaderLine("AD", 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a lci divided by interval size."));
        headerLines.add(new VCFInfoHeaderLine("Diagnose Targets", 0, VCFHeaderLineType.Flag, "DiagnoseTargets mode"));

        // FORMAT fields for each genotype
        headerLines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total depth in the site. Sum of the depth of all pools"));
        headerLines.add(new VCFFormatHeaderLine("AD", 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a lci divided by interval size."));

        // FILTER fields

        for (CallableStatus stat : CallableStatus.values()) {
            headerLines.add(new VCFHeaderLine(stat.name(), stat.description));
        }

        return headerLines;
    }

    @Override
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc refLocus = ref.getLocus();
        while (currentInterval == null || currentInterval.isBefore(refLocus)) {                                         // do this for first time and while currentInterval is behind current locus
            if (!intervalListIterator.hasNext())
                return 0L;

            if (currentInterval != null)
                processIntervalStats(currentInterval, Allele.create(ref.getBase(), true));

            currentInterval = intervalListIterator.next();
            addAndExpandIntervalToMap(currentInterval);
            currentIntervalStatistics = intervalMap.get(currentInterval);
        }

        if (currentInterval.isPast(refLocus))                                                                           // skip if we are behind the current interval
            return 0L;

        currentIntervalStatistics.addLocus(context);                                                                    // Add current locus to stats

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

    @Override
    public void onTraversalDone(Long result) {
        for (GenomeLoc interval : intervalMap.keySet()) 
            processIntervalStats(interval, Allele.create("<DT>", true));
    }

    @Override
    public RodBinding<VariantContext> getSnpEffRodBinding() {return null;}

    @Override
    public RodBinding<VariantContext> getDbsnpRodBinding() {return null;}

    @Override
    public List<RodBinding<VariantContext>> getCompRodBindings() {return null;}

    @Override
    public List<RodBinding<VariantContext>> getResourceRodBindings() {return null;}

    @Override
    public boolean alwaysAppendDbsnpId() {return false;}

    private GenomeLoc createIntervalBefore(GenomeLoc interval) {
        int start = Math.max(interval.getStart() - expandInterval, 0);
        int stop = Math.max(interval.getStart() - 1, 0);
        return parser.createGenomeLoc(interval.getContig(), interval.getContigIndex(), start, stop);
    }

    private GenomeLoc createIntervalAfter(GenomeLoc interval) {
        int contigLimit = getToolkit().getSAMFileHeader().getSequenceDictionary().getSequence(interval.getContigIndex()).getSequenceLength();
        int start = Math.min(interval.getStop() + 1, contigLimit);
        int stop = Math.min(interval.getStop() + expandInterval, contigLimit);
        return parser.createGenomeLoc(interval.getContig(), interval.getContigIndex(), start, stop);
    }

    /**
     * Takes an interval and commits it to memory.
     * It will expand it if so told by the -exp command line argument
     *
     * @param interval The new interval to process
     */
    private void addAndExpandIntervalToMap(GenomeLoc interval) {
        if (expandInterval > 0) {
            GenomeLoc before = createIntervalBefore(interval);
            GenomeLoc after = createIntervalAfter(interval);
            intervalList.add(before);
            intervalList.add(after);
            intervalMap.put(before, new IntervalStatistics(samples, before, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality));
            intervalMap.put(after, new IntervalStatistics(samples, after, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality));
        }
        if (!intervalList.contains(interval))
            intervalList.add(interval);
        intervalMap.put(interval, new IntervalStatistics(samples, interval, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality));
    }

    /**
     * Takes the interval, finds it in the stash, prints it to the VCF, and removes it
     *
     * @param interval The interval in memory that you want to write out and clear
     * @param allele the allele
     */
    private void processIntervalStats(GenomeLoc interval, Allele allele) {
        IntervalStatistics stats = intervalMap.get(interval);

        List<Allele> alleles = new ArrayList<Allele>();
        Map<String, Object> attributes = new HashMap<String, Object>();
        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();

        alleles.add(allele);
        VariantContextBuilder vcb = new VariantContextBuilder("DiagnoseTargets", interval.getContig(), interval.getStart(), interval.getStop(), alleles);

        vcb = vcb.log10PError(VariantContext.NO_LOG10_PERROR);                                                          // QUAL field makes no sense in our VCF
        vcb.filters(statusesToStrings(stats.callableStatuses()));

        attributes.put(VCFConstants.END_KEY, interval.getStop());
        attributes.put(VCFConstants.DEPTH_KEY, stats.totalCoverage());
        attributes.put("AV", stats.averageCoverage());

        vcb = vcb.attributes(attributes);

        for (String sample : samples) {
            Map<String, Object> infos = new HashMap<String, Object>();
            infos.put("DP", stats.getSample(sample).totalCoverage());
            infos.put("AV", stats.getSample(sample).averageCoverage());

            Set<String> filters = new HashSet<String>();
            filters.addAll(statusesToStrings(stats.getSample(sample).getCallableStatuses()));


            genotypes.add(new Genotype(sample, alleles, VariantContext.NO_LOG10_PERROR, filters, infos, false));
        }
        vcb = vcb.genotypes(genotypes);

        vcfWriter.add(vcb.make());

        intervalMap.remove(interval);
    }

    private static Set<String> statusesToStrings(Set<CallableStatus> statuses) {
        Set<String> output = new HashSet<String>(statuses.size());

        for (CallableStatus status : statuses)
            output.add(status.name());

        return output;
    }
}
