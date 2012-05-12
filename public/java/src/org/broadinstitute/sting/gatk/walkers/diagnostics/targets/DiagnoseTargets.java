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
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.codecs.vcf.writer.VCFWriter;
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
    private VCFWriter vcfWriter = null;

    @Argument(fullName = "minimum_base_quality", shortName = "mbq", doc = "", required = false)
    private int minimumBaseQuality = 20;

    @Argument(fullName = "minimum_mapping_quality", shortName = "mmq", doc = "", required = false)
    private int minimumMappingQuality = 20;

    @Argument(fullName = "minimum_coverage", shortName = "mincov", doc = "", required = false)
    private int minimumCoverage = 5;

    @Argument(fullName = "maximum_coverage", shortName = "maxcov", doc = "", required = false)
    private int maximumCoverage = 700;

    private HashMap<GenomeLoc, IntervalStatistics> intervalMap = null;                                                  // interval => statistics
    private PeekableIterator<GenomeLoc> intervalListIterator;                                                           // an iterator to go over all the intervals provided as we traverse the genome
    private Set<String> samples = null;                                                                                 // all the samples being processed

    private final Allele SYMBOLIC_ALLELE = Allele.create("<DT>", false);                                                // avoid creating the symbolic allele multiple times

    @Override
    public void initialize() {
        super.initialize();

        if (intervalTrack == null)
            throw new UserException("This tool currently only works if you provide an interval track");

        intervalMap = new HashMap<GenomeLoc, IntervalStatistics>();
        intervalListIterator = new PeekableIterator<GenomeLoc>(intervalTrack.getIntervals(getToolkit()).listIterator());

        samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());                                       // get all of the unique sample names for the VCF Header
        vcfWriter.writeHeader(new VCFHeader(getHeaderInfo(), samples));                                                 // initialize the VCF header
    }

    @Override
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc refLocus = ref.getLocus();

        removePastIntervals(refLocus, ref.getBase());                                                                   // process and remove any intervals in the map that are don't overlap the current locus anymore
        addNewOverlappingIntervals(refLocus);                                                                           // add all new intervals that may overlap this reference locus    

        for (IntervalStatistics intervalStatistics : intervalMap.values())
            intervalStatistics.addLocus(context);                                                                       // Add current locus to stats

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
            processIntervalStats(intervalMap.get(interval), Allele.create("A"));
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

    /**
     * Removes all intervals that are behind the current reference locus from the intervalMap
     *
     * @param refLocus the current reference locus
     * @param refBase the reference allele                
     */
    private void removePastIntervals(GenomeLoc refLocus, byte refBase) {
        List<GenomeLoc> toRemove = new LinkedList<GenomeLoc>();
        for (GenomeLoc interval : intervalMap.keySet())
            if (interval.isBefore(refLocus)) {
                processIntervalStats(intervalMap.get(interval), Allele.create(refBase, true));
                toRemove.add(interval);
            }

        for (GenomeLoc interval : toRemove)
            intervalMap.remove(interval);

        GenomeLoc interval = intervalListIterator.peek();                                                               // clean up all intervals that we might have skipped because there was no data
        while(interval != null && interval.isBefore(refLocus)) {
            interval = intervalListIterator.next();
            processIntervalStats(createIntervalStatistic(interval), Allele.create(refBase, true));
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
           System.out.println("LOCUS : " + refLocus + " -- " + interval);
            intervalMap.put(interval, createIntervalStatistic(interval));
            intervalListIterator.next();                                                                                // discard the interval (we've already added it to the map)
            interval = intervalListIterator.peek();
        }
    }

    /**
     * Takes the interval, finds it in the stash, prints it to the VCF, and removes it
     *
     * @param stats The statistics of the interval
     * @param refAllele the reference allele
     */
    private void processIntervalStats(IntervalStatistics stats, Allele refAllele) {
        GenomeLoc interval = stats.getInterval();
        
        List<Allele> alleles = new ArrayList<Allele>();
        Map<String, Object> attributes = new HashMap<String, Object>();
        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();

        alleles.add(refAllele);
        alleles.add(SYMBOLIC_ALLELE);
        VariantContextBuilder vcb = new VariantContextBuilder("DiagnoseTargets", interval.getContig(), interval.getStart(), interval.getStart(), alleles);

        vcb = vcb.log10PError(VariantContext.NO_LOG10_PERROR);                                                          // QUAL field makes no sense in our VCF
        vcb.filters(statusesToStrings(stats.callableStatuses()));

        attributes.put(VCFConstants.END_KEY, interval.getStop());
        attributes.put(VCFConstants.DEPTH_KEY, stats.averageCoverage());

        vcb = vcb.attributes(attributes);

        for (String sample : samples) {
            Map<String, Object> infos = new HashMap<String, Object>();
            infos.put(VCFConstants.DEPTH_KEY, stats.getSample(sample).averageCoverage());

            Set<String> filters = new HashSet<String>();
            filters.addAll(statusesToStrings(stats.getSample(sample).getCallableStatuses()));


            genotypes.add(new Genotype(sample, null, VariantContext.NO_LOG10_PERROR, filters, infos, false));
        }
        vcb = vcb.genotypes(genotypes);

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
        headerLines.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a lci divided by interval size."));

        // FILTER fields
        for (CallableStatus stat : CallableStatus.values())
            headerLines.add(new VCFHeaderLine(stat.name(), stat.description));

        return headerLines;
    }


    private static Set<String> statusesToStrings(Set<CallableStatus> statuses) {
        Set<String> output = new HashSet<String>(statuses.size());

        for (CallableStatus status : statuses)
            output.add(status.name());

        return output;
    }

    private IntervalStatistics createIntervalStatistic(GenomeLoc interval) {
        return new IntervalStatistics(samples, interval, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality);
    }
}
