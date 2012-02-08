package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.IntervalBinding;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocComparator;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

/**
 * Short one line description of the walker.
 *
 * <p>
 * [Long description of the walker]
 * </p>
 *
 *
 * <h2>Input</h2>
 * <p>
 * [Description of the Input]
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * [Description of the Output]
 * </p>
 *
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
public class DiagnoseTargets extends LocusWalker<Long, Long> {
    @Input(fullName = "interval_track", shortName = "int", doc = "", required = true)
    private IntervalBinding<Feature> intervalTrack = null;

    @Output
    private PrintStream out = System.out;

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

    private TreeSet<GenomeLoc> intervalList = null;                     // The list of intervals of interest (plus expanded intervals if user wants them)
    private HashMap<GenomeLoc, IntervalStatistics> intervalMap = null;  // interval => statistics
    private Iterator<GenomeLoc> intervalListIterator;                   // An iterator to go over all the intervals provided as we traverse the genome
    private GenomeLoc currentInterval = null;                           // The "current" interval loaded and being filled with statistics
    private IntervalStatistics currentIntervalStatistics = null;                 // The "current" interval loaded and being filled with statistics

    private GenomeLocParser parser;                                     // just an object to allow us to create genome locs (for the expanded intervals)

    @Override
    public void initialize() {
        super.initialize();

        if (intervalTrack == null)
            throw new UserException("This tool currently only works if you provide an interval track");

        parser = new GenomeLocParser(getToolkit().getMasterSequenceDictionary());       // Important to initialize the parser before creating the intervals below

        List<GenomeLoc> originalList = intervalTrack.getIntervals(getToolkit());        // The original list of targets provided by the user that will be expanded or not depending on the options provided
        intervalList = new TreeSet<GenomeLoc>(new GenomeLocComparator());
        intervalMap = new HashMap<GenomeLoc, IntervalStatistics>(originalList.size() * 2);
        for (GenomeLoc interval : originalList)
            addAndExpandIntervalToLists(interval);

        intervalListIterator = intervalList.iterator();
    }

    @Override
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc refLocus = ref.getLocus();
        while (currentInterval == null || currentInterval.isBefore(refLocus)) {
            if (!intervalListIterator.hasNext())
                return 0L;

            currentInterval = intervalListIterator.next();
            currentIntervalStatistics = intervalMap.get(currentInterval);
        }

        if (currentInterval.isPast(refLocus))
            return 0L;

        byte[] mappingQualities = context.getBasePileup().getMappingQuals();
        byte[] baseQualities = context.getBasePileup().getQuals();
        int coverage = context.getBasePileup().getBaseAndMappingFilteredPileup(minimumBaseQuality, minimumMappingQuality).depthOfCoverage();
        int rawCoverage = context.size();

        IntervalStatisticLocus locusData = new IntervalStatisticLocus(mappingQualities, baseQualities, coverage, rawCoverage);
        currentIntervalStatistics.addLocus(refLocus, locusData);

        return 1L;
    }

    @Override
    public Long reduceInit() {
        return 0L;
    }

    @Override
    public Long reduce(Long value, Long sum) {
        return sum + value;
    }

    @Override
    public void onTraversalDone(Long result) {
        super.onTraversalDone(result);
        out.println("Interval\tCallStatus\tCOV\tAVG");
        for (GenomeLoc interval : intervalList) {
            IntervalStatistics stats = intervalMap.get(interval);
            out.println(String.format("%s\t%s\t%d\t%f", interval, stats.callableStatus(), stats.totalCoverage(), stats.averageCoverage()));
        }
    }

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

    private void addAndExpandIntervalToLists(GenomeLoc interval) {
        if (expandInterval > 0) {
            GenomeLoc before = createIntervalBefore(interval);
            GenomeLoc after = createIntervalAfter(interval);
            intervalList.add(before);
            intervalList.add(after);
            intervalMap.put(before, new IntervalStatistics(before, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality));
            intervalMap.put(after, new IntervalStatistics(after, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality));
        }
        intervalList.add(interval);
        intervalMap.put(interval, new IntervalStatistics(interval, minimumCoverage, maximumCoverage, minimumMappingQuality, minimumBaseQuality));
    }
}
