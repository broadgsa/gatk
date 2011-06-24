/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.oneoffprojects.walkers.phasing;

import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.PrintStream;
import java.util.*;

/**
 * Walks along all variant ROD loci and verifies the phasing from the reads for user-defined pairs of sites.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = "variant", type = ReferenceOrderedDatum.class), @RMD(name = CountHetPhasingInIntervalWalker.INTERVALS_ROD_NAME, type = ReferenceOrderedDatum.class)})

@ReadFilters({ZeroMappingQualityReadFilter.class})
// Filter out all reads with zero mapping quality

public class CountHetPhasingInIntervalWalker extends RodWalker<Integer, Integer> {
    private LinkedList<String> rodNames = null;

    private GenomeLoc prevInterval = null;

    private MultiSampleIntervalStats intervalStats = null;

    @Output
    protected PrintStream out;

    @Argument(fullName = "perIntervalOut", shortName = "perIntervalOut", doc = "File to which to write per-sample, per-interval phased het statistics", required = false)
    protected PrintStream perIntervalOut = null;

    public final static String INTERVALS_ROD_NAME = "intervals";

    public void initialize() {
        rodNames = new LinkedList<String>();
        rodNames.add("variant");

        intervalStats = new MultiSampleIntervalStats(perIntervalOut);
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        int processed = 1;

        List<GATKFeature> interval = tracker.getGATKFeatureMetaData(INTERVALS_ROD_NAME, true);
        if (interval.size() != 1) {
            String error = "At " + ref.getLocus() + " : Must provide a track named '"+ INTERVALS_ROD_NAME  +"' with exactly ONE interval per locus in -L argument!";
            if (interval.size() < 1)
                throw new UserException(error);
            else // interval.size() > 1
                logger.warn(error);
        }
        // Take the FIRST interval covering this locus, and WARN about multiple intervals (above):
        GenomeLoc curInterval = interval.get(0).getLocation();
        logger.debug("refLocus: " + ref.getLocus() + "\tcurInterval = " + curInterval);

        boolean isNewInterval = (prevInterval == null || !curInterval.equals(prevInterval));
        if (isNewInterval)
            intervalStats.startNewInterval(curInterval);

        boolean requireStartHere = true; // only see each VariantContext once
        boolean takeFirstOnly = false; // take as many entries as the VCF file has
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly)) {
            Map<String, Genotype> sampToGenotypes = vc.getGenotypes();
            for (Map.Entry<String, Genotype> sampEntry : sampToGenotypes.entrySet()) {
                Genotype gt = sampEntry.getValue();
                intervalStats.processHetSiteInInterval(sampEntry.getKey(), gt.isHet(), gt.isPhased());
            }
        }

        prevInterval = curInterval;

        return processed;
    }

    public Integer reduce(Integer addIn, Integer runningCount) {
        if (addIn == null)
            addIn = 0;

        return runningCount + addIn;
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(Integer result) {
        intervalStats.finalizeStats();

        System.out.println("Processed " + result + " sites.");

        for (Map.Entry<String, MultiSampleIntervalStats.SingleSampleIntervalStats> sampleEntry : intervalStats.entrySet()) {
            out.println("Sample:\t" + sampleEntry.getKey());
            out.println(sampleEntry.getValue() + "\n");
        }
    }

    class MultiSampleIntervalStats {
        private Map<String, SingleSampleIntervalStats> sampleToStat;
        protected int numIntervals;

        private PrintStream perIntervalOut;
        private GenomeLoc curInterval;

        public MultiSampleIntervalStats(PrintStream perIntervalOut) {
            this.sampleToStat = new HashMap<String, SingleSampleIntervalStats>();
            this.numIntervals = 0;
            this.perIntervalOut = perIntervalOut;
        }

        public void processHetSiteInInterval(String sample, boolean isHet, boolean isPhased) {
            SingleSampleIntervalStats sampleStats = sampleToStat.get(sample);
            if (sampleStats == null) {
                sampleStats = new SingleSampleIntervalStats();
                sampleToStat.put(sample, sampleStats);
            }

            sampleStats.updateHetStats(isHet, isPhased);
        }

        public void finalizeStats() {
            if (curInterval == null)
                return;

            for (Map.Entry<String, SingleSampleIntervalStats> sampleStatEntry : sampleToStat.entrySet()) {
                SingleSampleIntervalStats stats = sampleStatEntry.getValue();
                if (perIntervalOut != null && stats.numHetsInCurrentInterval > 0) {
                    String sample = sampleStatEntry.getKey();
                    perIntervalOut.println(sample + "\t" + curInterval + "\t" + stats.numPhasedInCurrentInterval + "\t" + stats.numHetsInCurrentInterval + "\t" + stats.firstHetIsPhasedInCurrentInterval);
                }
                stats.finalizeStats(); // now, can reset the counters [after print-out]
            }
        }

        public void startNewInterval(GenomeLoc curInterval) {
            finalizeStats();
            numIntervals++;
            this.curInterval = curInterval;
        }

        public Set<Map.Entry<String, SingleSampleIntervalStats>> entrySet() {
            return sampleToStat.entrySet();
        }

        private class SingleSampleIntervalStats {
            public Map<PhasedHetsStat, Integer> hetStatInIntervalToCount;
            public int firstHetIsPhasedCount;

            private int numHetsInCurrentInterval;
            private int numPhasedInCurrentInterval;
            private boolean firstHetIsPhasedInCurrentInterval;

            public SingleSampleIntervalStats() {
                this.hetStatInIntervalToCount = new TreeMap<PhasedHetsStat, Integer>(); // implemented PhasedHetsStat.compareTo()
                this.firstHetIsPhasedCount = 0;

                resetCurrentIntervalCounters();
            }

            private void resetCurrentIntervalCounters() {
                this.numHetsInCurrentInterval = 0;
                this.numPhasedInCurrentInterval = 0;
                this.firstHetIsPhasedInCurrentInterval = false;
            }

            public void updateHetStats(boolean isHet, boolean isPhased) {
                if (isHet) {
                    numHetsInCurrentInterval++;

                    if (isPhased) {
                        numPhasedInCurrentInterval++;

                        if (numHetsInCurrentInterval == 1)
                            firstHetIsPhasedInCurrentInterval = true;
                    }
                }
            }

            public void finalizeStats() {
                if (numIntervals == 0) // have not yet seen any intervals
                    return;

                PhasedHetsStat hetsAndPhased = new PhasedHetsStat(numHetsInCurrentInterval, numPhasedInCurrentInterval);
                Integer cnt = hetStatInIntervalToCount.get(hetsAndPhased);
                if (cnt == null)
                    cnt = 0;
                hetStatInIntervalToCount.put(hetsAndPhased, cnt + 1);

                if (firstHetIsPhasedInCurrentInterval)
                    firstHetIsPhasedCount++;

                resetCurrentIntervalCounters();
            }

            public String toString() {
                StringBuilder sb = new StringBuilder();

                sb.append("# of intervals: " + numIntervals + "\n");
                sb.append("First het is phased: " + firstHetIsPhasedCount + "\n");

                sb.append("Distribution of number of phased / hets per interval:" + "\n");
                for (Map.Entry<PhasedHetsStat, Integer> hetStatEntry : hetStatInIntervalToCount.entrySet())
                    sb.append(hetStatEntry.getKey() + "\t" + hetStatEntry.getValue() + "\n");

                return sb.toString();
            }
        }
    }

    class PhasedHetsStat implements Comparable<PhasedHetsStat> {
        public int numHets;
        public int numPhased;

        public PhasedHetsStat(int numHets, int numPhased) {
            this.numHets = numHets;
            this.numPhased = numPhased;
        }

        public int compareTo(PhasedHetsStat that) {
            if (this.numHets != that.numHets)
                return this.numHets - that.numHets;

            return this.numPhased - that.numPhased;
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append(numPhased + " / " + numHets);

            return sb.toString();
        }
    }
}