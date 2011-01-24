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

package org.broadinstitute.sting.oneoffprojects.phasing;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
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
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = @RMD(name = "variant", type = ReferenceOrderedDatum.class))

@ReadFilters({ZeroMappingQualityReadFilter.class})
// Filter out all reads with zero mapping quality

public class CountHetPhasingInIntervalWalker extends RodWalker<Integer, Integer> {
    private LinkedList<String> rodNames = null;

    private GenomeLoc prevInterval = null;

    private MultiSampleIntervalStats intervalStats = null;

    @Output
    protected PrintStream out;

    public void initialize() {
        rodNames = new LinkedList<String>();
        rodNames.add("variant");

        intervalStats = new MultiSampleIntervalStats();
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

        List<GATKFeature> interval = tracker.getGATKFeatureMetaData("intervals", true);
        if (interval.size() != 1) {
            String error = "At " + ref.getLocus() + " : Must provide a track named 'intervals' with exactly ONE interval per locus in -L argument!";
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
            intervalStats.startNewInterval();

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

        public MultiSampleIntervalStats() {
            this.sampleToStat = new HashMap<String, SingleSampleIntervalStats>();
            this.numIntervals = 0;
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
            for (SingleSampleIntervalStats stats : sampleToStat.values())
                stats.finalizeStats();
        }

        public Set<Map.Entry<String, SingleSampleIntervalStats>> entrySet() {
            return sampleToStat.entrySet();
        }

        public void startNewInterval() {
            finalizeStats();
            numIntervals++;
        }

        private class SingleSampleIntervalStats {
            public Map<PhasedHetsStat, Integer> hetStatInIntervalToCount;
            public int firstHetIsPhased;

            private int numHetsInCurrentInterval;
            private int numPhasedInCurrentInterval;

            public SingleSampleIntervalStats() {
                this.hetStatInIntervalToCount = new TreeMap<PhasedHetsStat, Integer>(); // implemented PhasedHetsStat.compareTo()
                this.firstHetIsPhased = 0;

                this.numHetsInCurrentInterval = 0;
                this.numPhasedInCurrentInterval = 0;
            }

            public void updateHetStats(boolean isHet, boolean isPhased) {
                if (isHet) {
                    numHetsInCurrentInterval++;

                    if (isPhased) {
                        numPhasedInCurrentInterval++;

                        if (numHetsInCurrentInterval == 1)
                            firstHetIsPhased++;
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

                numHetsInCurrentInterval = 0;
                numPhasedInCurrentInterval = 0;
            }

            public String toString() {
                StringBuilder sb = new StringBuilder();

                sb.append("# of intervals: " + numIntervals + "\n");
                sb.append("First het is phased: " + firstHetIsPhased + "\n");

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