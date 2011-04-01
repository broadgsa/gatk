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

package org.broadinstitute.sting.oneoffprojects.walkers.CNV;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.fasta.FastaSequence;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Walks along reference and prints the reference sequence (as FASTA) for the BED file intervals ("intervals" ROD).
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = ReferenceFASTAforBedIntervalsWalker.INTERVALS_ROD_NAME, type = ReferenceOrderedDatum.class)})

public class ReferenceFASTAforBedIntervalsWalker extends RodWalker<Integer, Integer> {
    @Output
    protected PrintStream out;

    public final static String INTERVALS_ROD_NAME = "intervals";

    private Map<GenomeLoc, FastaSequence> intervalSequences;

    private final static int LINE_WIDTH = 60;

    public void initialize() {
        this.intervalSequences = new HashMap<GenomeLoc, FastaSequence>();
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
     * @return number of interval sequences printed
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc curLoc = ref.getLocus();
        int curPos = curLoc.getStart();
        int entries = 0;

        List<GATKFeature> intervals = tracker.getGATKFeatureMetaData(INTERVALS_ROD_NAME, true);
        for (GATKFeature interval : intervals) {
            GenomeLoc loc = interval.getLocation();
            /* TODO: note that an interval may actually start BEFORE here, but not be covered, but would need to cache the remappings
               of origLoc -> newLoc, and then setName(newLoc.toString()) */

            FastaSequence seq = null;
            if (loc.getStart() == curPos) { // at the start of this interval:
                seq = new FastaSequence(out, LINE_WIDTH, false);
                seq.setName(loc.toString());
                intervalSequences.put(loc, seq);
            }
            else {
                seq = intervalSequences.get(loc);
            }

            seq.append(String.valueOf((char) ref.getBase()));

            if (loc.getStop() == curPos) { // at the end of this interval:
                intervalSequences.remove(loc);
                seq.flush();
                entries++;
            }
        }

        return entries;
    }

    public Integer reduce(Integer add, Integer runningCount) {
        if (add == null)
            add = 0;

        return runningCount + add;
    }

    /**
     * @param result the genes found in each interval.
     */
    public void onTraversalDone(Integer result) {
        result += intervalSequences.size();

        for (Map.Entry<GenomeLoc, FastaSequence> locSeqEntry : intervalSequences.entrySet()) {
            GenomeLoc interval = locSeqEntry.getKey();
            FastaSequence seq = locSeqEntry.getValue();

            int actualStop = interval.getStart() + (int)seq.getCurrentCount() - 1;
            GenomeLoc actualInterval = getToolkit().getGenomeLocParser().setStop(interval, actualStop);
            seq.setName(actualInterval.toString());

            seq.flush();
        }

        System.out.println("Printed out " + result + " sequence entries.");
    }
}

