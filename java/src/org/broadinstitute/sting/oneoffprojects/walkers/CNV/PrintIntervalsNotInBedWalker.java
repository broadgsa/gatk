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
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.PrintStream;
import java.util.List;

/**
 * Walks along reference and prints intervals of sequence not covered in ANY interval in "intervals" ROD.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = {@RMD(name = PrintIntervalsNotInBedWalker.INTERVALS_ROD_NAME, type = ReferenceOrderedDatum.class)})
@By(DataSource.REFERENCE) // So that we will actually enter loci with no ROD on them

public class PrintIntervalsNotInBedWalker extends RodWalker<Integer, Integer> {
    @Output
    protected PrintStream out;

    public final static String INTERVALS_ROD_NAME = "intervals";

    private GenomeLoc waitingInterval = null;

    public void initialize() {
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
     * @return number of intervals printed.
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc curLoc = ref.getLocus();
        int curPos = curLoc.getStart();
        int printed = 0;

        List<GATKFeature> intervals = tracker.getGATKFeatureMetaData(INTERVALS_ROD_NAME, true);
        if (intervals.isEmpty()) {
            if (waitingInterval != null && curLoc.compareContigs(waitingInterval) == 0 && curPos == waitingInterval.getStop() + 1) {
                waitingInterval = getToolkit().getGenomeLocParser().setStop(waitingInterval, curPos);
            }
            else {
                printed += printWaitingIntervalAsBed();
                waitingInterval = ref.getLocus().clone();
            }
        }
        else {
            printed += printWaitingIntervalAsBed();
        }

        return printed;
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
        result += printWaitingIntervalAsBed();

        System.out.println("Printed out " + result + " intervals.");
    }

    private int printWaitingIntervalAsBed() {
        if (waitingInterval == null)
            return 0;

        out.println(waitingInterval.getContig() + "\t" + (waitingInterval.getStart() - 1) + "\t" + waitingInterval.getStop());
        waitingInterval = null;

        return 1;
    }
}

