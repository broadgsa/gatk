
package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.IntervalWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import net.sf.samtools.*;

@WalkerName("IntervalCleaner")
public class IntervalCleanerWalker extends IntervalWalker<Integer, Integer> {
    @Argument(fullName="maxReadLength", shortName="maxRead", required=false, defaultValue="-1")
    public int maxReadLength;

    public void initialize() {}

    public Integer map(RefMetaDataTracker tracker, String ref, LocusContext context) {
        out.println(context.getLocation() + " (" + context.getReads().size() + ") => " + ref);
        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        out.println("Saw " + result + " intervals");
    }
}