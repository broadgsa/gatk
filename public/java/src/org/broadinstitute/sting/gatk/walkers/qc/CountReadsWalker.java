package org.broadinstitute.sting.gatk.walkers.qc;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

/**
 * Walks over the input data set, calculating the number of reads seen for diagnostic purposes.
 * Can also count the number of reads matching a given criterion using read filters (see the
 * --read-filter command line argument).  Simplest example of a read-backed analysis.
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountReadsWalker extends ReadWalker<Integer, Integer> {
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker tracker) {

        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
