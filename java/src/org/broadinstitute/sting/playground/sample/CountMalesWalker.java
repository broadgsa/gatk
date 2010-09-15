package org.broadinstitute.sting.playground.sample;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;

/**
 * Walks over the input data set, calculating the number of reads seen for diagnostic purposes.
 * Can also count the number of reads matching a given criterion using read filters (see the
 * --read-filter command line argument).  Simplest example of a read-backed analysis.
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountMalesWalker extends ReadWalker<Integer, Integer> {
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker tracker) {
        Sample sample = getToolkit().getSampleByRead(read);
        return sample.isMale() ? 1 : 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}                                       