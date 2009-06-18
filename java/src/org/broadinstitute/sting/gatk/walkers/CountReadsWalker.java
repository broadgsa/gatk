package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;

@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountReadsWalker extends ReadWalker<Integer, Integer> {
    public Integer map(char[] ref, SAMRecord read) {
        //System.out.println(read.format());
        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
