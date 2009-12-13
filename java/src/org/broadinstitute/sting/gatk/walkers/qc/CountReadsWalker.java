package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;

/**
 * Walks over the input data set, calculating the number of reads seen for diagnostic purposes.
 * Can also count the number of reads matching a given criterion using read filters (see the
 * --read-filter command line argument).  Simplest example of a read-backed analysis.  
 */
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
