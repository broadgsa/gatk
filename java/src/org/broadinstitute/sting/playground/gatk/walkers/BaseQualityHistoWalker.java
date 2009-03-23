package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.BasicReadWalker;
import org.broadinstitute.sting.gatk.LocusContext;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class BaseQualityHistoWalker extends BasicReadWalker<Integer, Integer> {
    long[] qualCounts = new long[100];

    public String getName() {
        return "Base_Quality_Histogram";
    }

    public void initialize() {
        for ( int i = 0; i < this.qualCounts.length; i++ ) {
            this.qualCounts[i] = 0;
        }
    }

    // Do we actually want to operate on the context?
    public boolean filter(LocusContext context, SAMRecord read) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public Integer map(LocusContext context, SAMRecord read) {
        for ( byte qual : read.getBaseQualities() ) {
            //System.out.println(qual);
            this.qualCounts[qual]++;
        }
        //System.out.println(read.getReadName());
        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone() {
        int lastNonZero = -1;
        for ( int i = this.qualCounts.length-1; i >= 0; i-- ) {
            if ( this.qualCounts[i] > 0 ) {
                lastNonZero = i;
                break;
            }
        }

        for ( int i = 0; i < lastNonZero+1; i++ ) {
            System.out.printf("%3d : %10d%n", i, this.qualCounts[i]);
        }
    }
}
