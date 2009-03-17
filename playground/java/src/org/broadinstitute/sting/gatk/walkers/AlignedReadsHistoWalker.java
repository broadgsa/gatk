package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.LocusContext;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class AlignedReadsHistoWalker extends BasicReadWalker<Integer, Integer> {
    long[] alignCounts = new long[51];

    public String getName() {
        return "Aligned_Reads_Histogram";
    }

    public void initialize() {
        for ( int i = 0; i < this.alignCounts.length; i++ ) {
            this.alignCounts[i] = 0;
        }
    }

    public String walkerType() { return "ByRead"; }

    // Do we actually want to operate on the context?
    public boolean filter(LocusContext context, SAMRecord read) {
	    // we only want aligned reads
	    return !read.getReadUnmappedFlag();
    }

    // Map over the org.broadinstitute.sting.atk.LocusContext
    public Integer map(LocusContext context, SAMRecord read) {
        //System.out.println(read.getAttribute("NM"));
        int editDist = Integer.parseInt(read.getAttribute("NM").toString());
        this.alignCounts[editDist]++;
        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone() {
        int curTotal = 0;
        for ( int i = 0; i < this.alignCounts.length; i++ ) {
            curTotal += alignCounts[i];
            System.out.printf("%3d %10d%n", i, curTotal);
        }
    }
}
