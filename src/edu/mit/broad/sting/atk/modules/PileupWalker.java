package edu.mit.broad.sting.atk.modules;

import edu.mit.broad.sting.atk.LocusWalker;
import edu.mit.broad.sting.atk.LocusIterator;
import edu.mit.broad.sam.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class PileupWalker implements LocusWalker<Integer, Integer> {
    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(char ref, LocusIterator context) {
        return true;    // We are keeping all the reads
    }

    // Map over the edu.mit.broad.sting.atk.LocusContext
    public Integer map(char ref, LocusIterator context) {
        //System.out.printf("Reads %s:%d %d%n", context.getContig(), context.getPosition(), context.getReads().size());
        //for ( SAMRecord read : context.getReads() ) {
        //    System.out.println("  -> " + read.getReadName());
        //}

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        String bases = "";
        String quals = "";
        //String offsetString = "";
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            //if ( offset >= read.getReadString().length() )
            //    System.out.printf("  [%2d] [%s] %s%n", offset, read.format(), read.getReadString());

            bases += read.getReadString().charAt(offset);
            quals += read.getBaseQualityString().charAt(offset);
            //offsetString += i;
            //System.out.printf("  [%2d] [%s] %s%n", offset, read.getReadString().charAt(offset), read.getReadString());
        }

        if ( context.getPosition() % 10 == 0 )
            System.out.printf("%s:%d: %s %s %s%n", context.getContig(), context.getPosition(), ref, bases, quals);

        //for ( int offset : context.getOffsets() ) {
        //    System.out.println("  -> " + read.getReadName());
        //}
        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraveralDone() {
    }
}