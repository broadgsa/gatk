package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class PileupWalker extends BasicLociWalker<Integer, Integer> {
    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
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

        String rodString = "";
        for ( ReferenceOrderedDatum datum : rodData ) {
            if ( datum != null ) {
                if ( datum instanceof rodDbSNP) {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    //System.out.printf("  DBSNP %s on %s => %s%n", dbsnp.toSimpleString(), dbsnp.strand, Utils.join("/", dbsnp.getAllelesFWD()));
                    rodString += dbsnp.toMediumString();
                }
                else {
                    rodString += datum.toSimpleString();
                }
            }
        }
        if ( rodString != "" )
            rodString = "[ROD: " + rodString + "]";

        if ( context.getLocation().getStart() % 1 == 0 ) {
            System.out.printf("%s: %s %s %s %s%n", context.getLocation(), ref, bases, quals, rodString);
        }

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
}
