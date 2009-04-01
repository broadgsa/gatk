package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class PileupWalker extends LocusWalker<Integer, Integer> {
    public boolean FLAG_UNCOVERED_BASES = true;     // todo: how do I make this a command line argument?
    
    public void initialize() {
    }

    // Do we actually want to operate on the context?
    public boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        String bases = Utils.basePileupAsString(reads, offsets);
        String quals = Utils.qualPileupAsString(reads, offsets);

        if ( bases.equals("") && FLAG_UNCOVERED_BASES ) {
            bases = "*** UNCOVERED SITE ***";
        }

        String rodString = "";
        for ( ReferenceOrderedDatum datum : rodData ) {
            if ( datum != null ) {
                if ( datum instanceof rodDbSNP) {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    rodString += dbsnp.toMediumString();
                }
                else {
                    rodString += datum.toSimpleString();
                }
            }
        }
        if ( rodString != "" )
            rodString = "[ROD: " + rodString + "]";

        //if ( context.getLocation().getStart() % 1 == 0 ) {
        out.printf("%s: %s %s %s %s%n", context.getLocation(), ref, bases, quals, rodString);
        //}

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
