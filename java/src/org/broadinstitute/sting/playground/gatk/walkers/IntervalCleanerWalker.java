
package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWindowWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.indels.*;
import org.broadinstitute.sting.playground.utils.CountedObject;
import org.broadinstitute.sting.playground.utils.CountedObjectComparatorAdapter;

import net.sf.samtools.*;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;
import java.io.File;

@WalkerName("IntervalCleaner")
public class IntervalCleanerWalker extends LocusWindowWalker<Integer, Integer> {
    @Argument(fullName="maxReadLength", shortName="maxRead", doc="max read length", required=false, defaultValue="-1")
    public int maxReadLength;
    @Argument(fullName="OutputCleaned", shortName="O", required=true, doc="Output file (sam or bam) for improved (realigned) reads")
    public String OUT;

    private SAMFileWriter writer;

    public void initialize() {
        SAMFileHeader header = getToolkit().getSamReader().getFileHeader();
        writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, new File(OUT));
    }

    public Integer map(RefMetaDataTracker tracker, String ref, LocusContext context) {
        List<SAMRecord> reads = context.getReads();
        ArrayList<SAMRecord> goodReads = new ArrayList<SAMRecord>();
        long leftmostIndex = context.getLocation().getStart();
        for ( SAMRecord read : reads ) {
            if ( read.getReadLength() <= maxReadLength )
                goodReads.add(read);
        }
        clean(goodReads, ref, context.getLocation().getStart());

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

    private void clean(List<SAMRecord> reads, String reference, long leftmostIndex) {
        // total mismatches across all reads
        //int totalMismatches = 0;
        //TreeSet< CountedObject<Indel> > all_indels = new TreeSet< CountedObject<Indel> >(
        //        new CountedObjectComparatorAdapter<Indel>(new IntervalComparator()));

        for ( SAMRecord read : reads ) {
            System.out.println(read.getReadString());
            System.out.println(reference.substring(read.getAlignmentStart()-(int)leftmostIndex, read.getAlignmentEnd()-(int)leftmostIndex+1));
            //totalMismatches += AlignmentUtils.numMismatches(read, reference);
            //System.out.println(totalMismatches + "\n");
        }



    }
}