package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * At each locus in the input data set, prints the reference base, genomic location, and
 * all aligning reads in a compact but human-readable form. 
 */
public class PrintLocusContextWalker extends LocusWalker<AlignmentContext, Integer> implements TreeReducible<Integer> {
    @Output
    private PrintStream out;

    public AlignmentContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        out.printf( "In map: ref = %s, loc = %s %s, reads = %s%n", ref.getBaseAsChar(),
                                                                   context.getLocation(),
                                                                   context.hasExtendedEventPileup() ? "[extended]" : "",
                                                                   Arrays.deepToString( getReadNames(context.getReads()) ) );
        return context;
    }



    public Integer reduceInit() { return 0; }

    public Integer reduce(AlignmentContext context, Integer sum) {
        return sum + 1;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    private String[] getReadNames( List<GATKSAMRecord> reads ) {
        String[] readNames = new String[ reads.size() ];
        for( int i = 0; i < reads.size(); i++ ) {
            readNames[i] = String.format("%nname = %s, start = %d, end = %d", reads.get(i).getReadName(), reads.get(i).getAlignmentStart(), reads.get(i).getAlignmentEnd());
        }
        //Arrays.sort(readNames);
        return readNames;
    }

    @Override
    public boolean generateExtendedEvents() {
        return true;
    }    
}
