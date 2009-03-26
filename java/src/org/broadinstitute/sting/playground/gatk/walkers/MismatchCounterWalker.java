package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.Utils;
import edu.mit.broad.picard.reference.ReferenceSequence;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MismatchCounterWalker extends ReadWalker<Integer, Integer> {
    public String getName() {
        return "CountMismatches";
    }

    public Integer map(LocusContext context, SAMRecord read) {
        int nMismatches = 0;
        
        ReferenceSequence refseq = context.getReferenceContig();

        int start = read.getAlignmentStart()-1;
        int stop  = read.getAlignmentEnd();
        if ( read.getAlignmentBlocks().size() == 1 ) {
            // No indels
            List<Byte> refSeq = Utils.subseq(context.getReferenceContig().getBases(), start, stop);
            List<Byte> readBases = Utils.subseq(read.getReadBases());

            assert(refSeq.size() == readBases.size());

            System.out.printf("start, stop = %d %d%n", start, stop);
            System.out.println(read.format());
            System.out.println(Utils.baseList2string(refSeq));
            System.out.println(Utils.baseList2string(readBases));
        }

        return nMismatches;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
