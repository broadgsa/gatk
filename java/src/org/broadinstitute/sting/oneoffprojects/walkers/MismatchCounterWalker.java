package org.broadinstitute.sting.oneoffprojects.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;

@WalkerName("CountMismatches")
public class MismatchCounterWalker extends ReadWalker<Integer, Integer> {
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        int nMismatches = 0;
        
        int start = read.getAlignmentStart()-1;
        int stop  = read.getAlignmentEnd();
        // sometimes BWA outputs screwy reads
        if ( stop - start > ref.getBases().length )
            return 0;

        if ( read.getAlignmentBlocks().size() == 1 ) {
            // No indels
            List<Byte> refSeq = Utils.subseq(ref.getBases());
            List<Byte> readBases = Utils.subseq(read.getReadBases());

            assert(refSeq.size() == readBases.size());

            out.printf("start, stop = %d %d%n", start, stop);
            out.println(read.format());
            out.println(Utils.baseList2string(refSeq));
            out.println(Utils.baseList2string(readBases));
            for ( int i = 0; i < refSeq.size(); i++) {
                if ( refSeq.get(i) != readBases.get(i) )
                    nMismatches++;
            }
            out.println(nMismatches);
        }

        return nMismatches;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
