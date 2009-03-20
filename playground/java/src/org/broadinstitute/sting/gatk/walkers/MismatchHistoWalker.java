package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.Utils;
import edu.mit.broad.picard.reference.ReferenceSequence;

import java.util.Iterator;
import java.util.List;
import static java.lang.reflect.Array.*;

public class MismatchHistoWalker extends BasicReadWalker<Integer, Integer> {

    protected long[] mismatchCounts = new long[0];
    protected final int MIN_TARGET_EDIT_DISTANCE = 5;
    protected final int MAX_TARGET_EDIT_DISTANCE = 10;

    public String getName() {
        return "Mismatch_Histogram";
    }

    // Do we actually want to operate on the context?
    public boolean filter(LocusContext context, SAMRecord read) {
	    // we only want aligned reads
	    return !read.getReadUnmappedFlag();
    }

    public Integer map(LocusContext context, SAMRecord read) {

        int editDist = Integer.parseInt(read.getAttribute("NM").toString());

        // ignore alignments with indels for now
        if ( read.getAlignmentBlocks().size() == 1 &&
             editDist >= MIN_TARGET_EDIT_DISTANCE &&
             editDist <= MAX_TARGET_EDIT_DISTANCE ) {

            ReferenceSequence refseq = context.getReferenceContig();

            int start = read.getAlignmentStart()-1;
            int stop  = read.getAlignmentEnd();

            List<Byte> refSeq = Utils.subseq(context.getReferenceContig().getBases(), start, stop);
            List<Byte> readBases = Utils.subseq(read.getReadBases());

            assert(refSeq.size() == readBases.size());

            // it's actually faster to reallocate a resized array than to use ArrayLists...
            if ( refSeq.size() > mismatchCounts.length ) {
                int oldLength = mismatchCounts.length;
                mismatchCounts = (long[])resizeArray(mismatchCounts, refSeq.size());
                for ( int i = oldLength; i < refSeq.size(); i++ )
                    mismatchCounts[i] = 0;
            }

            String refStr = Utils.baseList2string(refSeq).toUpperCase();
            String readStr = Utils.baseList2string(readBases).toUpperCase();

            boolean reverseFlag = read.getReadNegativeStrandFlag();
            for ( int i = 0; i < refStr.length(); i++) {
                if ( refStr.charAt(i) != readStr.charAt(i) )
                    mismatchCounts[(reverseFlag ? (refStr.length()-1-i) : i)]++;
            }
        }

        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone() {
        for ( int i = 0; i < mismatchCounts.length; i++ )
            System.out.println((i+1) + "\t" + mismatchCounts[i]);
    }

    private static Object resizeArray (Object oldArray, int newSize) {
        int oldSize = getLength(oldArray);
        Class elementType = oldArray.getClass().getComponentType();
        Object newArray = newInstance(elementType,newSize);
        int preserveLength = Math.min(oldSize,newSize);
        if (preserveLength > 0)
            System.arraycopy (oldArray,0,newArray,0,preserveLength);
        return newArray;
    }
}