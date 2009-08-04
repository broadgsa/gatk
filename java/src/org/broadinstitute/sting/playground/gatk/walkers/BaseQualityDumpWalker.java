package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;

@WalkerName("Base_Quality_Dump")
public class BaseQualityDumpWalker extends ReadWalker<Integer, Integer> {

    protected final int MIN_TARGET_EDIT_DISTANCE = 0; //5;
    protected final int MAX_TARGET_EDIT_DISTANCE = 4; //10;

    // Do we actually want to operate on the context?
    public boolean filter(char[] ref, SAMRecord read) {
	    // we only want aligned reads
	    return !read.getReadUnmappedFlag();
    }

    public Integer map(char[] ref, SAMRecord read) {

        int editDist = Integer.parseInt(read.getAttribute("NM").toString());

        // ignore alignments with indels for now
        if ( read.getAlignmentBlocks().size() == 1 &&
             editDist >= MIN_TARGET_EDIT_DISTANCE &&
             editDist <= MAX_TARGET_EDIT_DISTANCE ) {

            String qualStr = read.getBaseQualityString();
            int[] scores = new int[qualStr.length()];
            boolean reverseFlag = read.getReadNegativeStrandFlag();
            for ( int i = 0; i < qualStr.length(); i++)
                scores[(reverseFlag ? (qualStr.length()-1-i) : i)] += (int)qualStr.charAt(i) - 33;
            for ( int i = 0; i < scores.length; i++ )
                out.print(scores[i] + " ");
            out.println("");
        }

        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
