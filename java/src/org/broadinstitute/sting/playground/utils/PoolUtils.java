package org.broadinstitute.sting.playground.utils;

import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

import org.broadinstitute.sting.utils.Pair;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 27, 2009
 * Time: 12:31:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class PoolUtils {

    private PoolUtils () {}

    public static Pair<Pair<List<SAMRecord>, List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> splitReadsByReadDirection(List<SAMRecord> reads, List<Integer> offsets) {
        ArrayList<SAMRecord> forwardReads;
        ArrayList<SAMRecord> reverseReads;
        ArrayList<Integer> forwardOffsets;
        ArrayList<Integer> reverseOffsets;

        if ( reads == null) {
            forwardReads = null;
            reverseReads = null;
            forwardOffsets = null;
            reverseOffsets = null;
        } else {
            forwardReads = new ArrayList();
            reverseReads = new ArrayList();
            forwardOffsets = new ArrayList();
            reverseOffsets = new ArrayList();

            for ( int readNo = 0; readNo < reads.size(); readNo ++ ) {
                if ( reads.get(readNo).getReadNegativeStrandFlag() ) {
                    forwardReads.add(reads.get(readNo));
                    forwardOffsets.add(offsets.get(readNo));
                } else {
                    reverseReads.add(reads.get(readNo));
                    reverseOffsets.add(offsets.get(readNo));
                }
            }
        }

        return new Pair(new Pair(forwardReads,reverseReads), new Pair(forwardOffsets,reverseOffsets));
    }
    
}
