package org.broadinstitute.sting.playground.contexts;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 9, 2009
 * Time: 11:01:53 AM
 * To change this template use File | Settings | File Templates.
 */
public class ForwardReadsContext extends FilteredAlignmentContext {

    public Pair<List<SAMRecord>,List<Integer>> filter(AlignmentContext context) {
        List<SAMRecord> inReads = context.getReads();
        List<Integer> inOffsets = context.getOffsets();
        List<SAMRecord> filteredReads = new ArrayList<SAMRecord>();
        List<Integer> filteredOffsets = new ArrayList<Integer>();

        for( int i = 0; i < inReads.size(); i++ ) {
            if( ! inReads.get(i).getReadNegativeStrandFlag() ) {
                filteredReads.add(inReads.get(i));
                filteredOffsets.add(inOffsets.get(i));
            }
        }

        return new Pair<List<SAMRecord>,List<Integer>>(filteredReads,filteredOffsets);
    }
}
