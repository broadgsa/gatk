package org.broadinstitute.sting.playground.contexts;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.Pair;

import java.util.List;
import java.util.ArrayList;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 9, 2009
 * Time: 11:17:30 AM
 * To change this template use File | Settings | File Templates.
 */
public class QualityScoreThresholdContext extends FilteredAlignmentContext{
    /*
     * @Param: qThresh - default value for thresholding
     */
    protected byte qThresh = 22;

    public QualityScoreThresholdContext(AlignmentContext context, byte qThresh) {
        this.qThresh = qThresh;
        Pair<List<SAMRecord>, List<Integer>> filteredRO = filter(context);
        this.reads = filteredRO.getFirst();
        this.offsets = filteredRO.getSecond();
        this.loc = context.getLocation();
    }
    
    public byte getQualityScoreThreshold() {
        return this.qThresh;
    }

    public Pair<List<SAMRecord>,List<Integer>> filter(AlignmentContext context) {
        List<SAMRecord> inReads = context.getReads();
        List<Integer> inOffsets = context.getOffsets();
        List<SAMRecord> outReads = new ArrayList<SAMRecord>();
        List<Integer> outOffsets = new ArrayList<Integer>();

        for( int i = 0; i < inReads.size(); i++) {
            if(inReads.get(i).getBaseQualities()[inOffsets.get(i)] >= this.qThresh) {
                outReads.add(inReads.get(i));
                outOffsets.add(inOffsets.get(i));
            }
        }

        return new Pair<List<SAMRecord>,List<Integer>>(outReads,outOffsets);
    }
}
