package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.HashMap;
import java.util.Map;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Mar 3, 2010
 */
public class CoverageUtils {

    public enum PartitionType { BY_READ_GROUP, BY_SAMPLE }

    /**
     * Returns the counts of bases from reads with MAPQ > minMapQ and base quality > minBaseQ in the context
     * as an array of ints, indexed by the index fields of BaseUtils
     *
     * @param context
     * @param minMapQ
     * @param minBaseQ
     * @return
     */
    public static int[] getBaseCounts(AlignmentContext context, int minMapQ, int minBaseQ) {
        int[] counts = new int[6];

        for (PileupElement e : context.getBasePileup()) {
            if ( e.getMappingQual() >= minMapQ && ( e.getQual() >= minBaseQ || e.isDeletion() ) ) {
                updateCounts(counts,e);
            }
        }

        return counts;
    }


    public static Map<String,int[]> getBaseCountsBySample(AlignmentContext context, int minMapQ, int minBaseQ, PartitionType type) {
        Map<String,int[]> samplesToCounts = new HashMap<String,int[]>();

        for (PileupElement e : context.getBasePileup()) {
            if ( e.getMappingQual() >= minMapQ && ( e.getQual() >= minBaseQ || e.isDeletion() ) ) {
                String sample = type == PartitionType.BY_SAMPLE ? e.getRead().getReadGroup().getSample() :
                        e.getRead().getReadGroup().getSample()+"_rg_"+e.getRead().getReadGroup().getReadGroupId();
                if ( ! samplesToCounts.keySet().contains(sample) )
                    samplesToCounts.put(sample,new int[6]);
                updateCounts(samplesToCounts.get(sample),e);
            }
        }

        return samplesToCounts;
    }

    private static void updateCounts(int[] counts, PileupElement e) {
        if ( e.isDeletion() ) {
            counts[BaseUtils.DELETION_INDEX]++;
        } else if ( BaseUtils.basesAreEqual((byte) 'N', e.getBase()) ) {
            counts[BaseUtils.NO_CALL_INDEX]++;
        } else {
            try {
                counts[BaseUtils.simpleBaseToBaseIndex(e.getBase())]++;
            } catch (ArrayIndexOutOfBoundsException exc) {
                throw new StingException("Expected a simple base, but actually received"+(char)e.getBase());
            }
        }
    }
}
