package org.broadinstitute.sting.gatk.walkers.coverage;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Mar 3, 2010
 */
public class CoverageUtils {

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

    public static String getTypeID( SAMRecord r, CoverageAggregator.AggregationType type ) {
        if ( type == CoverageAggregator.AggregationType.SAMPLE ) {
            return r.getReadGroup().getSample();
        } else if ( type == CoverageAggregator.AggregationType.READGROUP ) {
            return String.format("%s_rg_%s",r.getReadGroup().getSample(),r.getReadGroup().getReadGroupId());
        } else if ( type == CoverageAggregator.AggregationType.LIBRARY ) {
            return r.getReadGroup().getLibrary();
        } else {
            throw new StingException("Invalid type ID sent to getTypeID. This is a BUG!");
        }
    }
    public static Map<CoverageAggregator.AggregationType,Map<String,int[]>>
                    getBaseCountsByPartition(AlignmentContext context, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ, List<CoverageAggregator.AggregationType> types) {

        Map<CoverageAggregator.AggregationType,Map<String,int[]>> countsByIDByType = new HashMap<CoverageAggregator.AggregationType,Map<String,int[]>>();
        for (CoverageAggregator.AggregationType t : types ) {
            countsByIDByType.put(t,new HashMap<String,int[]>());
        }
        for (PileupElement e : context.getBasePileup()) {
            if ( e.getMappingQual() >= minMapQ && e.getMappingQual() <= maxMapQ && ( ( e.getQual() >= minBaseQ && e.getQual() <= maxBaseQ ) || e.isDeletion() ) ) {
                for (CoverageAggregator.AggregationType t : types ) {
                    String id = getTypeID(e.getRead(),t);
                    if ( ! countsByIDByType.get(t).keySet().contains(id) ) {
                        countsByIDByType.get(t).put(id,new int[6]);
                    }
                    updateCounts(countsByIDByType.get(t).get(id),e);
                }
            }
        }

        return countsByIDByType;
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
