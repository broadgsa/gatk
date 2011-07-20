package org.broadinstitute.sting.gatk.walkers.coverage;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

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

    public static String getTypeID( SAMReadGroupRecord r, DoCOutputType.Partition type ) {
        if ( type == DoCOutputType.Partition.sample ) {
            return r.getSample();
        } else if ( type == DoCOutputType.Partition.readgroup ) {
            return String.format("%s_rg_%s",r.getSample(),r.getReadGroupId());
        } else if ( type == DoCOutputType.Partition.library ) {
            return r.getLibrary();
        } else if ( type == DoCOutputType.Partition.center ) {
            return r.getSequencingCenter();
        } else if ( type == DoCOutputType.Partition.platform ) {
            return r.getPlatform();
        } else if ( type == DoCOutputType.Partition.sample_by_center ) {
            return String.format("%s_cn_%s",r.getSample(),r.getSequencingCenter());
        } else if ( type == DoCOutputType.Partition.sample_by_platform) {
            return String.format("%s_pl_%s",r.getSample(),r.getPlatform());
        } else if ( type == DoCOutputType.Partition.sample_by_platform_by_center ) {
            return String.format("%s_pl_%s_cn_%s",r.getSample(),r.getPlatform(),r.getSequencingCenter());
        } else {
            throw new ReviewedStingException("Invalid type ID sent to getTypeID. This is a BUG!");
        }
    }

    public static Map<DoCOutputType.Partition,Map<String,int[]>>
                    getBaseCountsByPartition(AlignmentContext context, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ, Collection<DoCOutputType.Partition> types) {

        Map<DoCOutputType.Partition,Map<String,int[]>> countsByIDByType = new HashMap<DoCOutputType.Partition,Map<String,int[]>>();
        Map<SAMReadGroupRecord,int[]> countsByRG = getBaseCountsByReadGroup(context,minMapQ,maxMapQ,minBaseQ,maxBaseQ);
        for (DoCOutputType.Partition t : types ) {
            // iterate through the read group counts and build the type associations
            for ( Map.Entry<SAMReadGroupRecord,int[]> readGroupCountEntry : countsByRG.entrySet() ) {
                String typeID = getTypeID(readGroupCountEntry.getKey(),t);

                if ( ! countsByIDByType.keySet().contains(t) ) {
                    countsByIDByType.put(t,new HashMap<String,int[]>());
                }

                if ( ! countsByIDByType.get(t).keySet().contains(typeID) ) {
                    countsByIDByType.get(t).put(typeID,readGroupCountEntry.getValue().clone());
                } else {
                    addCounts(countsByIDByType.get(t).get(typeID),readGroupCountEntry.getValue());
                }
            }
        }


        return countsByIDByType;
    }

    public static void addCounts(int[] updateMe, int[] leaveMeAlone ) {
        for ( int index = 0; index < leaveMeAlone.length; index++ ) {
            updateMe[index] += leaveMeAlone[index];
        }
    }

    public static Map<SAMReadGroupRecord,int[]> getBaseCountsByReadGroup(AlignmentContext context, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ) {
        Map<SAMReadGroupRecord, int[]> countsByRG = new HashMap<SAMReadGroupRecord,int[]>();
        for ( PileupElement e : context.getBasePileup() ) {
            if ( e.getMappingQual() >= minMapQ && e.getMappingQual() <= maxMapQ && ( e.getQual() >= minBaseQ && e.getQual() <= maxBaseQ || e.isDeletion() ) ) {
                SAMReadGroupRecord readGroup = getReadGroup(e.getRead());
                if ( ! countsByRG.keySet().contains(readGroup) ) {
                    countsByRG.put(readGroup,new int[6]);
                    updateCounts(countsByRG.get(readGroup),e);
                } else {
                    updateCounts(countsByRG.get(readGroup),e);
                }
            }
        }

        return countsByRG;
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
                throw new ReviewedStingException("Expected a simple base, but actually received"+(char)e.getBase());
            }
        }
    }

    private static SAMReadGroupRecord getReadGroup(SAMRecord r) {
        SAMReadGroupRecord rg = r.getReadGroup();
        if ( rg == null ) {
            String msg = "Read "+r.getReadName()+" lacks read group information; Please associate all reads with read groups";
            throw new UserException.MalformedBAM(r, msg);
        }

        return rg;
    }
}
