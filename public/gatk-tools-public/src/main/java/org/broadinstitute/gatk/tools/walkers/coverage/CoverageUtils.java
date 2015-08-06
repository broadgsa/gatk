/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.coverage;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fragments.FragmentCollection;
import org.broadinstitute.gatk.utils.pileup.PileupElement;

import java.util.*;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Mar 3, 2010
 */
public class CoverageUtils {

    public enum CountPileupType {
        /**
         * Count all reads independently (even if from the same fragment).
         */
        COUNT_READS,
        /**
         * Count all fragments (even if the reads that compose the fragment are not consistent at that base).
         */
        COUNT_FRAGMENTS,
        /**
         * Count all fragments (but only if the reads that compose the fragment are consistent at that base).
         */
        COUNT_FRAGMENTS_REQUIRE_SAME_BASE
    }

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
            throw new ReviewedGATKException("Invalid type ID sent to getTypeID. This is a BUG!");
        }
    }

    public static Map<DoCOutputType.Partition,Map<String,int[]>>
                    getBaseCountsByPartition(AlignmentContext context, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ, CountPileupType countType, Collection<DoCOutputType.Partition> types) {

        Map<DoCOutputType.Partition,Map<String,int[]>> countsByIDByType = new HashMap<DoCOutputType.Partition,Map<String,int[]>>();
        Map<SAMReadGroupRecord,int[]> countsByRG = getBaseCountsByReadGroup(context,minMapQ,maxMapQ,minBaseQ,maxBaseQ,countType);
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

    public static Map<SAMReadGroupRecord,int[]> getBaseCountsByReadGroup(AlignmentContext context, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ, CountPileupType countType) {
        Map<SAMReadGroupRecord, int[]> countsByRG = new HashMap<SAMReadGroupRecord,int[]>();

        Map<String, int[]> countsByRGName = new HashMap<String, int[]>();
        Map<String, SAMReadGroupRecord> RGByName = new HashMap<String, SAMReadGroupRecord>();

        List<PileupElement> countPileup = new LinkedList<PileupElement>();
        FragmentCollection<PileupElement> fpile;

        switch (countType) {

            case COUNT_READS:
                for (PileupElement e : context.getBasePileup())
                    if (countElement(e, minMapQ, maxMapQ, minBaseQ, maxBaseQ))
                        countPileup.add(e);
                break;

            case COUNT_FRAGMENTS: // ignore base identities and put in FIRST base that passes filters:
                fpile = context.getBasePileup().getStartSortedPileup().toFragments();

                for (PileupElement e : fpile.getSingletonReads())
                    if (countElement(e, minMapQ, maxMapQ, minBaseQ, maxBaseQ))
                        countPileup.add(e);

                for (List<PileupElement> overlappingPair : fpile.getOverlappingPairs()) {
                    // iterate over all elements in fragment:
                    for (PileupElement e : overlappingPair) {
                        if (countElement(e, minMapQ, maxMapQ, minBaseQ, maxBaseQ)) {
                            countPileup.add(e); // add the first passing element per fragment
                            break;
                        }
                    }
                }
                break;

            case COUNT_FRAGMENTS_REQUIRE_SAME_BASE:
                fpile = context.getBasePileup().getStartSortedPileup().toFragments();

                for (PileupElement e : fpile.getSingletonReads())
                    if (countElement(e, minMapQ, maxMapQ, minBaseQ, maxBaseQ))
                        countPileup.add(e);

                for (List<PileupElement> overlappingPair : fpile.getOverlappingPairs()) {
                    PileupElement firstElem = null;
                    PileupElement addElem = null;

                    // iterate over all elements in fragment:
                    for (PileupElement e : overlappingPair) {
                        if (firstElem == null)
                            firstElem = e;
                        else if (e.getBase() != firstElem.getBase()) {
                            addElem = null;
                            break;
                        }

                        // will add the first passing element per base-consistent fragment:
                        if (addElem == null && countElement(e, minMapQ, maxMapQ, minBaseQ, maxBaseQ))
                            addElem = e;
                    }

                    if (addElem != null)
                        countPileup.add(addElem);
                }
                break;

            default:
                throw new UserException("Must use valid CountPileupType");
        }

        for (PileupElement e : countPileup) {
            SAMReadGroupRecord readGroup = getReadGroup(e.getRead());

            String readGroupId = readGroup.getSample() + "_" + readGroup.getReadGroupId();
            int[] counts = countsByRGName.get(readGroupId);
            if (counts == null) {
                counts = new int[6];
                countsByRGName.put(readGroupId, counts);
                RGByName.put(readGroupId, readGroup);
            }

            updateCounts(counts, e);
        }

        for (String readGroupId : RGByName.keySet()) {
            countsByRG.put(RGByName.get(readGroupId), countsByRGName.get(readGroupId));
        }

        return countsByRG;
    }

    private static boolean countElement(PileupElement e, int minMapQ, int maxMapQ, byte minBaseQ, byte maxBaseQ) {
        return (e.getMappingQual() >= minMapQ && e.getMappingQual() <= maxMapQ && ( e.getQual() >= minBaseQ && e.getQual() <= maxBaseQ || e.isDeletion() ));
    }

    private static void updateCounts(int[] counts, PileupElement e) {
        if ( e.isDeletion() ) {
            counts[BaseUtils.Base.D.ordinal()]++;
        } else if ( BaseUtils.basesAreEqual(BaseUtils.Base.N.base, e.getBase()) ) {
            counts[BaseUtils.Base.N.ordinal()]++;
        } else {
            try {
                counts[BaseUtils.simpleBaseToBaseIndex(e.getBase())]++;
            } catch (ArrayIndexOutOfBoundsException exc) {
                throw new ReviewedGATKException("Expected a simple base, but actually received"+(char)e.getBase());
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
