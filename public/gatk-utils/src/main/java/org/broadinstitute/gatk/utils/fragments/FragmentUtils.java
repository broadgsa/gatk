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

package org.broadinstitute.gatk.utils.fragments;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.recalibration.EventType;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.util.*;

/**
 * An easy to access fragment-based pileup, which contains two separate pileups.  The first
 * is a regular collection of PileupElements containing all of the reads in the original RBP
 * that uniquely info about a fragment.  The second are TwoReadPileupElements that, as the
 * name suggests, contain two reads that are sequenced from the same underlying fragment.
 *
 * Based on the original code by E. Banks
 *
 * Oct 21: note that the order of the oneReadPileup and twoReadPileups are not
 * defined.  The algorithms that produce these lists are in fact producing
 * lists of Pileup elements *NOT* sorted by alignment start position of the underlying
 * reads.
 *
 * User: depristo
 * Date: 3/26/11
 * Time: 10:09 PM
 */
public final class FragmentUtils {

    public final static double DEFAULT_PCR_ERROR_RATE = 1e-4;
    public final static int DEFAULT_PCR_ERROR_QUAL = QualityUtil.getPhredScoreFromErrorProbability(DEFAULT_PCR_ERROR_RATE);
    public final static int HALF_OF_DEFAULT_PCR_ERROR_QUAL = DEFAULT_PCR_ERROR_QUAL / 2;

    protected final static byte MIN_QUAL_BAD_OVERLAP = 16;
    private FragmentUtils() {} // private constructor

    /**
     * A getter function that takes an Object of type T and returns its associated SAMRecord.
     *
     * Allows us to write a generic T -> Fragment algorithm that works with any object containing
     * a read.
     *
     * @param <T> The type of the object that contains a GATKSAMRecord
     */
    public interface ReadGetter<T> {
        /**
         * Get the GATKSAMRecord associated with object
         *
         * @param object the thing that contains the read
         * @return a non-null GATKSAMRecord read
         */
        public GATKSAMRecord get(T object);
    }

    /**
     * Identify getter for SAMRecords themselves
     */
    private final static ReadGetter<GATKSAMRecord> SamRecordGetter = new ReadGetter<GATKSAMRecord>() {
        @Override public GATKSAMRecord get(final GATKSAMRecord object) { return object; }
    };

    /**
     * Gets the SAMRecord in a PileupElement
     */
    private final static ReadGetter<PileupElement> PileupElementGetter = new ReadGetter<PileupElement>() {
        @Override public GATKSAMRecord get(final PileupElement object) { return object.getRead(); }
    };


    /**
     * Generic algorithm that takes an iterable over T objects, a getter routine to extract the reads in T,
     * and returns a FragmentCollection that contains the T objects whose underlying reads either overlap (or
     * not) with their mate pairs.
     *
     * @param readContainingObjects An iterator of objects that contain GATKSAMRecords
     * @param nElements the number of elements to be provided by the iterator, which is usually known upfront and
     *                  greatly improves the efficiency of the fragment calculation
     * @param getter a helper function that takes an object of type T and returns is associated GATKSAMRecord
     * @param <T>
     * @return a fragment collection
     */
    @Requires({
            "readContainingObjects != null",
            "nElements >= 0",
            "getter != null"
    })
    @Ensures("result != null")
    private static <T> FragmentCollection<T> create(final Iterable<T> readContainingObjects, final int nElements, final ReadGetter<T> getter) {
        Collection<T> singletons = null;
        Collection<List<T>> overlapping = null;
        Map<String, T> nameMap = null;

        int lastStart = -1;

        // build an initial map, grabbing all of the multi-read fragments
        for ( final T p : readContainingObjects ) {
            final SAMRecord read = getter.get(p);

            if ( read.getAlignmentStart() < lastStart ) {
                throw new IllegalArgumentException(String.format(
                        "FragmentUtils.create assumes that the incoming objects are ordered by " +
                                "SAMRecord alignment start, but saw a read %s with alignment start " +
                                "%d before the previous start %d", read.getSAMString(), read.getAlignmentStart(), lastStart));
            }
            lastStart = read.getAlignmentStart();

            final int mateStart = read.getMateAlignmentStart();
            if ( mateStart == 0 || mateStart > read.getAlignmentEnd() ) {
                // if we know that this read won't overlap its mate, or doesn't have one, jump out early
                if ( singletons == null ) singletons = new ArrayList<T>(nElements); // lazy init
                singletons.add(p);
            } else {
                // the read might overlap it's mate, or is the rightmost read of a pair
                final String readName = read.getReadName();
                final T pe1 = nameMap == null ? null : nameMap.get(readName);
                if ( pe1 != null ) {
                    // assumes we have at most 2 reads per fragment
                    if ( overlapping == null ) overlapping = new ArrayList<List<T>>(); // lazy init
                    overlapping.add(Arrays.asList(pe1, p));
                    nameMap.remove(readName);
                } else {
                    if ( nameMap == null ) nameMap = new HashMap<String, T>(nElements); // lazy init
                    nameMap.put(readName, p);
                }
            }
        }

        // add all of the reads that are potentially overlapping but whose mate never showed
        // up to the oneReadPile
        if ( nameMap != null && ! nameMap.isEmpty() ) {
            if ( singletons == null )
                singletons = nameMap.values();
            else
                singletons.addAll(nameMap.values());
        }

        return new FragmentCollection<T>(singletons, overlapping);
    }

    /**
     * Create a FragmentCollection containing PileupElements from the ReadBackedPileup rbp
     * @param rbp a non-null read-backed pileup.  The elements in this ReadBackedPileup must be ordered
     * @return a non-null FragmentCollection
     */
    @Ensures("result != null")
    public static FragmentCollection<PileupElement> create(final ReadBackedPileup rbp) {
        if ( rbp == null ) throw new IllegalArgumentException("Pileup cannot be null");
        return create(rbp, rbp.getNumberOfElements(), PileupElementGetter);
    }

    /**
     * Create a FragmentCollection containing GATKSAMRecords from a list of reads
     *
     * @param reads a non-null list of reads, ordered by their start location
     * @return a non-null FragmentCollection
     */
    @Ensures("result != null")
    public static FragmentCollection<GATKSAMRecord> create(final List<GATKSAMRecord> reads) {
        if ( reads == null ) throw new IllegalArgumentException("Pileup cannot be null");
        return create(reads, reads.size(), SamRecordGetter);
    }

    public static void adjustQualsOfOverlappingPairedFragments( final List<GATKSAMRecord> overlappingPair ) {
        if( overlappingPair.size() != 2 ) { throw new ReviewedGATKException("Found overlapping pair with " + overlappingPair.size() + " reads, but expecting exactly 2."); }

        final GATKSAMRecord firstRead = overlappingPair.get(0);
        final GATKSAMRecord secondRead = overlappingPair.get(1);

        if ( secondRead.getSoftStart() < firstRead.getSoftStart() ) {
            adjustQualsOfOverlappingPairedFragments(secondRead, firstRead);
        } else {
            adjustQualsOfOverlappingPairedFragments(firstRead, secondRead);
        }
    }

    /**
     * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible
     *
     * firstRead and secondRead must be part of the same fragment (though this isn't checked).  Looks
     * at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
     * are not treated independently.
     *
     * Assumes that firstRead starts before secondRead (according to their soft clipped starts)
     *
     * @param clippedFirstRead the left most read
     * @param clippedSecondRead the right most read
     *
     * @return a strandless merged read of first and second, or null if the algorithm cannot create a meaningful one
     */
    public static void adjustQualsOfOverlappingPairedFragments(final GATKSAMRecord clippedFirstRead, final GATKSAMRecord clippedSecondRead) {
        if ( clippedFirstRead == null ) throw new IllegalArgumentException("clippedFirstRead cannot be null");
        if ( clippedSecondRead == null ) throw new IllegalArgumentException("clippedSecondRead cannot be null");
        if ( ! clippedFirstRead.getReadName().equals(clippedSecondRead.getReadName()) ) throw new IllegalArgumentException("attempting to merge two reads with different names " + clippedFirstRead + " and " + clippedSecondRead);

        // don't adjust fragments that do not overlap
        if ( clippedFirstRead.getAlignmentEnd() < clippedSecondRead.getAlignmentStart() || !clippedFirstRead.getReferenceIndex().equals(clippedSecondRead.getReferenceIndex()) )
            return;

        final Pair<Integer, Boolean> pair = ReadUtils.getReadCoordinateForReferenceCoordinate(clippedFirstRead, clippedSecondRead.getAlignmentStart());
        final int firstReadStop = ( pair.getSecond() ? pair.getFirst() + 1 : pair.getFirst() );
        final int numOverlappingBases = Math.min(clippedFirstRead.getReadLength() - firstReadStop, clippedSecondRead.getReadLength());

        final byte[] firstReadBases = clippedFirstRead.getReadBases();
        final byte[] firstReadQuals = clippedFirstRead.getBaseQualities();
        final byte[] secondReadBases = clippedSecondRead.getReadBases();
        final byte[] secondReadQuals = clippedSecondRead.getBaseQualities();

        for ( int i = 0; i < numOverlappingBases; i++ ) {
            final int firstReadIndex = firstReadStop + i;
            final byte firstReadBase = firstReadBases[firstReadIndex];
            final byte secondReadBase = secondReadBases[i];

            if ( firstReadBase == secondReadBase ) {
                firstReadQuals[firstReadIndex] = (byte) Math.min(firstReadQuals[firstReadIndex], HALF_OF_DEFAULT_PCR_ERROR_QUAL);
                secondReadQuals[i] = (byte) Math.min(secondReadQuals[i], HALF_OF_DEFAULT_PCR_ERROR_QUAL);
            } else {
                // TODO -- use the proper statistical treatment of the quals from DiploidSNPGenotypeLikelihoods.java
                firstReadQuals[firstReadIndex] = 0;
                secondReadQuals[i] = 0;
            }
        }

        clippedFirstRead.setBaseQualities(firstReadQuals);
        clippedSecondRead.setBaseQualities(secondReadQuals);
    }

    public static List<GATKSAMRecord> mergeOverlappingPairedFragments( final List<GATKSAMRecord> overlappingPair ) {
        if( overlappingPair.size() != 2 ) { throw new ReviewedGATKException("Found overlapping pair with " + overlappingPair.size() + " reads, but expecting exactly 2."); }

        final GATKSAMRecord firstRead = overlappingPair.get(0);
        final GATKSAMRecord secondRead = overlappingPair.get(1);

        final GATKSAMRecord merged;
        if( !(secondRead.getSoftStart() <= firstRead.getSoftEnd() && secondRead.getSoftStart() >= firstRead.getSoftStart() && secondRead.getSoftEnd() >= firstRead.getSoftEnd()) ) {
            merged = mergeOverlappingPairedFragments(secondRead, firstRead);
        } else {
            merged = mergeOverlappingPairedFragments(firstRead, secondRead);
        }

        return merged == null ? overlappingPair : Collections.singletonList(merged);
    }

    /**
     * Merge two overlapping reads from the same fragment into a single super read, if possible
     *
     * firstRead and secondRead must be part of the same fragment (though this isn't checked).  Looks
     * at the bases and alignment, and tries its best to create a meaningful synthetic single super read
     * that represents the entire sequenced fragment.
     *
     * Assumes that firstRead starts before secondRead (according to their soft clipped starts)
     *
     * @param unclippedFirstRead the left most read
     * @param unclippedSecondRead the right most read
     *
     * @return a strandless merged read of first and second, or null if the algorithm cannot create a meaningful one
     */
    public static GATKSAMRecord mergeOverlappingPairedFragments(final GATKSAMRecord unclippedFirstRead, final GATKSAMRecord unclippedSecondRead) {
        if ( unclippedFirstRead == null ) throw new IllegalArgumentException("unclippedFirstRead cannot be null");
        if ( unclippedSecondRead == null ) throw new IllegalArgumentException("unclippedSecondRead cannot be null");
        if ( ! unclippedFirstRead.getReadName().equals(unclippedSecondRead.getReadName()) ) throw new IllegalArgumentException("attempting to merge two reads with different names " + unclippedFirstRead + " and " + unclippedSecondRead);

        if( unclippedFirstRead.getCigarString().contains("I") || unclippedFirstRead.getCigarString().contains("D") || unclippedSecondRead.getCigarString().contains("I") || unclippedSecondRead.getCigarString().contains("D") ) {
            return null; // fragments contain indels so don't merge them
        }

        final GATKSAMRecord firstRead = ReadClipper.hardClipAdaptorSequence(ReadClipper.revertSoftClippedBases(unclippedFirstRead));
        final GATKSAMRecord secondRead = ReadClipper.hardClipAdaptorSequence(ReadClipper.revertSoftClippedBases(unclippedSecondRead));

        if( !(secondRead.getSoftStart() <= firstRead.getSoftEnd() && secondRead.getSoftStart() >= firstRead.getSoftStart() && secondRead.getSoftEnd() >= firstRead.getSoftEnd()) ) {
            return null; // can't merge them, yet:  AAAAAAAAAAA-BBBBBBBBBBB-AAAAAAAAAAAAAA, B is contained entirely inside A
        }

        final Pair<Integer, Boolean> pair = ReadUtils.getReadCoordinateForReferenceCoordinate(firstRead, secondRead.getAlignmentStart());

        final int firstReadStop = ( pair.getSecond() ? pair.getFirst() + 1 : pair.getFirst() );
        final int numBases = firstReadStop + secondRead.getReadLength();
        final byte[] bases = new byte[numBases];
        final byte[] quals = new byte[numBases];
        final byte[] insertionQuals = new byte[numBases];
        final byte[] deletionQuals = new byte[numBases];
        final byte[] firstReadBases = firstRead.getReadBases();
        final byte[] firstReadQuals = firstRead.getBaseQualities();
        final byte[] secondReadBases = secondRead.getReadBases();
        final byte[] secondReadQuals = secondRead.getBaseQualities();

        for(int iii = 0; iii < firstReadStop; iii++) {
            bases[iii] = firstReadBases[iii];
            quals[iii] = firstReadQuals[iii];
        }
        for(int iii = firstReadStop; iii < firstRead.getReadLength(); iii++) {
            if( firstReadQuals[iii] > MIN_QUAL_BAD_OVERLAP && secondReadQuals[iii-firstReadStop] > MIN_QUAL_BAD_OVERLAP && firstReadBases[iii] != secondReadBases[iii-firstReadStop] ) {
                return null; // high qual bases don't match exactly, probably indel in only one of the fragments, so don't merge them
            }
            if( firstReadQuals[iii] < MIN_QUAL_BAD_OVERLAP && secondReadQuals[iii-firstReadStop] < MIN_QUAL_BAD_OVERLAP ) {
                return null; // both reads have low qual bases in the overlap region so don't merge them because don't know what is going on
            }
            bases[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadBases[iii] : secondReadBases[iii-firstReadStop] );
            quals[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadQuals[iii] : secondReadQuals[iii-firstReadStop] );
        }
        for(int iii = firstRead.getReadLength(); iii < numBases; iii++) {
            bases[iii] = secondReadBases[iii-firstReadStop];
            quals[iii] = secondReadQuals[iii-firstReadStop];
        }

        final GATKSAMRecord returnRead = new GATKSAMRecord( firstRead.getHeader() );
        returnRead.setIsStrandless(true);
        returnRead.setAlignmentStart( firstRead.getAlignmentStart() );
        returnRead.setReadBases( bases );
        returnRead.setBaseQualities( quals );
        returnRead.setReadGroup( firstRead.getReadGroup() );
        returnRead.setReferenceName( firstRead.getReferenceName() );
        returnRead.setReadName( firstRead.getReadName() );
        final CigarElement c = new CigarElement(bases.length, CigarOperator.M);
        final ArrayList<CigarElement> cList = new ArrayList<CigarElement>();
        cList.add(c);
        returnRead.setCigar( new Cigar( cList ));
        returnRead.setMappingQuality( firstRead.getMappingQuality() );

        if( firstRead.hasBaseIndelQualities() || secondRead.hasBaseIndelQualities() ) {
            final byte[] firstReadInsertionQuals = firstRead.getBaseInsertionQualities();
            final byte[] firstReadDeletionQuals = firstRead.getBaseDeletionQualities();
            final byte[] secondReadInsertionQuals = secondRead.getBaseInsertionQualities();
            final byte[] secondReadDeletionQuals = secondRead.getBaseDeletionQualities();
            for(int iii = 0; iii < firstReadStop; iii++) {
                insertionQuals[iii] = firstReadInsertionQuals[iii];
                deletionQuals[iii] = firstReadDeletionQuals[iii];
            }
            for(int iii = firstReadStop; iii < firstRead.getReadLength(); iii++) {
                insertionQuals[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadInsertionQuals[iii] : secondReadInsertionQuals[iii-firstReadStop] ); // Purposefully checking the highest *base* quality score
                deletionQuals[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadDeletionQuals[iii] : secondReadDeletionQuals[iii-firstReadStop] ); // Purposefully checking the highest *base* quality score
            }
            for(int iii = firstRead.getReadLength(); iii < numBases; iii++) {
                insertionQuals[iii] = secondReadInsertionQuals[iii-firstReadStop];
                deletionQuals[iii] = secondReadDeletionQuals[iii-firstReadStop];
            }
            returnRead.setBaseQualities( insertionQuals, EventType.BASE_INSERTION );
            returnRead.setBaseQualities( deletionQuals, EventType.BASE_DELETION );
        }

        return returnRead;
    }
}
