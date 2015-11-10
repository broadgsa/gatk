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

package org.broadinstitute.gatk.utils;

import java.math.BigInteger;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Sep 10, 2010
 * Time: 1:56:24 PM
 *
 * A series of utility functions that enable the GATK to compare two sequence dictionaries -- from the reference,
 * from BAMs, or from RODs -- for consistency.  The system supports two basic modes: get an enum state that
 * describes at a high level the consistency between two dictionaries, or a validateDictionaries that will
 * blow up with a UserException if the dicts are too incompatible.
 *
 * Dictionaries are tested for contig name overlaps, consistency in ordering in these overlap set, and length,
 * if available.  Examines the Engine arguments to decided if the -U option to allow danger seq dict inconsistency
 * is enabled before it blows up.
 */
public class SequenceDictionaryUtils {
    //
    // for detecting lexicographically sorted human references
    //
    private static final boolean ENABLE_LEXICOGRAPHIC_REQUIREMENT_FOR_HUMAN = true;

    // hg18
    protected static final SAMSequenceRecord CHR1_HG18 = new SAMSequenceRecord("chr1", 247249719);
    protected static final SAMSequenceRecord CHR2_HG18 = new SAMSequenceRecord("chr2", 242951149);
    protected static final SAMSequenceRecord CHR10_HG18 = new SAMSequenceRecord("chr10", 135374737);

    // hg19
    protected static final SAMSequenceRecord CHR1_HG19 = new SAMSequenceRecord("chr1", 249250621);
    protected static final SAMSequenceRecord CHR2_HG19 = new SAMSequenceRecord("chr2", 243199373);
    protected static final SAMSequenceRecord CHR10_HG19 = new SAMSequenceRecord("chr10", 135534747);

    // b36
    protected static final SAMSequenceRecord CHR1_B36 = new SAMSequenceRecord("1", 247249719);
    protected static final SAMSequenceRecord CHR2_B36 = new SAMSequenceRecord("2", 242951149);
    protected static final SAMSequenceRecord CHR10_B36 = new SAMSequenceRecord("10", 135374737);

    // b37
    protected static final SAMSequenceRecord CHR1_B37 = new SAMSequenceRecord("1", 249250621);
    protected static final SAMSequenceRecord CHR2_B37 = new SAMSequenceRecord("2", 243199373);
    protected static final SAMSequenceRecord CHR10_B37 = new SAMSequenceRecord("10", 135534747);


    public enum SequenceDictionaryCompatibility {
        IDENTICAL,                      // the dictionaries are identical
        COMMON_SUBSET,                  // there exists a common subset of equivalent contigs
        NO_COMMON_CONTIGS,              // no overlap between dictionaries
        UNEQUAL_COMMON_CONTIGS,         // common subset has contigs that have the same name but different lengths and/or MD5s
        NON_CANONICAL_HUMAN_ORDER,      // human reference detected but the order of the contigs is non-standard (lexicographic, for examine)
        OUT_OF_ORDER,                   // the two dictionaries overlap but the overlapping contigs occur in different
                                        // orders with respect to each other
        DIFFERENT_INDICES               // the two dictionaries overlap and the overlapping contigs occur in the same
                                        // order with respect to each other, but one or more of them have different
                                        // indices in the two dictionaries. Eg., { chrM, chr1, chr2 } vs. { chr1, chr2 }
    }

    /**
     * @param validationExclusion exclusions to validation
     * @return Returns true if the engine is in tolerant mode and we'll let through dangerous but not fatal dictionary inconsistency
     */
    private static boolean allowNonFatalIncompabilities(final ValidationExclusion.TYPE validationExclusion) {
        return ( validationExclusion == ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY ||
                        validationExclusion == ValidationExclusion.TYPE.ALL );
    }

    /**
     * Tests for compatibility between two sequence dictionaries.  If the dictionaries are incompatible, then
     * UserExceptions are thrown with detailed error messages.  If the engine is in permissive mode, then
     * logger warnings are generated instead.
     *
     * @param logger for warnings
     * @param validationExclusion exclusions to validation
     * @param name1 name associated with dict1
     * @param dict1 the sequence dictionary dict1
     * @param name2 name associated with dict2
     * @param dict2 the sequence dictionary dict2
     * @param isReadsToReferenceComparison true if one of the dictionaries comes from a reads data source (eg., a BAM),
     *                                     and the other from a reference data source
     * @param intervals the user-specified genomic intervals: only required when isReadsToReferenceComparison is true,
     *                  otherwise can be null
     */
    public static void validateDictionaries( final Logger logger,
                                             final ValidationExclusion.TYPE validationExclusion,
                                             final String name1,
                                             final SAMSequenceDictionary dict1,
                                             final String name2,
                                             final SAMSequenceDictionary dict2,
                                             final boolean isReadsToReferenceComparison,
                                             final GenomeLocSortedSet intervals ) {

        final SequenceDictionaryCompatibility type = compareDictionaries(dict1, dict2);

        switch ( type ) {
            case IDENTICAL:
                return;
            case COMMON_SUBSET:
                 return;
            case NO_COMMON_CONTIGS:
                throw new UserException.IncompatibleSequenceDictionaries("No overlapping contigs found", name1, dict1, name2, dict2);

            case UNEQUAL_COMMON_CONTIGS: {
                final List<SAMSequenceRecord> x = findNotEqualCommonContigs(getCommonContigsByName(dict1, dict2), dict1, dict2);
                final SAMSequenceRecord elt1 = x.get(0);
                final SAMSequenceRecord elt2 = x.get(1);

                String msg = "Found contigs with the same name but different lengths";
                String contig1  = "  contig  " + name1 + " is named " + elt1.getSequenceName()  + " with length " + Integer.toString(elt1.getSequenceLength());
                if ( elt1.getMd5() != null )
                    contig1 += " and MD5 " + elt1.getMd5();
                String contig2  = "  contig  " + name2 + " is named " + elt2.getSequenceName()  + " with length " + Integer.toString(elt2.getSequenceLength());
                if ( elt2.getMd5() != null )
                    contig2 += " and MD5 " + elt2.getMd5();
                if ( elt1.getMd5() != null ||  elt2.getMd5() != null )
                    msg += " or MD5s:";
                msg += "\n" + contig1 + "\n" + contig2;

                // todo -- replace with toString when SAMSequenceRecord has a nice toString routine
                final UserException ex = new UserException.IncompatibleSequenceDictionaries(msg, name1, dict1, name2, dict2);

                if ( allowNonFatalIncompabilities(validationExclusion) )
                    logger.warn(ex.getMessage());
                else
                    throw ex;
                break;
            }

            case NON_CANONICAL_HUMAN_ORDER: {
                UserException ex;
                if ( nonCanonicalHumanContigOrder(dict1) )
                    ex = new UserException.LexicographicallySortedSequenceDictionary(name1, dict1);
                else
                    ex = new UserException.LexicographicallySortedSequenceDictionary(name2, dict2);
                
                if ( allowNonFatalIncompabilities(validationExclusion) )
                    logger.warn(ex.getMessage());
                else
                    throw ex;
                break;
            }

            case OUT_OF_ORDER: {
                UserException ex = new UserException.IncompatibleSequenceDictionaries(
			"The contig order in " + name1 + " and " + name2 + "is not "
			+ "the same; to fix this please see: "
			+ "(https://www.broadinstitute.org/gatk/guide/article?id=1328), "
			+ " which describes reordering contigs in BAM and VCF files.",
			name1, dict1, name2, dict2);
                if ( allowNonFatalIncompabilities(validationExclusion) )
                    logger.warn(ex.getMessage());
                else
                    throw ex;
                break;
            }

            case DIFFERENT_INDICES: {
                // This is currently only known to be problematic when the index mismatch is between a bam and the
                // reference AND when the user's intervals actually include one or more of the contigs that are
                // indexed differently from the reference. In this case, the engine will fail to correctly serve
                // up the reads from those contigs, so throw an exception unless unsafe operations are enabled.
                if ( isReadsToReferenceComparison && intervals != null ) {

                     final Set<String> misindexedContigs = findMisindexedContigsInIntervals(intervals, dict1, dict2);

                     if ( ! misindexedContigs.isEmpty() ) {
                         final String msg = String.format("The following contigs included in the intervals to process have " +
                                                          "different indices in the sequence dictionaries for the reads vs. " +
                                                          "the reference: %s.  As a result, the GATK engine will not correctly " +
                                                          "process reads from these contigs. You should either fix the sequence " +
                                                          "dictionaries for your reads so that these contigs have the same indices " +
                                                          "as in the sequence dictionary for your reference, or exclude these contigs " +
                                                          "from your intervals. This error can be disabled via -U %s, " +
                                                          "however this is not recommended as the GATK engine will not behave correctly.",
                                                          misindexedContigs, ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY);
                         final UserException ex = new UserException.IncompatibleSequenceDictionaries(msg, name1, dict1, name2, dict2);

                         if ( allowNonFatalIncompabilities(validationExclusion) )
                             logger.warn(ex.getMessage());
                         else
                             throw ex;
                     }
                }
                break;
            }

            default:
                throw new ReviewedGATKException("Unexpected SequenceDictionaryComparison type: " + type);
        }
    }

    /**
     * Workhorse routine that takes two dictionaries and returns their compatibility.
     *
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return A SequenceDictionaryCompatibility enum value describing the compatibility of the two dictionaries
     */
    public static SequenceDictionaryCompatibility compareDictionaries( final SAMSequenceDictionary dict1, final SAMSequenceDictionary dict2) {
        if ( nonCanonicalHumanContigOrder(dict1) || nonCanonicalHumanContigOrder(dict2) )
            return SequenceDictionaryCompatibility.NON_CANONICAL_HUMAN_ORDER;

        final Set<String> commonContigs = getCommonContigsByName(dict1, dict2);

        if (commonContigs.isEmpty())
            return SequenceDictionaryCompatibility.NO_COMMON_CONTIGS;
        else if ( ! commonContigsHaveSameLengths(commonContigs, dict1, dict2) )
            return SequenceDictionaryCompatibility.UNEQUAL_COMMON_CONTIGS;
        else if ( ! commonContigsAreInSameRelativeOrder(commonContigs, dict1, dict2) )
            return SequenceDictionaryCompatibility.OUT_OF_ORDER;
        else if ( commonContigs.size() == dict1.size() && commonContigs.size() == dict2.size() )
            return SequenceDictionaryCompatibility.IDENTICAL;
        else if ( ! commonContigsAreAtSameIndices(commonContigs, dict1, dict2) )
            return SequenceDictionaryCompatibility.DIFFERENT_INDICES;
        else {
            return SequenceDictionaryCompatibility.COMMON_SUBSET;
        }
    }

    /**
     * Utility function that tests whether the commonContigs in both dicts are equivalent.  Equivalence means
     * that the seq records have the same length, if both are non-zero.
     *
     * @param commonContigs
     * @param dict1
     * @param dict2
     * @return true if all of the common contigs are equivalent
     */
    private static boolean commonContigsHaveSameLengths(final Set<String> commonContigs, final SAMSequenceDictionary dict1, final SAMSequenceDictionary dict2) {
        return findNotEqualCommonContigs(commonContigs, dict1, dict2) == null;
    }

    /**
     * Returns a List(x,y) that contains two sequence records that are not equal among the common contigs in both dicts.  Returns
     * null if all common contigs are equivalent
     *
     * @param commonContigs
     * @param dict1
     * @param dict2
     * @return
     */
    private static List<SAMSequenceRecord> findNotEqualCommonContigs(final Set<String> commonContigs, final SAMSequenceDictionary dict1, final SAMSequenceDictionary dict2) {
        for ( String name : commonContigs ) {
            SAMSequenceRecord elt1 = dict1.getSequence(name);
            SAMSequenceRecord elt2 = dict2.getSequence(name);
            if ( ! sequenceRecordsAreEquivalent(elt1, elt2) )
                return Arrays.asList(elt1,elt2);
        }

        return null;
    }

    /**
     * Helper routine that determines if two sequence records are equivalent, defined as having the same name,
     * lengths (if both are non-zero) and MD5 (if present)
     *
     * @param record1  a SAMSequenceRecord
     * @param record2  a SAMSequenceRecord
     * @return true if the records are equivalent, false otherwise
     */
    private static boolean sequenceRecordsAreEquivalent(final SAMSequenceRecord record1, final SAMSequenceRecord record2) {
        if ( record1 == record2 ) return true;
        if ( record1 == null || record2 == null ) return false;

        // compare length
        if ( record1.getSequenceLength() != 0 && record2.getSequenceLength() != 0 && record1.getSequenceLength() != record2.getSequenceLength() )
            return false;

        // compare name
        if ( !record1.getSequenceName().equals(record2.getSequenceName() ))
            return false;

         // compare MD5
         if ( record1.getMd5() != null && record2.getMd5() != null ){
            final BigInteger firstMd5 = new BigInteger(record1.getMd5(), 16);
            final BigInteger secondMd5 = new BigInteger(record2.getMd5(), 16);
            if ( !firstMd5.equals(secondMd5) )
                return false;
       }

        return true;
    }

    /**
     * A very simple (and naive) algorithm to determine (1) if the dict is a human reference (hg18/hg19) and if it's
     * lexicographically sorted.  Works by matching lengths of the static chr1, chr10, and chr2, and then if these
     * are all matched, requiring that the order be chr1, chr2, chr10.
     *
     * @param dict
     * @return
     */
    private static boolean nonCanonicalHumanContigOrder(final SAMSequenceDictionary dict) {
        if ( ! ENABLE_LEXICOGRAPHIC_REQUIREMENT_FOR_HUMAN ) // if we don't want to enable this test, just return false
            return false;

        SAMSequenceRecord chr1 = null, chr2 = null, chr10 = null;

        for ( final SAMSequenceRecord elt : dict.getSequences() ) {
            if ( isHumanSeqRecord(elt, CHR1_HG18, CHR1_HG19 ) ) chr1 = elt;
            if ( isHumanSeqRecord(elt, CHR2_HG18, CHR2_HG19 ) ) chr2 = elt;
            if ( isHumanSeqRecord(elt, CHR10_HG18, CHR10_HG19 ) ) chr10 = elt;
        }

        if ( chr1 != null && chr2 != null && chr10 != null) {
            // we found them all
            return ! ( chr1.getSequenceIndex() < chr2.getSequenceIndex() && chr2.getSequenceIndex() < chr10.getSequenceIndex() );
        } else {
            return false;
        }
    }

    /**
     * Trivial helper that returns true if elt has the same length as rec1 or rec2
     * @param elt record to test
     * @param rec1 first record to test for length equivalence
     * @param rec2 first record to test for length equivalence
     * @return true if elt has the same length as either rec1 or rec2
     */
    private static boolean isHumanSeqRecord(SAMSequenceRecord elt, SAMSequenceRecord rec1, SAMSequenceRecord rec2 ) {
        return elt.getSequenceLength() == rec1.getSequenceLength() || elt.getSequenceLength() == rec2.getSequenceLength();
    }

    /**
     * Returns true if the common contigs in dict1 and dict2 are in the same relative order, without regard to
     * absolute index position. This is accomplished by getting the common contigs in both dictionaries, sorting
     * these according to their indices, and then walking through the sorted list to ensure that each ordered contig
     * is equivalent
     *
     * @param commonContigs names of the contigs common to both dictionaries
     * @param dict1 first SAMSequenceDictionary
     * @param dict2 second SAMSequenceDictionary
     * @return true if the common contigs occur in the same relative order in both dict1 and dict2, otherwise false
     */
    private static boolean commonContigsAreInSameRelativeOrder(final Set<String> commonContigs, final SAMSequenceDictionary dict1, final SAMSequenceDictionary dict2) {
        List<SAMSequenceRecord> list1 = sortSequenceListByIndex(getSequencesOfName(commonContigs, dict1));
        List<SAMSequenceRecord> list2 = sortSequenceListByIndex(getSequencesOfName(commonContigs, dict2));

        for ( int i = 0; i < list1.size(); i++ ) {
            SAMSequenceRecord elt1 = list1.get(i);
            SAMSequenceRecord elt2 = list2.get(i);
            if ( ! elt1.getSequenceName().equals(elt2.getSequenceName()) )
                return false;
        }

        return true;
    }

    /**
     * Gets the subset of SAMSequenceRecords in commonContigs in dict
     *
     * @param commonContigs
     * @param dict
     * @return
     */
    private static List<SAMSequenceRecord> getSequencesOfName(final Set<String> commonContigs, final SAMSequenceDictionary dict) {
        final List<SAMSequenceRecord> l = new ArrayList<SAMSequenceRecord>(commonContigs.size());
        for ( String name : commonContigs ) {
            l.add(dict.getSequence(name) );
        }

        return l;
    }

    /**
     * Compares sequence records by their order
     */
    private static class CompareSequenceRecordsByIndex implements Comparator<SAMSequenceRecord> {
        public int compare(SAMSequenceRecord x, SAMSequenceRecord y) {
            return Integer.valueOf(x.getSequenceIndex()).compareTo(y.getSequenceIndex());
        }
    }

    /**
     * Returns a sorted list of SAMSequenceRecords sorted by their indices.  Note that the
     * list is modified in place, so the returned list is == to the unsorted list.
     *
     * @param unsorted
     * @return
     */
    private static List<SAMSequenceRecord> sortSequenceListByIndex(final List<SAMSequenceRecord> unsorted) {
        Collections.sort(unsorted, new CompareSequenceRecordsByIndex());
        return unsorted;
    }

    /**
     * Checks whether the common contigs in the given sequence dictionaries occur at the same indices
     * in both dictionaries
     *
     * @param commonContigs Set of names of the contigs that occur in both dictionaries
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return true if the contigs common to dict1 and dict2 occur at the same indices in both dictionaries,
     *         otherwise false
     */
    private static boolean commonContigsAreAtSameIndices( final Set<String> commonContigs, final SAMSequenceDictionary dict1, final SAMSequenceDictionary dict2 ) {
        for ( String commonContig : commonContigs ) {
            final SAMSequenceRecord dict1Record = dict1.getSequence(commonContig);
            final SAMSequenceRecord dict2Record = dict2.getSequence(commonContig);

            // Each common contig must have the same index in both dictionaries
            if ( dict1Record.getSequenceIndex() != dict2Record.getSequenceIndex() ) {
                return false;
            }
        }

        return true;
    }

    /**
     * Gets the set of names of the contigs found in both sequence dictionaries that have different indices
     * in the two dictionaries.
     *
     * @param commonContigs Set of names of the contigs common to both dictionaries
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return a Set containing the names of the common contigs indexed differently in dict1 vs. dict2,
     *         or an empty Set if there are no such contigs
     */
    private static Set<String> getDifferentlyIndexedCommonContigs( final Set<String> commonContigs,
                                                                   final SAMSequenceDictionary dict1,
                                                                   final SAMSequenceDictionary dict2 ) {

        final Set<String> differentlyIndexedCommonContigs = new LinkedHashSet<String>(Utils.optimumHashSize(commonContigs.size()));

        for ( String commonContig : commonContigs ) {
            if ( dict1.getSequence(commonContig).getSequenceIndex() != dict2.getSequence(commonContig).getSequenceIndex() ) {
                differentlyIndexedCommonContigs.add(commonContig);
            }
        }

        return differentlyIndexedCommonContigs;
    }

    /**
     * Finds the names of any contigs indexed differently in the two sequence dictionaries that also
     * occur in the provided set of intervals.
     *
     * @param intervals GenomeLocSortedSet containing the intervals to check
     * @param dict1 first sequence dictionary
     * @param dict2 second sequence dictionary
     * @return a Set of the names of the contigs indexed differently in dict1 vs dict2 that also
     *         occur in the provided intervals, or an empty Set if there are no such contigs
     */
    private static Set<String> findMisindexedContigsInIntervals( final GenomeLocSortedSet intervals,
                                                                 final SAMSequenceDictionary dict1,
                                                                 final SAMSequenceDictionary dict2 ) {

        final Set<String> differentlyIndexedCommonContigs = getDifferentlyIndexedCommonContigs(getCommonContigsByName(dict1, dict2), dict1, dict2);
        final Set<String> misindexedContigsInIntervals = new LinkedHashSet<String>(Utils.optimumHashSize(differentlyIndexedCommonContigs.size()));

        // We know differentlyIndexedCommonContigs is a HashSet, so this loop is O(intervals)
        for ( GenomeLoc interval : intervals ) {
            if ( differentlyIndexedCommonContigs.contains(interval.getContig()) ) {
                misindexedContigsInIntervals.add(interval.getContig());
            }
        }

        return misindexedContigsInIntervals;
    }

    /**
     * Returns the set of contig names found in both dicts.
     * @param dict1
     * @param dict2
     * @return
     */
    public static Set<String> getCommonContigsByName(SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        final Set<String> intersectingSequenceNames = getContigNames(dict1);
        intersectingSequenceNames.retainAll(getContigNames(dict2));
        return intersectingSequenceNames;
    }

    public static Set<String> getContigNames(SAMSequenceDictionary dict) {
        final Set<String> contigNames = new HashSet<String>(Utils.optimumHashSize(dict.size()));
        for (SAMSequenceRecord dictionaryEntry : dict.getSequences())
            contigNames.add(dictionaryEntry.getSequenceName());
        return contigNames;
    }

    /**
     * Returns a compact String representation of the sequence dictionary it's passed
     *
     * The format of the returned String is:
     * [ contig1Name(length: contig1Length) contig2Name(length: contig2Length) ... ]
     *
     * @param dict a non-null SAMSequenceDictionary
     * @return A String containing all of the contig names and lengths from the sequence dictionary it's passed
     */
    public static String getDictionaryAsString( final SAMSequenceDictionary dict ) {
        if ( dict == null ) {
            throw new IllegalArgumentException("Sequence dictionary must be non-null");
        }

        final StringBuilder s = new StringBuilder("[ ");

        for ( SAMSequenceRecord dictionaryEntry : dict.getSequences() ) {
            s.append(dictionaryEntry.getSequenceName());
            s.append("(length:");
            s.append(dictionaryEntry.getSequenceLength());
            s.append(") ");
        }

        s.append("]");

        return s.toString();
    }
}
