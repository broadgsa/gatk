/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

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

    private static boolean ENABLE_LEXICOGRAPHIC_REQUIREMENT_FOR_HUMAN = true;

    // hg18
    private static SAMSequenceRecord CHR1_HG18 = new SAMSequenceRecord("chr1", 247249719);
    private static SAMSequenceRecord CHR2_HG18 = new SAMSequenceRecord("chr2", 242951149);
    private static SAMSequenceRecord CHR10_HG18 = new SAMSequenceRecord("chr10", 135374737);

    // hg19
    private static SAMSequenceRecord CHR1_HG19 = new SAMSequenceRecord("chr1", 249250621);
    private static SAMSequenceRecord CHR2_HG19 = new SAMSequenceRecord("chr2", 243199373);
    private static SAMSequenceRecord CHR10_HG19 = new SAMSequenceRecord("chr10", 135534747);

    public enum SequenceDictionaryCompatibility {
        IDENTICAL,                      // the dictionaries are identical
        COMMON_SUBSET,                  // there exists a common subset of equivalent contigs
        NO_COMMON_CONTIGS,              // no overlap between dictionaries
        UNEQUAL_COMMON_CONTIGS,         // common subset has contigs that have the same name but aren't equivalent
        NON_CANONICAL_HUMAN_ORDER,      // human reference detected but the order of the contigs is non-standard (lexicographic, for examine)
        OUT_OF_ORDER                    // the two dictionaries overlap but the contigs occur out of order w.r.t each other
    }

    /**
     * @param validationExclusion exclusions to validation
     * @return Returns true if the engine is in tolerant mode and we'll let through dangerous but not fatal dictionary inconsistency
     */
    public static boolean allowNonFatalIncompabilities(ValidationExclusion.TYPE validationExclusion) {
        return ( validationExclusion == ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY ||
                        validationExclusion == ValidationExclusion.TYPE.ALL );
    }

    /**
     * Testings for compatbility between dict1 and dict2.  If the dictionaries are incompatible, then UserExceptions are
     * thrown with detailed error messages.  If the engine is in permissive mode, then logger.warnings of generated instead
     *
     * @param logger for warnings
     * @param validationExclusion exclusions to validation
     * @param name1 name associated with dict1
     * @param dict1 the sequence dictionary dict1
     * @param name2 name associated with dict2
     * @param dict2 the sequence dictionary dict2
     */
    public static void validateDictionaries(Logger logger, ValidationExclusion.TYPE validationExclusion, String name1, SAMSequenceDictionary dict1, String name2, SAMSequenceDictionary dict2) {
        SequenceDictionaryCompatibility type = compareDictionaries(dict1, dict2);
        switch ( type ) {
            case IDENTICAL:
                return;
            case COMMON_SUBSET:
                 return;
            case NO_COMMON_CONTIGS:
                throw new UserException.IncompatibleSequenceDictionaries("No overlapping contigs found", name1, dict1, name2, dict2);

            case UNEQUAL_COMMON_CONTIGS: {
                List<SAMSequenceRecord> x = findDisequalCommonContigs(getCommonContigsByName(dict1, dict2), dict1, dict2);
                SAMSequenceRecord elt1 = x.get(0);
                SAMSequenceRecord elt2 = x.get(1);

                // todo -- replace with toString when SAMSequenceRecord has a nice toString routine
                UserException ex = new UserException.IncompatibleSequenceDictionaries(String.format("Found contigs with the same name but different lengths:\n  contig %s = %s / %d\n  contig %s = %s / %d",
                        name1, elt1.getSequenceName(), elt1.getSequenceLength(),
                        name2, elt2.getSequenceName(), elt2.getSequenceLength()),
                        name1, dict1, name2, dict2);

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
            }

            case OUT_OF_ORDER: {
                UserException ex = new UserException.IncompatibleSequenceDictionaries("Order of contigs differences, which is unsafe", name1, dict1, name2, dict2);
                if ( allowNonFatalIncompabilities(validationExclusion) )
                    logger.warn(ex.getMessage());
                else
                    throw ex;
            } break;
            default:
                throw new ReviewedStingException("Unexpected SequenceDictionaryComparison type: " + type);
        }
    }

    /**
     * Workhorse routine that takes two dictionaries and returns their compatibility.
     *
     * @param dict1
     * @param dict2
     * @return
     */
    public static SequenceDictionaryCompatibility compareDictionaries(SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        // If there's no overlap between reads and reference, data will be bogus.  Throw an exception.
        if ( nonCanonicalHumanContigOrder( dict1 ) || nonCanonicalHumanContigOrder(dict2) )
            return SequenceDictionaryCompatibility.NON_CANONICAL_HUMAN_ORDER;

        Set<String> commonContigs = getCommonContigsByName(dict1, dict2);

        if (commonContigs.size() == 0)
            return SequenceDictionaryCompatibility.NO_COMMON_CONTIGS;
        else if ( ! commonContigsAreEquivalent( commonContigs, dict1, dict2 ) )
            return SequenceDictionaryCompatibility.UNEQUAL_COMMON_CONTIGS;
        else if ( ! commonContigsAreInOrder( commonContigs, dict1, dict2 ) )
            return SequenceDictionaryCompatibility.OUT_OF_ORDER;
        else if ( commonContigs.size() == dict1.size() && commonContigs.size() == dict2.size() )
            return SequenceDictionaryCompatibility.IDENTICAL;
        else {
            return SequenceDictionaryCompatibility.COMMON_SUBSET;
        }
    }

    /**
     * Utility function that tests whether the commonContigs in both dicts are equivalent.  Equivalece means
     * that the seq records have the same length, if both are non-zero.
     *
     * @param commonContigs
     * @param dict1
     * @param dict2
     * @return true if all of the common contigs are equivalent
     */
    private static boolean commonContigsAreEquivalent(Set<String> commonContigs, SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        return findDisequalCommonContigs(commonContigs, dict1, dict2) == null;
    }

    /**
     * Returns a List(x,y) that contains two disequal sequence records among the common contigs in both dicts.  Returns
     * null if all common contigs are equivalent
     *
     * @param commonContigs
     * @param dict1
     * @param dict2
     * @return
     */
    private static List<SAMSequenceRecord> findDisequalCommonContigs(Set<String> commonContigs, SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        for ( String name : commonContigs ) {
            SAMSequenceRecord elt1 = dict1.getSequence(name);
            SAMSequenceRecord elt2 = dict2.getSequence(name);
            if ( ! SequenceRecordsAreEquivalent(elt1, elt2) )
                return Arrays.asList(elt1,elt2);
        }

        return null;
    }

    /**
     * Helper routine that returns two sequence records are equivalent, defined as having the same name and
     * lengths, if both are non-zero
     *
     * @param me
     * @param that
     * @return
     */
    private static boolean SequenceRecordsAreEquivalent(final SAMSequenceRecord me, final SAMSequenceRecord that) {
        if (me == that) return true;
        if (that == null) return false;

        // I don't care if the indices are difference
        //if (me.getSequenceIndex() != that.getSequenceIndex()) return false;
        if (me.getSequenceLength() != 0 && that.getSequenceLength() != 0 && me.getSequenceLength() != that.getSequenceLength())
            return false;

            // todo -- reenable if we want to be really strict here
//        if (me.getAttribute(SAMSequenceRecord.MD5_TAG) != null && that.getAttribute(SAMSequenceRecord.MD5_TAG) != null) {
//            final BigInteger thisMd5 = new BigInteger((String)me.getAttribute(SAMSequenceRecord.MD5_TAG), 16);
//            final BigInteger thatMd5 = new BigInteger((String)that.getAttribute(SAMSequenceRecord.MD5_TAG), 16);
//            if (!thisMd5.equals(thatMd5)) {
//                return false;
//            }
//        }
//        else {
        if (me.getSequenceName() != that.getSequenceName())
            return false; // Compare using == since we intern() the Strings
//        }

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
    private static boolean nonCanonicalHumanContigOrder(SAMSequenceDictionary dict) {
        if ( ! ENABLE_LEXICOGRAPHIC_REQUIREMENT_FOR_HUMAN ) // if we don't want to enable this test, just return false
            return false;

        SAMSequenceRecord chr1 = null, chr2 = null, chr10 = null;

        for ( SAMSequenceRecord elt : dict.getSequences() ) {
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
     * Returns true if the common contigs in dict1 and dict2 are in the same order.  This is accomplished by getting the
     * common contigs in both dictionaries, sorting these according to their indices, and the walking through
     * the sorted list to ensure that each ordered contig is equivalent
     *
     * @param commonContigs
     * @param dict1
     * @param dict2
     * @return
     */
    public static boolean commonContigsAreInOrder(Set<String> commonContigs, SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
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
    private static List<SAMSequenceRecord> getSequencesOfName(Set<String> commonContigs, SAMSequenceDictionary dict) {
        List<SAMSequenceRecord> l = new ArrayList<SAMSequenceRecord>(commonContigs.size());
        for ( String name : commonContigs ) {
            l.add(dict.getSequence(name) );
        }

        return l;
    }

    // --------------------------------------------------------------------------------------------------------------
    // Utilities for comparing the order of sequence records
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Compares sequence records by their order
     */
    private static class CompareSequenceRecordsByIndex implements Comparator<SAMSequenceRecord> {
        public int compare(SAMSequenceRecord x, SAMSequenceRecord y) {
            return new Integer(x.getSequenceIndex()).compareTo(y.getSequenceIndex());
        }
    }

    /**
     * Returns a sorted list of SAMSequenceRecords sorted by their indices.  Note that the
     * list is modified in place, so the returned list is == to the unsorted list.
     *
     * @param unsorted
     * @return
     */
    private static List<SAMSequenceRecord> sortSequenceListByIndex(List<SAMSequenceRecord> unsorted) {
        Collections.sort(unsorted, new CompareSequenceRecordsByIndex());
        return unsorted;
    }


    /**
     * Returns the set of contig names found in both dicts.
     * @param dict1
     * @param dict2
     * @return
     */
    public static Set<String> getCommonContigsByName(SAMSequenceDictionary dict1, SAMSequenceDictionary dict2) {
        Set<String> intersectingSequenceNames = getContigNames(dict1);
        intersectingSequenceNames.retainAll(getContigNames(dict2));
        return intersectingSequenceNames;
    }

    public static Set<String> getContigNames(SAMSequenceDictionary dict) {
        Set<String> contigNames = new HashSet<String>((int)(dict.size() / 0.75f) + 1, 0.75f);
        for (SAMSequenceRecord dictionaryEntry : dict.getSequences())
            contigNames.add(dictionaryEntry.getSequenceName());
        return contigNames;
    }
}
