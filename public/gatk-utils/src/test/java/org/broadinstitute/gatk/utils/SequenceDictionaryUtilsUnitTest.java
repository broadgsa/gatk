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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import static org.broadinstitute.gatk.utils.SequenceDictionaryUtils.*;
import static org.broadinstitute.gatk.utils.SequenceDictionaryUtils.SequenceDictionaryCompatibility.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SequenceDictionaryUtilsUnitTest extends BaseTest {

    private static Logger logger = Logger.getLogger(SequenceDictionaryUtilsUnitTest.class);


    @DataProvider( name = "SequenceDictionaryDataProvider" )
    public Object[][] generateSequenceDictionaryTestData() {
        final SAMSequenceRecord CHRM_HG19 = new SAMSequenceRecord("chrM", 16571);
        final SAMSequenceRecord CHR_NONSTANDARD1 = new SAMSequenceRecord("NonStandard1", 8675309);
        final SAMSequenceRecord CHR_NONSTANDARD2 = new SAMSequenceRecord("NonStandard2", 8675308);

        final Class NO_COMMON_CONTIGS_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;
        final Class UNEQUAL_COMMON_CONTIGS_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;
        final Class NON_CANONICAL_HUMAN_ORDER_EXCEPTION = UserException.LexicographicallySortedSequenceDictionary.class;
        final Class OUT_OF_ORDER_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;
        final Class DIFFERENT_INDICES_EXCEPTION = UserException.IncompatibleSequenceDictionaries.class;

        final List<SAMSequenceRecord> hg19Sequences = Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19, CHR10_HG19);
        final GenomeLocParser hg19GenomeLocParser = new GenomeLocParser(new SAMSequenceDictionary(hg19Sequences));
        final List<GenomeLoc> hg19AllContigsIntervals = Arrays.asList(hg19GenomeLocParser.createGenomeLoc("chrM", 0, 1),
                                                                      hg19GenomeLocParser.createGenomeLoc("chr1", 0, 1),
                                                                      hg19GenomeLocParser.createGenomeLoc("chr2", 0, 1),
                                                                      hg19GenomeLocParser.createGenomeLoc("chr10", 0, 1));
        final List<GenomeLoc> hg19PartialContigsIntervals = Arrays.asList(hg19GenomeLocParser.createGenomeLoc("chrM", 0, 1),
                                                                          hg19GenomeLocParser.createGenomeLoc("chr1", 0, 1));
        final GenomeLocSortedSet hg19AllContigsIntervalSet = new GenomeLocSortedSet(hg19GenomeLocParser, hg19AllContigsIntervals);
        final GenomeLocSortedSet hg19PartialContigsIntervalSet = new GenomeLocSortedSet(hg19GenomeLocParser, hg19PartialContigsIntervals);

        return new Object[][]  {
            // Identical dictionaries:
            { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG19),                        null, IDENTICAL, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), null, IDENTICAL, null, false, null },
            { Arrays.asList(CHR1_B37),                         Arrays.asList(CHR1_B37),                         null, IDENTICAL, null, false, null },
            { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),    Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),    null, IDENTICAL, null, false, null },

            // Dictionaries with a common subset:
            { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                                   null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1),                        Arrays.asList(CHR1_HG19, CHR_NONSTANDARD2),                                   null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19),                        Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19),                                   null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19),                        Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHRM_HG19),                        null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD2),            null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                              null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1),            null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHR2_HG19, CHR10_HG19),            null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG19, CHR2_HG19, CHR10_HG19, CHRM_HG19), null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37, CHR_NONSTANDARD1),    Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37, CHR_NONSTANDARD2),               null, COMMON_SUBSET, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                              null, COMMON_SUBSET, null, false, null },

            // Dictionaries with no common contigs:
            { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                     null, NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_B37),                      null, NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37), null, NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),             Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37), null, NO_COMMON_CONTIGS, NO_COMMON_CONTIGS_EXCEPTION, false, null },

            // Dictionaries with unequal common contigs:
            { Arrays.asList(CHR1_HG19),                                          Arrays.asList(CHR1_HG18),                                          null, UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_B36),                                           Arrays.asList(CHR1_B37),                                           null, UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19),                   Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18),                   null, UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_B37, CHR2_B37, CHR10_B37),                      Arrays.asList(CHR1_B36, CHR2_B36, CHR10_B36),                      null, UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19, CHR_NONSTANDARD1), Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18, CHR_NONSTANDARD2), null, UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR_NONSTANDARD2, CHR1_HG18, CHR2_HG18, CHR10_HG18), null, UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),                               Arrays.asList(CHR1_HG18, CHR2_HG18, CHR10_HG18),                   null, UNEQUAL_COMMON_CONTIGS, UNEQUAL_COMMON_CONTIGS_EXCEPTION, false, null },

            // One or both dictionaries in non-canonical human order:
            { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), null, NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG18, CHR10_HG18, CHR2_HG18), Arrays.asList(CHR1_HG18, CHR10_HG18, CHR2_HG18), null, NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), null, NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), null, NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHR1_B37, CHR10_B37, CHR2_B37),    Arrays.asList(CHR1_B37, CHR10_B37, CHR2_B37),    null, NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHR1_B36, CHR10_B36, CHR2_B36),    Arrays.asList(CHR1_B36, CHR10_B36, CHR2_B36),    null, NON_CANONICAL_HUMAN_ORDER, NON_CANONICAL_HUMAN_ORDER_EXCEPTION, false, null },

            // Dictionaries with a common subset, but different relative ordering within that subset:
            { Arrays.asList(CHR1_HG19, CHR2_HG19),            Arrays.asList(CHR2_HG19, CHR1_HG19),                              null, OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHR1_HG19, CHRM_HG19),                   null, OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHRM_HG19, CHR2_HG19, CHR1_HG19),                   null, OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19), Arrays.asList(CHR2_HG19, CHRM_HG19, CHR1_HG19),                   null, OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, null },
            { Arrays.asList(CHR1_B37, CHR2_B37),              Arrays.asList(CHR2_B37, CHR1_B37),                                null, OUT_OF_ORDER, OUT_OF_ORDER_EXCEPTION, false, null },


            // Dictionaries with a common subset in the same relative order, but with different indices.
            // This will only throw an exception during validation if isReadsToReferenceComparison is true,
            // and there are intervals overlapping the misindexed contigs:

            // These have isReadsToReferenceComparison == true and overlapping intervals, so we expect an exception:
            { Arrays.asList(CHRM_HG19, CHR1_HG19),                                                 Arrays.asList(CHR1_HG19),                                          null, DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true, hg19AllContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19),                    null, DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true, hg19AllContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1),  null, DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true, hg19AllContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),                                                 Arrays.asList(CHRM_HG19, CHR_NONSTANDARD1, CHR1_HG19, CHR2_HG19),  null, DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true, hg19AllContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19 ),                   Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19),                    null, DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true, hg19AllContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR2_HG19, CHR_NONSTANDARD1, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), null, DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true, hg19AllContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19, CHR_NONSTANDARD2 ), null, DIFFERENT_INDICES, DIFFERENT_INDICES_EXCEPTION, true, hg19AllContigsIntervalSet },

            // These have isReadsToReferenceComparison == true but no overlapping intervals, so we don't expect an exception:
            { Arrays.asList(CHR2_HG19, CHR10_HG19),                              Arrays.asList(CHR10_HG19),                       null, DIFFERENT_INDICES, null, true, hg19PartialContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19),             Arrays.asList(CHR1_HG19, CHR2_HG19),             null, DIFFERENT_INDICES, null, true, hg19PartialContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHR10_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHR10_HG19), null, DIFFERENT_INDICES, null, true, hg19PartialContigsIntervalSet },

            // These have isReadsToReferenceComparison == false, so we don't expect an exception:
            { Arrays.asList(CHRM_HG19, CHR1_HG19),                              Arrays.asList(CHR1_HG19),                       null, DIFFERENT_INDICES, null, false, hg19AllContigsIntervalSet },
            { Arrays.asList(CHR1_HG19, CHR_NONSTANDARD1, CHR2_HG19, CHRM_HG19), Arrays.asList(CHR1_HG19, CHR2_HG19, CHRM_HG19), null, DIFFERENT_INDICES, null, false, hg19AllContigsIntervalSet },


            // Tests for validation exclusions. Note that errors resulting from NO_COMMON_CONTIGs cannot be suppressed
            { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                        ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY, NO_COMMON_CONTIGS,         NO_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR2_HG19),                        ValidationExclusion.TYPE.ALL,                            NO_COMMON_CONTIGS,         NO_COMMON_CONTIGS_EXCEPTION, false, null },
            { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG18),                        ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY, UNEQUAL_COMMON_CONTIGS,    null, false, null },
            { Arrays.asList(CHR1_HG19),                        Arrays.asList(CHR1_HG18),                        ValidationExclusion.TYPE.ALL,                            UNEQUAL_COMMON_CONTIGS,    null, false, null },
            { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY, NON_CANONICAL_HUMAN_ORDER, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), Arrays.asList(CHR1_HG19, CHR10_HG19, CHR2_HG19), ValidationExclusion.TYPE.ALL,                            NON_CANONICAL_HUMAN_ORDER, null, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),             Arrays.asList(CHR2_HG19, CHR1_HG19),             ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY, OUT_OF_ORDER,              null, false, null },
            { Arrays.asList(CHR1_HG19, CHR2_HG19),             Arrays.asList(CHR2_HG19, CHR1_HG19),             ValidationExclusion.TYPE.ALL,                            OUT_OF_ORDER,              null, false, null },
            { Arrays.asList(CHRM_HG19, CHR1_HG19),             Arrays.asList(CHR1_HG19),                        ValidationExclusion.TYPE.ALLOW_SEQ_DICT_INCOMPATIBILITY, DIFFERENT_INDICES,         null, true, hg19AllContigsIntervalSet },
            { Arrays.asList(CHRM_HG19, CHR1_HG19),             Arrays.asList(CHR1_HG19),                        ValidationExclusion.TYPE.ALL,                            DIFFERENT_INDICES,         null, true, hg19AllContigsIntervalSet }
        };
    }

    @Test( dataProvider = "SequenceDictionaryDataProvider" )
    public void testSequenceDictionaryValidation( final List<SAMSequenceRecord> firstDictionaryContigs,
                                                  final List<SAMSequenceRecord> secondDictionaryContigs,
                                                  final ValidationExclusion.TYPE validationExclusions,
                                                  final SequenceDictionaryUtils.SequenceDictionaryCompatibility dictionaryCompatibility,
                                                  final Class expectedExceptionUponValidation,
                                                  final boolean isReadsToReferenceComparison,
                                                  final GenomeLocSortedSet intervals ) {

        final SAMSequenceDictionary firstDictionary = createSequenceDictionary(firstDictionaryContigs);
        final SAMSequenceDictionary secondDictionary = createSequenceDictionary(secondDictionaryContigs);
        final String testDescription = String.format("First dictionary: %s  Second dictionary: %s  Validation exclusions: %s",
                                                     SequenceDictionaryUtils.getDictionaryAsString(firstDictionary),
                                                     SequenceDictionaryUtils.getDictionaryAsString(secondDictionary),
                                                     validationExclusions);

        Exception exceptionThrown = null;
        try {
            SequenceDictionaryUtils.validateDictionaries(logger,
                                                         validationExclusions,
                                                         "firstDictionary",
                                                         firstDictionary,
                                                         "secondDictionary",
                                                         secondDictionary,
                                                         isReadsToReferenceComparison,
                                                         intervals);
        }
        catch ( Exception e ) {
            exceptionThrown = e;
        }

        if ( expectedExceptionUponValidation != null ) {
            Assert.assertTrue(exceptionThrown != null && expectedExceptionUponValidation.isInstance(exceptionThrown),
                              String.format("Expected exception %s but saw %s instead. %s",
                                            expectedExceptionUponValidation.getSimpleName(),
                                            exceptionThrown == null ? "no exception" : exceptionThrown.getClass().getSimpleName(),
                                            testDescription));
        }
        else {
            Assert.assertTrue(exceptionThrown == null,
                              String.format("Expected no exception but saw exception %s instead. %s",
                                            exceptionThrown != null ? exceptionThrown.getClass().getSimpleName() : "none",
                                            testDescription));
        }
    }

    @Test( dataProvider = "SequenceDictionaryDataProvider" )
    public void testSequenceDictionaryComparison( final List<SAMSequenceRecord> firstDictionaryContigs,
                                                  final List<SAMSequenceRecord> secondDictionaryContigs,
                                                  final ValidationExclusion.TYPE validationExclusions,
                                                  final SequenceDictionaryUtils.SequenceDictionaryCompatibility dictionaryCompatibility,
                                                  final Class expectedExceptionUponValidation,
                                                  final boolean isReadsToReferenceComparison,
                                                  final GenomeLocSortedSet intervals ) {

        final SAMSequenceDictionary firstDictionary = createSequenceDictionary(firstDictionaryContigs);
        final SAMSequenceDictionary secondDictionary = createSequenceDictionary(secondDictionaryContigs);
        final String testDescription = String.format("First dictionary: %s  Second dictionary: %s",
                                                     SequenceDictionaryUtils.getDictionaryAsString(firstDictionary),
                                                     SequenceDictionaryUtils.getDictionaryAsString(secondDictionary));

        final SequenceDictionaryUtils.SequenceDictionaryCompatibility reportedCompatibility =
              SequenceDictionaryUtils.compareDictionaries(firstDictionary, secondDictionary);

        Assert.assertTrue(reportedCompatibility == dictionaryCompatibility,
                          String.format("Dictionary comparison should have returned %s but instead returned %s. %s",
                                        dictionaryCompatibility, reportedCompatibility, testDescription));
    }

    private SAMSequenceDictionary createSequenceDictionary( final List<SAMSequenceRecord> contigs ) {
        final List<SAMSequenceRecord> clonedContigs = new ArrayList<SAMSequenceRecord>(contigs.size());

        // Clone the individual SAMSequenceRecords to avoid contig-index issues with shared objects
        // across multiple dictionaries in tests
        for ( SAMSequenceRecord contig : contigs ) {
            clonedContigs.add(contig.clone());
        }

        return new SAMSequenceDictionary(clonedContigs);
    }
}
