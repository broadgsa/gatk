// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.


import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;

/**
 * Basic unit test for GenomeLoc
 */
public class GenomeLocUnitTest extends BaseTest {
    private static ReferenceSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    /**
     * Tests that we got a string parameter in correctly
     */
    @Test
    public void testIsBetween() {
        logger.warn("Executing testIsBetween");

        GenomeLoc locMiddle = genomeLocParser.createGenomeLoc("chr1", 3, 3);

        GenomeLoc locLeft = genomeLocParser.createGenomeLoc("chr1", 1, 1);
        GenomeLoc locRight = genomeLocParser.createGenomeLoc("chr1", 5, 5);

        Assert.assertTrue(locMiddle.isBetween(locLeft, locRight));
        Assert.assertFalse(locLeft.isBetween(locMiddle, locRight));
        Assert.assertFalse(locRight.isBetween(locLeft, locMiddle));

    }
    @Test
    public void testContigIndex() {
        logger.warn("Executing testContigIndex");
        GenomeLoc locOne = genomeLocParser.createGenomeLoc("chr1",1,1);
        Assert.assertEquals(1, locOne.getContigIndex());
        Assert.assertEquals("chr1", locOne.getContig());

        GenomeLoc locX = genomeLocParser.createGenomeLoc("chrX",1,1);
        Assert.assertEquals(23, locX.getContigIndex());
        Assert.assertEquals("chrX", locX.getContig());

        GenomeLoc locNumber = genomeLocParser.createGenomeLoc(seq.getSequenceDictionary().getSequence(1).getSequenceName(),1,1);
        Assert.assertEquals(1, locNumber.getContigIndex());
        Assert.assertEquals("chr1", locNumber.getContig());
        Assert.assertEquals(0, locOne.compareTo(locNumber));

    }

    @Test
    public void testCompareTo() {
        logger.warn("Executing testCompareTo");
        GenomeLoc twoOne = genomeLocParser.createGenomeLoc("chr2", 1);
        GenomeLoc twoFive = genomeLocParser.createGenomeLoc("chr2", 5);
        GenomeLoc twoOtherFive = genomeLocParser.createGenomeLoc("chr2", 5);
        Assert.assertEquals(twoFive.compareTo(twoOtherFive), 0);

        Assert.assertEquals(twoOne.compareTo(twoFive), -1);
        Assert.assertEquals(twoFive.compareTo(twoOne), 1);

        GenomeLoc oneOne = genomeLocParser.createGenomeLoc("chr1", 5);
        Assert.assertEquals(oneOne.compareTo(twoOne), -1);
        Assert.assertEquals(twoOne.compareTo(oneOne), 1);
    }


    @Test
    public void testUnmappedSort() {
        GenomeLoc chr1 = genomeLocParser.createGenomeLoc("chr1",1,10000000);
        GenomeLoc chr2 = genomeLocParser.createGenomeLoc("chr2",1,10000000);
        GenomeLoc unmapped = GenomeLoc.UNMAPPED;

        List<GenomeLoc> unmappedOnly = Arrays.asList(unmapped);
        Collections.sort(unmappedOnly);
        Assert.assertEquals(unmappedOnly.size(),1,"Wrong number of elements in unmapped-only list.");
        Assert.assertEquals(unmappedOnly.get(0),unmapped,"List sorted in wrong order");

        List<GenomeLoc> chr1Presorted = Arrays.asList(chr1,unmapped);
        Collections.sort(chr1Presorted);
        Assert.assertEquals(chr1Presorted.size(),2,"Wrong number of elements in chr1,unmapped list.");
        Assert.assertEquals(chr1Presorted,Arrays.asList(chr1,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1Inverted = Arrays.asList(unmapped,chr1);
        Collections.sort(chr1Inverted);
        Assert.assertEquals(chr1Inverted.size(),2,"Wrong number of elements in chr1,unmapped list.");
        Assert.assertEquals(chr1Inverted,Arrays.asList(chr1,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1and2Presorted = Arrays.asList(chr1,chr2,unmapped);
        Collections.sort(chr1and2Presorted);
        Assert.assertEquals(chr1and2Presorted.size(),3,"Wrong number of elements in chr1,chr2,unmapped list.");
        Assert.assertEquals(chr1and2Presorted,Arrays.asList(chr1,chr2,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1and2UnmappedInFront = Arrays.asList(unmapped,chr1,chr2);
        Collections.sort(chr1and2UnmappedInFront);
        Assert.assertEquals(chr1and2UnmappedInFront.size(),3,"Wrong number of elements in unmapped,chr1,chr2 list.");
        Assert.assertEquals(chr1and2UnmappedInFront,Arrays.asList(chr1,chr2,unmapped),"List sorted in wrong order");

        List<GenomeLoc> chr1and2UnmappedSandwiched = Arrays.asList(chr1,unmapped,chr2);
        Collections.sort(chr1and2UnmappedSandwiched);
        Assert.assertEquals(chr1and2UnmappedSandwiched.size(),3,"Wrong number of elements in chr1,unmapped,chr2 list.");
        Assert.assertEquals(chr1and2UnmappedSandwiched,Arrays.asList(chr1,chr2,unmapped),"List sorted in wrong order");
    }

    @Test
    public void testUnmappedMerge() {
        GenomeLoc chr1 = genomeLocParser.createGenomeLoc("chr1",1,10000000);
        GenomeLoc unmapped = GenomeLoc.UNMAPPED;

        List<GenomeLoc> oneUnmappedOnly = Arrays.asList(unmapped);
        oneUnmappedOnly = IntervalUtils.sortAndMergeIntervals(genomeLocParser,oneUnmappedOnly, IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(oneUnmappedOnly.size(),1,"Wrong number of elements in list.");
        Assert.assertEquals(oneUnmappedOnly.get(0),unmapped,"List sorted in wrong order");

        List<GenomeLoc> twoUnmapped = Arrays.asList(unmapped,unmapped);
        twoUnmapped = IntervalUtils.sortAndMergeIntervals(genomeLocParser,twoUnmapped,IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(twoUnmapped.size(),1,"Wrong number of elements in list.");
        Assert.assertEquals(twoUnmapped.get(0),unmapped,"List sorted in wrong order");

        List<GenomeLoc> twoUnmappedAtEnd = Arrays.asList(chr1,unmapped,unmapped);
        twoUnmappedAtEnd = IntervalUtils.sortAndMergeIntervals(genomeLocParser,twoUnmappedAtEnd,IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(twoUnmappedAtEnd.size(),2,"Wrong number of elements in list.");
        Assert.assertEquals(twoUnmappedAtEnd,Arrays.asList(chr1,unmapped),"List sorted in wrong order");

        List<GenomeLoc> twoUnmappedMixed = Arrays.asList(unmapped,chr1,unmapped);
        twoUnmappedMixed = IntervalUtils.sortAndMergeIntervals(genomeLocParser,twoUnmappedMixed,IntervalMergingRule.OVERLAPPING_ONLY).toList();
        Assert.assertEquals(twoUnmappedMixed.size(),2,"Wrong number of elements in list.");
        Assert.assertEquals(twoUnmappedMixed,Arrays.asList(chr1,unmapped),"List sorted in wrong order");
    }

    // -------------------------------------------------------------------------------------
    //
    // testing overlap detection
    //
    // -------------------------------------------------------------------------------------

    private class ReciprocalOverlapProvider extends TestDataProvider {
        GenomeLoc gl1, gl2;
        int overlapSize;
        double overlapFraction;

        private ReciprocalOverlapProvider(int start1, int stop1, int start2, int stop2) {
            super(ReciprocalOverlapProvider.class);
            gl1 = genomeLocParser.createGenomeLoc("chr1", start1, stop1);
            gl2 = genomeLocParser.createGenomeLoc("chr1", start2, stop2);

            int shared = 0;
            for ( int i = start1; i <= stop1; i++ ) {
                if ( i >= start2 && i <= stop2 )
                    shared++;
            }

            this.overlapSize = shared;
            this.overlapFraction = Math.min((1.0*shared)/gl1.size(), (1.0*shared)/gl2.size());
            super.setName(String.format("%d-%d / %d-%d overlap=%d / %.2f", start1, stop1, start2, stop2, overlapSize, overlapFraction));
        }
    }

    @DataProvider(name = "ReciprocalOverlapProvider")
    public Object[][] makeReciprocalOverlapProvider() {
        for ( int start1 = 1; start1 <= 10; start1++ ) {
            for ( int stop1 = start1; stop1 <= 10; stop1++ ) {
                new ReciprocalOverlapProvider(start1, stop1, 1, 10);
                new ReciprocalOverlapProvider(start1, stop1, 5, 10);
                new ReciprocalOverlapProvider(start1, stop1, 5, 7);
                new ReciprocalOverlapProvider(start1, stop1, 5, 15);
                new ReciprocalOverlapProvider(start1, stop1, 11, 20);

                new ReciprocalOverlapProvider(1, 10, start1, stop1);
                new ReciprocalOverlapProvider(5, 10, start1, stop1);
                new ReciprocalOverlapProvider(5, 7, start1, stop1);
                new ReciprocalOverlapProvider(5, 15, start1, stop1);
                new ReciprocalOverlapProvider(11, 20, start1, stop1);
            }
        }

        return ReciprocalOverlapProvider.getTests(ReciprocalOverlapProvider.class);
    }

    @Test(dataProvider = "ReciprocalOverlapProvider")
    public void testReciprocalOverlapProvider(ReciprocalOverlapProvider cfg) {
        if ( cfg.overlapSize == 0 ) {
            Assert.assertFalse(cfg.gl1.overlapsP(cfg.gl2));
        } else {
            Assert.assertTrue(cfg.gl1.overlapsP(cfg.gl2));
            Assert.assertEquals(cfg.gl1.intersect(cfg.gl2).size(), cfg.overlapSize);
            Assert.assertEquals(cfg.gl1.reciprocialOverlapFraction(cfg.gl2), cfg.overlapFraction);
        }
    }

    // -------------------------------------------------------------------------------------
    //
    // testing comparison, hashcode, and equals
    //
    // -------------------------------------------------------------------------------------

    @DataProvider(name = "GenomeLocComparisons")
    public Object[][] createGenomeLocComparisons() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final int start = 10;
        for ( int stop = start; stop < start + 3; stop++ ) {
            final GenomeLoc g1 = genomeLocParser.createGenomeLoc("chr2", start, stop);
            for ( final String contig : Arrays.asList("chr1", "chr2", "chr3")) {
                for ( int start2 = start - 1; start2 <= stop + 1; start2++ ) {
                    for ( int stop2 = start2; stop2 < stop + 2; stop2++ ) {
                        final GenomeLoc g2 = genomeLocParser.createGenomeLoc(contig, start2, stop2);

                        ComparisonResult cmp = ComparisonResult.EQUALS;
                        if ( contig.equals("chr3") ) cmp = ComparisonResult.LESS_THAN;
                        else if ( contig.equals("chr1") ) cmp = ComparisonResult.GREATER_THAN;
                        else if ( start < start2 ) cmp = ComparisonResult.LESS_THAN;
                        else if ( start > start2 ) cmp = ComparisonResult.GREATER_THAN;
                        else if ( stop < stop2 ) cmp = ComparisonResult.LESS_THAN;
                        else if ( stop > stop2 ) cmp = ComparisonResult.GREATER_THAN;

                        tests.add(new Object[]{g1, g2, cmp});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private enum ComparisonResult {
        LESS_THAN(-1),
        EQUALS(0),
        GREATER_THAN(1);

        final int cmp;

        private ComparisonResult(int cmp) {
            this.cmp = cmp;
        }
    }

    @Test(dataProvider = "GenomeLocComparisons")
    public void testGenomeLocComparisons(GenomeLoc g1, GenomeLoc g2, ComparisonResult expected) {
        Assert.assertEquals(g1.compareTo(g2), expected.cmp, "Comparing genome locs failed");
        Assert.assertEquals(g1.equals(g2), expected == ComparisonResult.EQUALS);
        if ( expected == ComparisonResult.EQUALS )
            Assert.assertEquals(g1.hashCode(), g2.hashCode(), "Equal genome locs don't have the same hash code");
    }
}
