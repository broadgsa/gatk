package org.broadinstitute.sting.gatk.filters;

import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Checks that the Bad Cigar filter works for all kinds of wonky cigars
 *
 * @author Mauricio Carneiro
 * @since 3/20/12
 */
public class BadCigarFilterUnitTest {

    BadCigarFilter filter;

    @BeforeClass
    public void init() {
        filter = new BadCigarFilter();
    }

    @Test
    public void testWonkyCigars () {
        byte[] bases = {'A', 'A', 'A', 'A'};
        byte[] quals = {30, 30, 30, 30};
        GATKSAMRecord read;
                                                                                                                        // starting with multiple deletions
        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "2D4M");
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "4M2D");                                           // ending with multiple deletions
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "3M1I1D");                                         // adjacent indels AND ends in deletion
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "1M1I1D2M");                                       // adjacent indels I->D
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "1M1D2I1M");                                       // adjacent indels D->I
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "1M1I2M1D");                                       // ends in single deletion with insertion in the middle
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "4M1D");                                           // ends in single deletion
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "1D4M");                                           // starts with single deletion
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "2M1D1D2M");                                       // adjacent D's
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());

        read = ArtificialSAMUtils.createArtificialRead(bases, quals, "1M1I1I1M");                                       // adjacent I's
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());
    }
}
