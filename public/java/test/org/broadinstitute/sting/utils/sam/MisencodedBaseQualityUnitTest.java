package org.broadinstitute.sting.utils.sam;


import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Basic unit test for misencoded quals
 */
public class MisencodedBaseQualityUnitTest extends BaseTest {

    private static final String readBases = "AAAAAAAAAA";
    private static final byte[] badQuals = { 'Z', '[', 'c', 'd', 'e', 'a', 'b', 'Z', 'Y', 'X' };
    private static final byte[] goodQuals = { '[', '[', '[', '[', '[', '[', '[', '[', '[', '[' };
    private static final byte[] fixedQuals = { ';', '<', 'D', 'E', 'F', 'B', 'C', ';', ':', '9' };
    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
    }

    private GATKSAMRecord createRead(final boolean useGoodBases) {
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 10, readBases.getBytes(), useGoodBases ? goodQuals : badQuals);
        read.setCigarString("10M");
        return read;
    }

    @Test(enabled = true)
    public void testGoodQuals() {
        final List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>(10000);
        for ( int i = 0; i < 10000; i++ )
            reads.add(createRead(true));

        testEncoding(reads);
    }

    @Test(enabled = true, expectedExceptions = {UserException.class})
    public void testBadQualsThrowsError() {
        final List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>(10000);
        for ( int i = 0; i < 10000; i++ )
            reads.add(createRead(false));

        testEncoding(reads);
    }

    @Test(enabled = true)
    public void testFixBadQuals() {
        final GATKSAMRecord read = createRead(false);
        final GATKSAMRecord fixedRead = MisencodedBaseQualityReadTransformer.fixMisencodedQuals(read);
        for ( int i = 0; i < fixedQuals.length; i++ )
            Assert.assertEquals(fixedQuals[i], fixedRead.getBaseQualities()[i]);
    }

    private void testEncoding(final List<GATKSAMRecord> reads) {
        for ( final GATKSAMRecord read : reads )
            MisencodedBaseQualityReadTransformer.checkForMisencodedQuals(read);
    }
}