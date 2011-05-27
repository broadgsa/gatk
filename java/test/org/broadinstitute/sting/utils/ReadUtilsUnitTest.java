package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;


public class ReadUtilsUnitTest extends BaseTest {
    SAMRecord read;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?";

    @BeforeTest
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, BASES.length());
        read.setReadUnmappedFlag(true);
        read.setReadBases(new String(BASES).getBytes());
        read.setBaseQualityString(new String(QUALS));
    }

    private void testReadBasesAndQuals(SAMRecord read, int expectedStart, int expectedStop) {
        SAMRecord clipped = ReadUtils.hardClipBases(read, expectedStart, expectedStop - 1, null);
        String expectedBases = BASES.substring(expectedStart, expectedStop);
        String expectedQuals = QUALS.substring(expectedStart, expectedStop);
        Assert.assertEquals(clipped.getReadBases(), expectedBases.getBytes(), "Clipped bases not those expected");
        Assert.assertEquals(clipped.getBaseQualityString(), expectedQuals, "Clipped quals not those expected");
    }

    @Test public void testNoClip() { testReadBasesAndQuals(read, 0, 4); }
    @Test public void testClip1Front() { testReadBasesAndQuals(read, 1, 4); }
    @Test public void testClip2Front() { testReadBasesAndQuals(read, 2, 4); }
    @Test public void testClip1Back() { testReadBasesAndQuals(read, 0, 3); }
    @Test public void testClip2Back() { testReadBasesAndQuals(read, 0, 2); }
}
