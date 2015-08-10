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

package org.broadinstitute.gatk.engine.iterators;


import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.engine.iterators.MisencodedBaseQualityReadTransformer;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Basic unit test for misencoded quals
 */
public class MisencodedBaseQualityUnitTest extends BaseTest {

    private static final String readBases = "AAAAAAAAAA";
    private static final byte[] badQuals = { 59, 60, 62, 63, 64, 61, 62, 58, 57, 56 };
    private static final byte[] goodQuals = { 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 };
    private static final byte[] fixedQuals = { 28, 29, 31, 32, 33, 30, 31, 27, 26, 25 };
    private SAMFileHeader header;

    @BeforeMethod
    public void before() {
        // reset the read counter so that we are deterministic
        MisencodedBaseQualityReadTransformer.currentReadCounter = 0;
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
    }

    private GATKSAMRecord createRead(final boolean useGoodBases) {
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 10, readBases.getBytes(),
                                                                     useGoodBases ? Arrays.copyOf(goodQuals, goodQuals.length) :
                                                                                    Arrays.copyOf(badQuals, badQuals.length));
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