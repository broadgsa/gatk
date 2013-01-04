/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import net.sf.samtools.SAMUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.recalibration.EventType;
import org.broadinstitute.sting.utils.recalibration.ReadCovariates;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;

public final class ReadRecalibrationInfoUnitTest extends BaseTest {
    @DataProvider(name = "InfoProvider")
    public Object[][] createCombineTablesProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int readLength: Arrays.asList(10, 100, 1000) ) {
            for ( final boolean includeIndelErrors : Arrays.asList(true, false) ) {
                tests.add(new Object[]{readLength, includeIndelErrors});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "InfoProvider")
    public void testReadInfo(final int readLength, final boolean includeIndelErrors) {
        final ReadCovariates covariates = new ReadCovariates(readLength, 2);

        final byte[] bases = new byte[readLength];
        final byte[] baseQuals = new byte[readLength];
        final byte[] insertionQuals = new byte[readLength];
        final byte[] deletionQuals = new byte[readLength];
        final boolean[] skips = new boolean[readLength];
        final double[] snpErrors = new double[readLength];
        final double[] insertionErrors = new double[readLength];
        final double[] deletionsErrors = new double[readLength];
        for ( int i = 0; i < readLength; i++ ) {
            bases[i] = 'A';
            baseQuals[i] = (byte)(i % SAMUtils.MAX_PHRED_SCORE);
            insertionQuals[i] = (byte)((i+1) % SAMUtils.MAX_PHRED_SCORE);
            deletionQuals[i] = (byte)((i+2) % SAMUtils.MAX_PHRED_SCORE);
            skips[i] = i % 2 == 0;
            snpErrors[i] = 1.0 / (i+1);
            insertionErrors[i] = 0.5 / (i+1);
            deletionsErrors[i] = 0.3 / (i+1);
        }

        final EnumMap<EventType, double[]> errors = new EnumMap<EventType, double[]>(EventType.class);
        errors.put(EventType.BASE_SUBSTITUTION, snpErrors);
        errors.put(EventType.BASE_INSERTION, insertionErrors);
        errors.put(EventType.BASE_DELETION, deletionsErrors);

        final EnumMap<EventType, byte[]> quals = new EnumMap<EventType, byte[]>(EventType.class);
        quals.put(EventType.BASE_SUBSTITUTION, baseQuals);
        quals.put(EventType.BASE_INSERTION, insertionQuals);
        quals.put(EventType.BASE_DELETION, deletionQuals);

        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, baseQuals, readLength + "M");
        if ( includeIndelErrors ) {
            read.setBaseQualities(insertionQuals, EventType.BASE_INSERTION);
            read.setBaseQualities(deletionQuals, EventType.BASE_DELETION);
        }

        final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, covariates, skips, snpErrors, insertionErrors, deletionsErrors);

        Assert.assertEquals(info.getCovariatesValues(), covariates);
        Assert.assertEquals(info.getRead(), read);

        for ( int i = 0; i < readLength; i++ ) {
            Assert.assertEquals(info.skip(i), skips[i]);
            for ( final EventType et : EventType.values() ) {
                Assert.assertEquals(info.getErrorFraction(et, i), errors.get(et)[i]);
                final byte expectedQual = et == EventType.BASE_SUBSTITUTION || includeIndelErrors ? quals.get(et)[i]: GATKSAMRecord.DEFAULT_INSERTION_DELETION_QUAL;
                Assert.assertEquals(info.getQual(et, i), expectedQual);
            }
        }
    }
}
