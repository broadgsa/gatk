/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import com.google.java.contract.Requires;
import org.broadinstitute.gatk.engine.recalibration.covariates.Covariate;
import org.broadinstitute.gatk.engine.recalibration.covariates.RepeatLengthCovariate;
import org.broadinstitute.gatk.engine.recalibration.covariates.RepeatUnitAndLengthCovariate;
import org.broadinstitute.gatk.engine.recalibration.covariates.RepeatUnitCovariate;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class RepeatCovariatesUnitTest {

    RepeatLengthCovariate rlCovariate;
    RepeatUnitCovariate ruCovariate;
    RepeatUnitAndLengthCovariate rurlCovariate;
    RecalibrationArgumentCollection RAC;



    @BeforeClass
    public void init() {
        RAC = new RecalibrationArgumentCollection();
        rlCovariate = new RepeatLengthCovariate();
        ruCovariate = new RepeatUnitCovariate();
        rurlCovariate = new RepeatUnitAndLengthCovariate();
        rlCovariate.initialize(RAC);
        ruCovariate.initialize(RAC);
        rurlCovariate.initialize(RAC);
    }

    @BeforeMethod
    public void initCache() {
        ReadCovariates.clearKeysCache();
    }


    @Test
    public void testFindNumberOfRepetitions() {
        // First, test logic to compute number of repetitions of a substring on a given string.
        int result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACAC".getBytes(), true);
        Assert.assertEquals(2,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACAC".getBytes(), true);
        Assert.assertEquals(4,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACACGT".getBytes(), true);
        Assert.assertEquals(4,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "GTACACACAC".getBytes(), true);
        Assert.assertEquals(0,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGT".getBytes(), true);
        Assert.assertEquals(0,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GCAGCAGTAGGGTGTACACACAC".getBytes(), true);
        Assert.assertEquals(1,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GTAGGGTGTACACACACGCAGCAT".getBytes(), true);
        Assert.assertEquals(0,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGTGTACACACACGCAGCAGCA".getBytes(), true);
        Assert.assertEquals(0,result);
        // Same tests but looking backward on string
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACAC".getBytes(), false);
        Assert.assertEquals(2,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACAC".getBytes(), false);
        Assert.assertEquals(4,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "ACACACACGT".getBytes(), false);
        Assert.assertEquals(0,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("AC".getBytes(), "GTACACACAC".getBytes(), false);
        Assert.assertEquals(4,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGT".getBytes(), false);
        Assert.assertEquals(0,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GCAGCAGTAGGGTGTACACACAC".getBytes(), false);
        Assert.assertEquals(0,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCAGCA".getBytes(), "GTAGGGTGTACACACACGCAGCAT".getBytes(), false);
        Assert.assertEquals(0,result);
        result = GATKVariantContextUtils.findNumberOfRepetitions("GCA".getBytes(), "GTAGGGTGTACACACACGCAGCAGCA".getBytes(), false);
        Assert.assertEquals(3,result);

        // test logic to get repeat unit and number of repeats from covariate value
        final String[] repUnits = new String[]{"AG","CCG","TCCA","T"};
        for (String ru : repUnits) {
            for (int k=1; k < 10; k++) {
                Pair<String,Integer> pair = RepeatLengthCovariate.getRUandNRfromCovariate(String.format("%s%d",ru,k));
                Assert.assertEquals(pair.second.intValue(),k);
                Assert.assertEquals(pair.first,ru);
            }
        }

    }

    /**
     * Build synthetic reads with random content made up of tandem repeats, record computed Repeat Unit and # repeats and see if
     * they match with read context
     */
    @Test
    public void testManyObservations() {
        final int NUM_UNITS = 10;
        final int MAX_REPEAT_UNIT_LENGTH = RAC.MAX_STR_UNIT_LENGTH;
        final int MAX_NUM_REPETITIONS = RAC.MAX_REPEAT_LENGTH;
        final int NUM_TEST_CASES = 100;

        Random random = new Random();

        for (int r = 0; r < NUM_TEST_CASES; r++) {
            final StringBuilder sb = new StringBuilder();
            // for each unit, generate a repeat unit at random with given random length
            final ArrayList<String> repeatUnits = new ArrayList<String>();
            final ArrayList<Integer> numsRepetitions = new ArrayList<Integer>();
            for (int n=0; n < NUM_UNITS; n++) {
                final int repLength = 1+random.nextInt(MAX_REPEAT_UNIT_LENGTH);
                final String repeatUnit = getRandomBases(repLength);
                final int numRepetitions = 1+random.nextInt(MAX_NUM_REPETITIONS);

                // log for comparison with covariate
                numsRepetitions.add(numRepetitions);
                repeatUnits.add(repeatUnit);

                for (int k=0; k < numRepetitions; k++)
                    sb.append(repeatUnit);

            }

            final String readBases = sb.toString();
            System.out.println(readBases);
            final int readLength = readBases.length();

            final byte[] readQuals = new byte[readLength];
            Arrays.fill(readQuals,(byte)30);
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(readBases.getBytes(),readQuals,readLength+"M");

            Covariate[] requestedCovariates = new Covariate[3];
            requestedCovariates[0] = rlCovariate;
            requestedCovariates[1] = ruCovariate;
            requestedCovariates[2] = rurlCovariate;
            ReadCovariates rc = RecalUtils.computeCovariates(read, requestedCovariates);

            // check that the length is correct
            Assert.assertEquals(rc.getMismatchesKeySet().length, readLength);
            Assert.assertEquals(rc.getInsertionsKeySet().length, readLength);
            Assert.assertEquals(rc.getDeletionsKeySet().length, readLength);

            for (int offset = 0; offset < readBases.length(); offset++) { // recalibrate all bases in the read
                // check RepeatLength
                final String rlValM = rlCovariate.formatKey(rc.getMismatchesKeySet(offset)[0]);
                final String rlValI = rlCovariate.formatKey(rc.getInsertionsKeySet(offset)[0]);
                final String rlValD = rlCovariate.formatKey(rc.getDeletionsKeySet(offset)[0]);
                // check RepeatUnit
                final String ruValM = ruCovariate.formatKey(rc.getMismatchesKeySet(offset)[1]);
                final String ruValI = ruCovariate.formatKey(rc.getInsertionsKeySet(offset)[1]);
                final String ruValD = ruCovariate.formatKey(rc.getDeletionsKeySet(offset)[1]);
                // check RepeatUnitAndLength
                final String rurlValM = rurlCovariate.formatKey(rc.getMismatchesKeySet(offset)[2]);
                final String rurlValI = rurlCovariate.formatKey(rc.getInsertionsKeySet(offset)[2]);
                final String rurlValD = rurlCovariate.formatKey(rc.getDeletionsKeySet(offset)[2]);
                // check all 3 values are identical
                Assert.assertEquals(rlValD,rlValI);
                Assert.assertEquals(rlValM,rlValI);
                Assert.assertEquals(ruValD,ruValI);
                Assert.assertEquals(ruValM,ruValI);
                Assert.assertEquals(rurlValD,rurlValI);
                Assert.assertEquals(rurlValM,rurlValI);


                int fw = GATKVariantContextUtils.findNumberOfRepetitions(ruValM.getBytes(), readBases.substring(offset + 1, readLength).getBytes(), true);
                int bw = GATKVariantContextUtils.findNumberOfRepetitions(ruValM.getBytes(), readBases.substring(0, offset + 1).getBytes(), false);
                Assert.assertEquals(Math.min(fw+bw,RAC.MAX_REPEAT_LENGTH),(int)Integer.valueOf(rlValM));
            }

        }






    }

    /**
     * Returns random bases of given length
     * @param length  required length
     * @return        given random string
     */
    @Requires("length > 0")
    String getRandomBases(final int length) {
        byte[] bases = new byte[length];
        Random ran = new Random();
        for (int i=0; i < length; i++ ) {
            int idx = ran.nextInt(4);
            bases[i] = BaseUtils.baseIndexToSimpleBase(idx);
        }
        return new String(bases);
    }


}
