/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.beagle;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;
 
public class BeagleIntegrationTest extends WalkerTest {

    private static final String beagleValidationDataLocation = validationDataLocation + "/Beagle/";
    @Test
    public void testBeagleOutput() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T BeagleOutputToVCF -R " + hg19Reference + " " +
                        "-B:variant,VCF " + beagleValidationDataLocation + "inttestbgl.input.vcf " +
                        "-B:beagleR2,BEAGLE " + beagleValidationDataLocation + "inttestbgl.r2 " +
                        "-B:beagleProbs,BEAGLE " + beagleValidationDataLocation + "inttestbgl.gprobs " +
                        "-B:beaglePhased,BEAGLE " + beagleValidationDataLocation + "inttestbgl.phased " +
                        "-o %s -NO_HEADER", 1, Arrays.asList("6bccee48ad2f06ba5a8c774fed444478"));
        executeTest("test BeagleOutputToVCF", spec);
    }
   
    @Test
    public void testBeagleInput() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T ProduceBeagleInput -R " + hg19Reference + " " +
                        "-B:variant,VCF " + beagleValidationDataLocation + "inttestbgl.input.vcf " +
                         "-o %s", 1, Arrays.asList("c0d30e5dbe903874f8422a0e63a5118e"));
        executeTest("test BeagleInput", spec);
    }

    @Test
    public void testBeagleInput2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T ProduceBeagleInput -B:variant,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878_HSQ_chr22_14-16m.vcf "+
                        "-B:validation,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878_OMNI_chr22_14-16m.vcf "+
                        "-L 22:14000000-16000000 -o %s -bvcf %s -bs 0.8 -valp 0.98 -R /humgen/1kg/reference/human_g1k_v37.fasta -NO_HEADER ",2,
                Arrays.asList("44d28b6b092d5f4c0ae59af442612ea3","223fb977e8db567dcaf632c6ee51f294"));
        executeTest("test BeagleInputWithBootstrap",spec);
    }

    @Test
    public void testBeagleOutput2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T BeagleOutputToVCF -R "+hg19Reference+" "+
                "-B:variant,VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.vcf "+
                "-B:beagleR2,beagle /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.r2 "+
                "-B:beagleProbs,beagle /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.gprobs.bgl "+
                "-B:beaglePhased,beagle /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.phased.bgl "+
                "-L 20:1-70000 -o %s -NO_HEADER ",1,Arrays.asList("e09d2cba18ab5c571e563085296cd514"));

        executeTest("testBeagleChangesSitesToRef",spec);
    }

}
