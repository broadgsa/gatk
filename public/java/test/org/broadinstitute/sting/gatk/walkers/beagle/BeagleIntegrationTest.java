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
                        "--variant:VCF3 " + beagleValidationDataLocation + "inttestbgl.input.vcf " +
                        "--beagleR2:BEAGLE " + beagleValidationDataLocation + "inttestbgl.r2 " +
                        "--beagleProbs:BEAGLE " + beagleValidationDataLocation + "inttestbgl.gprobs " +
                        "--beaglePhased:BEAGLE " + beagleValidationDataLocation + "inttestbgl.phased " +
                        "-o %s -NO_HEADER", 1, Arrays.asList("b445d280fd8fee1eeb4aacb3f5a54847"));
        executeTest("test BeagleOutputToVCF", spec);
    }
   
    @Test
    public void testBeagleInput() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T ProduceBeagleInput -R " + hg19Reference + " " +
                        "--variant:VCF3 " + beagleValidationDataLocation + "inttestbgl.input.vcf " +
                         "-o %s", 1, Arrays.asList("a01c704246f3dd1b9c65774007e51e69"));
        executeTest("test BeagleInput", spec);
    }

    @Test
    public void testBeagleInput2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T ProduceBeagleInput --variant:VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878_HSQ_chr22_14-16m.vcf "+
                        "--validation:VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878_OMNI_chr22_14-16m.vcf "+
                        "-L 22:14000000-16000000 -o %s -bvcf %s -bs 0.8 -valp 0.98 -R /humgen/1kg/reference/human_g1k_v37.fasta -NO_HEADER ",2,
                Arrays.asList("660986891b30cdc937e0f2a3a5743faa","e96ddd51da9f4a797b2aa8c20e404166"));
        executeTest("test BeagleInputWithBootstrap",spec);
    }

    @Test
    public void testBeagleOutput2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T BeagleOutputToVCF -R "+hg19Reference+" "+
                "--variant:VCF /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.vcf "+
                "--beagleR2:beagle /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.r2 "+
                "--beagleProbs:beagle /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.gprobs.bgl "+
                "--beaglePhased:beagle /humgen/gsa-hpprojects/GATK/data/Validation_Data/EUR_beagle_in_test.phased.bgl "+
                "-L 20:1-70000 -o %s -NO_HEADER ",1,Arrays.asList("51a57ea565176edd96d907906914b0ee"));

        executeTest("testBeagleChangesSitesToRef",spec);
    }

}
