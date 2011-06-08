/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class BatchMergeIntegrationTest extends WalkerTest {
    @Test
    public void testBatchMerge1() {
        String bam = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
        String alleles = validationDataLocation +  "batch.merge.alleles.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T UnifiedGenotyper -NO_HEADER -BTI alleles -stand_call_conf 0.0 -glm BOTH -G none -nsl -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -o %s -R " + b37KGReference
                        + " -B:alleles,VCF " + alleles
                        + " -I " + bam,
                1,
                Arrays.asList("b7839064dc4979400af4792460d9884b"));
        executeTest("testBatchMerge UG genotype given alleles:" + new File(bam).getName() + " with " + new File(alleles).getName(), spec);
    }
}