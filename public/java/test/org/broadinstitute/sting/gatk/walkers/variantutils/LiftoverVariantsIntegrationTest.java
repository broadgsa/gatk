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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Tests LiftoverVariants
 */
public class LiftoverVariantsIntegrationTest extends WalkerTest {

    @Test
    public void testb36Tohg19() {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T LiftoverVariants -o %s -R " + b36KGReference + " --variant " + privateTestDir + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.noheader.vcf -chain " + validationDataLocation + "b36ToHg19.broad.over.chain -dict /seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19.dict",
                 1,
                 Arrays.asList("b8f4171b0c39954b283dfed4afed87d7"));
         executeTest("test b36 to hg19", spec);
    }

    @Test
    public void testb36Tohg19UnsortedSamples() {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T LiftoverVariants -o %s -R " + b36KGReference + " --variant " + privateTestDir + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.noheader.unsortedSamples.vcf -chain " + validationDataLocation + "b36ToHg19.broad.over.chain -dict /seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19.dict",
                 1,
                 Arrays.asList("2fd072a88006e3d7cfea4cac5d7b0d98"));
         executeTest("test b36 to hg19, unsorted samples", spec);
    }

    @Test
    public void testhg18Tohg19Unsorted() {
         WalkerTestSpec spec = new WalkerTestSpec(
                 "-T LiftoverVariants -o %s -R " + hg18Reference + " --variant:vcf " + privateTestDir + "liftover_test.vcf -chain " + validationDataLocation + "hg18ToHg19.broad.over.chain -dict /seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19.dict",
                 1,
                 Arrays.asList("e0b813ff873185ab51995a151f80ec98"));
         executeTest("test hg18 to hg19, unsorted", spec);
    }
}
