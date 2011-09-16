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
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;
import java.io.File;

public class VariantsToTableIntegrationTest extends WalkerTest {
    private String variantsToTableCmd(String moreArgs) {
        return "-R " + hg18Reference +
                " --variant:vcf " + validationDataLocation + "/soap_gatk_annotated.vcf" +
                " -T VariantsToTable" +
                " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F TRANSITION -F DP -F SB -F set -F RankSumP -F refseq.functionalClass*" +
                " -L chr1 -KMA -o %s" + moreArgs;
    }

    @Test(enabled = true)
    public void testComplexVariantsToTable() {
        WalkerTestSpec spec = new WalkerTestSpec(variantsToTableCmd(" -AMD"),
                Arrays.asList("e8f771995127b727fb433da91dd4ee98"));
        executeTest("testComplexVariantsToTable", spec).getFirst();
    }

    @Test(enabled = true)
    public void testComplexVariantsToTableFail() {
        WalkerTestSpec spec = new WalkerTestSpec(variantsToTableCmd(""), 1, UserException.class);
        executeTest("testComplexVariantsToTable-FAIL", spec);
    }
}
