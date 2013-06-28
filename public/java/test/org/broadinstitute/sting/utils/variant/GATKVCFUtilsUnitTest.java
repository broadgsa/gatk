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

package org.broadinstitute.sting.utils.variant;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Set;

public class GATKVCFUtilsUnitTest extends BaseTest {
    public static class VCFHeaderTestWalker extends RodWalker<Integer, Integer> {
        public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return null; }
        public Integer reduceInit() { return 0; }
        public Integer reduce(Integer value, Integer sum) { return value + sum; }
    }

    public static class VCFHeaderTest2Walker extends VCFHeaderTestWalker {}

    @Test
    public void testAddingVCFHeaderInfo() {
        final VCFHeader header = new VCFHeader();

        final Walker walker1 = new VCFHeaderTestWalker();
        final Walker walker2 = new VCFHeaderTest2Walker();

        final GenomeAnalysisEngine testEngine1 = new GenomeAnalysisEngine();
        testEngine1.setWalker(walker1);

        final GenomeAnalysisEngine testEngine2 = new GenomeAnalysisEngine();
        testEngine2.setWalker(walker2);

        final VCFHeaderLine line1 = GATKVCFUtils.getCommandLineArgumentHeaderLine(testEngine1, Collections.EMPTY_LIST);
        logger.warn(line1);
        Assert.assertNotNull(line1);
        Assert.assertEquals(line1.getKey(), GATKVCFUtils.GATK_COMMAND_LINE_KEY);
        for ( final String field : Arrays.asList("Version", "ID", "Date", "CommandLineOptions"))
            Assert.assertTrue(line1.toString().contains(field), "Couldn't find field " + field + " in " + line1.getValue());
        Assert.assertTrue(line1.toString().contains("ID=" + testEngine1.getWalkerName()));

        final VCFHeaderLine line2 = GATKVCFUtils.getCommandLineArgumentHeaderLine(testEngine2, Collections.EMPTY_LIST);
        logger.warn(line2);

        header.addMetaDataLine(line1);
        final Set<VCFHeaderLine> lines1 = header.getMetaDataInInputOrder();
        Assert.assertTrue(lines1.contains(line1));

        header.addMetaDataLine(line2);
        final Set<VCFHeaderLine> lines2 = header.getMetaDataInInputOrder();
        Assert.assertTrue(lines2.contains(line1));
        Assert.assertTrue(lines2.contains(line2));
    }
}