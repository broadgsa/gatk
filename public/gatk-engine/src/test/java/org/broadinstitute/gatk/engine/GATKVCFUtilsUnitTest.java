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

package org.broadinstitute.gatk.engine;

import htsjdk.tribble.index.DynamicIndexCreator;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.interval.IntervalIndexCreator;
import htsjdk.tribble.index.linear.LinearIndexCreator;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.Walker;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.reflect.Method;
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

        final VCFHeaderLine line1 = GATKVCFUtils.getCommandLineArgumentHeaderLine(header, testEngine1, Collections.EMPTY_LIST);
        logger.warn(line1);
        Assert.assertNotNull(line1);
        // assert the key matches the expected format (GATKVCFUtils.GATK_COMMAND_LINE_KEY).(walker name)
        final String expectedLine1Key = String.format("%s.%s", GATKVCFUtils.GATK_COMMAND_LINE_KEY, testEngine1.getWalkerName());
        Assert.assertEquals(line1.getKey(), expectedLine1Key);

        for (final String field : Arrays.asList("Version", "ID", "Date", "CommandLineOptions"))
            Assert.assertTrue(line1.toString().contains(field), "Couldn't find field " + field + " in " + line1.getValue());
        Assert.assertTrue(line1.toString().contains("ID=" + testEngine1.getWalkerName()));

        final VCFHeaderLine line2 = GATKVCFUtils.getCommandLineArgumentHeaderLine(header, testEngine2, Collections.EMPTY_LIST);
        logger.warn(line2);


        header.addMetaDataLine(line1);
        final Set<VCFHeaderLine> lines1 = header.getMetaDataInInputOrder();
        Assert.assertTrue(lines1.contains(line1));

        header.addMetaDataLine(line2);
        final Set<VCFHeaderLine> lines2 = header.getMetaDataInInputOrder();
        Assert.assertTrue(lines2.contains(line1));
        Assert.assertTrue(lines2.contains(line2));

        // create a new header line using the same engine as used by line 1
        final VCFHeaderLine line3 = GATKVCFUtils.getCommandLineArgumentHeaderLine(header, testEngine1, Collections.EMPTY_LIST);
        logger.warn(line3);

        // ensure convention followed by getCommandLineArgumentHeaderLine is to append ".(number of duplicate engine runs)"
        // line3 uses the same walker as line1, whereas line2 uses a different walker.  line3 is the second occurrence of walker1
        // so a ".2" gets appended afterwards
        final String expectedLine3Key = String.format("%s.%s.2", GATKVCFUtils.GATK_COMMAND_LINE_KEY, testEngine1.getWalkerName());
        Assert.assertEquals(line3.getKey(), expectedLine3Key);

        header.addMetaDataLine(line3);

        final Set<VCFHeaderLine> lines3 = header.getMetaDataInInputOrder();
        Assert.assertTrue(lines3.contains(line1));
        Assert.assertTrue(lines3.contains(line2));
        Assert.assertTrue(lines3.contains(line3));
    }

    private class IndexCreatorTest extends TestDataProvider {
        private final GATKVCFIndexType type;
        private final int parameter;
        private final Class expectedClass;
        private final Integer expectedDimension;
        private final Method dimensionGetter;

        private IndexCreatorTest(GATKVCFIndexType type, int parameter, Class expectedClass, Integer expectedDimension,
                                 String dimensionGetterName) {
            super(IndexCreatorTest.class);

            this.type = type;
            this.parameter = parameter;
            this.expectedClass = expectedClass;
            this.expectedDimension = expectedDimension;
            try {
                // Conditional matches testGetIndexCreator's if-statement
                this.dimensionGetter = this.expectedDimension == null ? null : expectedClass.getDeclaredMethod(dimensionGetterName);
            } catch (NoSuchMethodException e) {
                throw new RuntimeException(e);
            }
        }
    }

    @DataProvider(name = "indexCreator")
    public Object[][] indexCreatorData() {
        new IndexCreatorTest(GATKVCFIndexType.DYNAMIC_SEEK, 0, DynamicIndexCreator.class, null, null);
        new IndexCreatorTest(GATKVCFIndexType.DYNAMIC_SIZE, 0, DynamicIndexCreator.class, null, null);
        new IndexCreatorTest(GATKVCFIndexType.LINEAR, 100, LinearIndexCreator.class, 100, "getBinSize");
        new IndexCreatorTest(GATKVCFIndexType.INTERVAL, 200, IntervalIndexCreator.class, 200, "getFeaturesPerInterval");

        return IndexCreatorTest.getTests(IndexCreatorTest.class);
    }

    @Test(dataProvider = "indexCreator")
    public void testGetIndexCreator(IndexCreatorTest spec) throws Exception{
        File dummy = new File("");
        IndexCreator ic = GATKVCFUtils.getIndexCreator(spec.type, spec.parameter, dummy);
        Assert.assertEquals(ic.getClass(), spec.expectedClass, "Wrong IndexCreator type");
        if (spec.expectedDimension != null) {
            Integer dimension = (int)spec.dimensionGetter.invoke(ic);
            Assert.assertEquals(dimension, spec.expectedDimension, "Wrong dimension");
        }
    }
}