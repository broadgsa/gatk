/*
 * Copyright (c) 2011, The Broad Institute
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

// our package
package org.broadinstitute.sting.gatk.walkers.diffengine;


// the imports for unit testing.


import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DiffableReaderUnitTest extends BaseTest {
    DiffEngine engine;

    File vcfFile = new File(testDir + "diffTestMaster.vcf");
    File bamFile = new File(testDir + "exampleBAM.bam");

    @BeforeClass(enabled = true)
    public void createDiffEngine() {
        engine = new DiffEngine();
    }

    @Test(enabled = true)
    public void testPluggableDiffableReaders() {
        logger.warn("testPluggableDiffableReaders");
        Map<String, DiffableReader> readers = engine.getReaders();
        Assert.assertNotNull(readers);
        Assert.assertTrue(readers.size() > 0);
        Assert.assertNotNull(readers.get("VCF"));
        for ( Map.Entry<String, DiffableReader> e : engine.getReaders().entrySet() ) {
            logger.warn("Found diffable reader: " + e.getKey());
            Assert.assertEquals(e.getValue().getName(), e.getKey());
            Assert.assertEquals(e.getValue(), engine.getReader(e.getKey()));
        }
    }

    private static void testLeaf(DiffNode rec, String field, Object expected) {
        DiffElement value = rec.getElement(field);
        Assert.assertNotNull(value, "Expected to see leaf named " + field + " in rec " + rec);
        Assert.assertEquals(value.getValue().getValue(), expected, "Expected to see leaf named " + field + " to have value " + expected + " in rec " + rec + " but got instead " + value.getValue().getValue());
    }

    @Test(enabled = true, dependsOnMethods = "testPluggableDiffableReaders")
    public void testVCF1() {
        logger.warn("testVCF1");
        DiffableReader vcfReader = engine.getReader("VCF");
        Assert.assertTrue(vcfReader.canRead(vcfFile));
        Assert.assertFalse(vcfReader.canRead(bamFile));

        DiffElement diff = vcfReader.readFromFile(vcfFile, -1);
        Assert.assertNotNull(diff);

        Assert.assertEquals(diff.getName(), vcfFile.getName());
        Assert.assertSame(diff.getParent(), DiffElement.ROOT);

        DiffNode node = diff.getValueAsNode();
        Assert.assertEquals(node.getElements().size(), 11);

        // chr1    2646    rs62635284      G       A       0.15    PASS    AC=2;AF=1.00;AN=2       GT:AD:DP:GL:GQ  1/1:53,75:3:-12.40,-0.90,-0.00:9.03
        DiffNode rec1 = node.getElement("chr1:2646").getValueAsNode();
        testLeaf(rec1, "CHROM", "chr1");
        testLeaf(rec1, "POS", 2646);
        testLeaf(rec1, "ID", "rs62635284");
        testLeaf(rec1, "REF", Allele.create("G", true));
        testLeaf(rec1, "ALT", Arrays.asList(Allele.create("A")));
        testLeaf(rec1, "QUAL", 0.15);
        testLeaf(rec1, "FILTER", Collections.<Object>emptySet());
        testLeaf(rec1, "AC", "2");
        testLeaf(rec1, "AF", "1.00");
        testLeaf(rec1, "AN", "2");
    }

    @Test(enabled = true, dependsOnMethods = "testPluggableDiffableReaders")
    public void testBAM() {
        logger.warn("testBAM");
        DiffableReader bamReader = engine.getReader("BAM");
        Assert.assertTrue(bamReader.canRead(bamFile));
        Assert.assertFalse(bamReader.canRead(vcfFile));

        DiffElement diff = bamReader.readFromFile(bamFile, -1);
        Assert.assertNotNull(diff);

        Assert.assertEquals(diff.getName(), bamFile.getName());
        Assert.assertSame(diff.getParent(), DiffElement.ROOT);

        DiffNode node = diff.getValueAsNode();
        Assert.assertEquals(node.getElements().size(), 33);

        // 30PPJAAXX090125:1:42:512:1817#0 99      chr1    200     0       76M     =
        // 255     -130    ACCCTAACCCTAACCCTAACCCTAACCATAACCCTAAGACTAACCCTAAACCTAACCCTCATAATCGAAATACAAC
        // BBBBC@C?AABCBB<63>=B@>+B9-9+)2B8,+@327B5A>90((>-+''3?(/'''A)(''19('7.,**%)3:
        // PG:Z:0  RG:Z:exampleBAM.bam     SM:Z:exampleBAM.bam

        DiffNode rec1 = node.getElement("30PPJAAXX090125:1:42:512:1817#0_1").getValueAsNode();
        testLeaf(rec1, "NAME", "30PPJAAXX090125:1:42:512:1817#0");
        testLeaf(rec1, "FLAGS", 99);
        testLeaf(rec1, "RNAME", "chr1");
        testLeaf(rec1, "POS", 200);
        testLeaf(rec1, "MAPQ", 0);
        testLeaf(rec1, "CIGAR", "76M");
        testLeaf(rec1, "RNEXT", "chr1");
        testLeaf(rec1, "PNEXT", 255);
        testLeaf(rec1, "TLEN", -130);
        testLeaf(rec1, "SEQ", "ACCCTAACCCTAACCCTAACCCTAACCATAACCCTAAGACTAACCCTAAACCTAACCCTCATAATCGAAATACAAC");
        testLeaf(rec1, "QUAL", "BBBBC@C?AABCBB<63>=B@>+B9-9+)2B8,+@327B5A>90((>-+''3?(/'''A)(''19('7.,**%)3:");
        testLeaf(rec1, "PG", "0");
        testLeaf(rec1, "RG", "exampleBAM.bam");
        testLeaf(rec1, "SM", "exampleBAM.bam");
    }
}