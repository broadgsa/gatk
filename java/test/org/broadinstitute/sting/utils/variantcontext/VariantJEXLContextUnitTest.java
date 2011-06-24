/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext;

import net.sf.samtools.SAMFileHeader;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class VariantJEXLContextUnitTest
 *
 * Test out parts of the VariantJEXLContext
 */
public class VariantJEXLContextUnitTest extends BaseTest {


    private static String expression = "QUAL > 500.0";
    private static VariantContextUtils.JexlVCMatchExp exp;

    Allele A, Aref, T, Tref;

    Allele del, delRef, ATC, ATCref;
    // A [ref] / T at 10

    GenomeLoc snpLoc;
    // - / ATC [ref] from 20-23

    private static int startingChr = 1;
    private static int endingChr = 2;
    private static int readCount = 100;
    private static int DEFAULT_READ_LENGTH = ArtificialSAMUtils.DEFAULT_READ_LENGTH;
    static SAMFileHeader header;

    private static GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        try {
            exp = new VariantContextUtils.JexlVCMatchExp("name", VariantContextUtils.engine.createExpression(expression));
        } catch (Exception e) {
            Assert.fail("Unable to create expression" + e.getMessage());
        }
        snpLoc = genomeLocParser.createGenomeLoc("chr1", 10, 10, true);
    }

    @BeforeMethod
    public void before() {
        del = Allele.create("-");
        delRef = Allele.create("-", true);

        A = Allele.create("A");
        Aref = Allele.create("A", true);
        T = Allele.create("T");
        Tref = Allele.create("T", true);

        ATC = Allele.create("ATC");
        ATCref = Allele.create("ATC", true);
    }


    @Test
    public void testGetValue() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        // make sure the context has a value
        Assert.assertTrue(!map.isEmpty());
        Assert.assertEquals(map.size(), 1);

        // eval our known expression
        Assert.assertTrue(!map.get(exp));
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testContainsValue() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.containsValue(exp);
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testRemove() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.remove(exp);
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testEntrySet() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.entrySet();
    }

    @Test(expectedExceptions=UnsupportedOperationException.class)
    public void testClear() {
        Map<VariantContextUtils.JexlVCMatchExp, Boolean> map = getVarContext();

        map.clear();
    }

    /**
     * helper method
     * @return a VariantJEXLContext
     */
    private JEXLMap getVarContext() {
        List<Allele> alleles = Arrays.asList(Aref, T);

        VariantContext vc = new VariantContext("test", snpLoc.getContig(), snpLoc.getStart(), snpLoc.getStop(), alleles);
        return new JEXLMap(Arrays.asList(exp),vc);
    }


}
