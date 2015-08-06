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

package org.broadinstitute.gatk.tools.walkers.filters;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.Utils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.testng.Assert;
import org.testng.annotations.*;

public class VariantFiltrationUnitTest extends BaseTest {

    private String chr1 = null;
    private GenomeLoc genomeLoc = null;
    private String vcFilter = "testFilter";

    @BeforeTest
    public void before() {
        // Create GenomeLoc
        IndexedFastaSequenceFile fasta = CachingIndexedFastaSequenceFile.checkAndCreate(new File(privateTestDir + "iupacFASTA.fasta"));
        GenomeLocParser genomeLocParser = new GenomeLocParser(fasta);
        chr1 = fasta.getSequenceDictionary().getSequence(0).getSequenceName();
        genomeLoc = genomeLocParser.createGenomeLoc(chr1, 5, 10);
    }

    @DataProvider(name = "VariantMaskData")
    public Object[][] DoesMaskCoverVariantTestData() {

        final String maskName = "testMask";

        List<Object[]> tests = Arrays.asList(new Object[]{chr1, 0, 0, maskName, 10, true, true},
                                             new Object[]{"chr2", 0, 0, maskName, 10, true, false},
                                             new Object[]{chr1, 0, 0, null, 10, true, true},
                                             new Object[]{chr1, 0, 0, maskName, 10, true, true},
                                             new Object[]{chr1, 0, 0, vcFilter, 10, true, false},
                                             new Object[]{chr1, 0, 0, maskName, 1, true, false},
                                             new Object[]{chr1, 15, 15, maskName, 10, false, true},
                                             new Object[]{chr1, 15, 15, maskName, 1, false, false}
                                            );
        return tests.toArray(new Object[][]{});
    }

    /**
     * Test doesMaskCoverVariant() logic
     *
     * @param contig chromosome or contig name
     * @param start  variant context start
     * @param stop variant context stop
     * @param maskName mask or filter name
     * @param maskExtension bases beyond the mask
     * @param vcBeforeLoc if true, variant context is before the genome location; if false, the converse is true.
     * @param expectedValue  return the expected return value from doesMaskCoverVariant()
     */
    @Test(dataProvider = "VariantMaskData")
    public void TestDoesMaskCoverVariant(final String contig, final int start, final int stop, final String maskName, final int maskExtension,
                                         final boolean vcBeforeLoc, final boolean expectedValue) {

        // Build VariantContext
        final byte[] allele1 = Utils.dupBytes((byte) 'A', 1);
        final byte[] allele2 = Utils.dupBytes((byte) 'T', 2);

        final List<Allele> alleles = new ArrayList<Allele>(2);
        final Allele ref = Allele.create(allele1, true);
        final Allele alt = Allele.create(allele2, false);
        alleles.add(ref);
        alleles.add(alt);

        final VariantContext vc = new VariantContextBuilder("test", contig, start, stop, alleles).filter(vcFilter).make();

        boolean coversVariant = VariantFiltration.doesMaskCoverVariant(vc, genomeLoc, maskName, maskExtension, vcBeforeLoc);
        Assert.assertEquals(coversVariant, expectedValue);
    }
}
