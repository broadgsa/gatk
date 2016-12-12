/*
* Copyright 2012-2016 Broad Institute, Inc.
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

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class VariantFiltrationUnitTest extends BaseTest {

    private String chr1 = null;
    private GenomeLoc genomeLoc = null;
    private String vcFilter = "testFilter";

    @BeforeTest
    public void before() {
        // Create GenomeLoc
        ReferenceSequenceFile fasta = CachingIndexedFastaSequenceFile.checkAndCreate(new File(privateTestDir + "iupacFASTA.fasta"));
        GenomeLocParser genomeLocParser = new GenomeLocParser(fasta);
        chr1 = fasta.getSequenceDictionary().getSequence(0).getSequenceName();
        genomeLoc = genomeLocParser.createGenomeLoc(chr1, 5, 10);
    }

    @DataProvider(name = "VariantMaskData")
    public Object[][] doesMaskCoverVariantTestData() {

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
    @Test(dataProvider = "doesMaskCoverVariantTestData")
    public void testDoesMaskCoverVariant(final String contig, final int start, final int stop, final String maskName, final int maskExtension,
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

    @Test
    public void testApplyGenotypeFilters(){

        final VariantContext vc = buildDataForFilters().make();

        final String filterName = "LowZ"; //an attribute that doesn't appear in the VariantContext, so there isn't any chance of confusion like with the INFO DP
        final String filterExpr = "Z < 10";

        final List<VariantContextUtils.JexlVCMatchExp> genotypeFilterExps = VariantContextUtils.initializeMatchExps(Arrays.asList(filterName), Arrays.asList(filterExpr));

        final VariantContextBuilder anotherVCBuilder = VariantFiltration.applyGenotypeFilters(vc, genotypeFilterExps, false, false, false);
        final VariantContext anotherVC = anotherVCBuilder.filters().make();

        Assert.assertEquals(anotherVC.getGenotype("one").isFiltered(), true);
        Assert.assertTrue(anotherVC.getGenotype("one").getFilters().equals(filterName));

        Assert.assertEquals(anotherVC.getGenotype("two").isFiltered(), false);

        Assert.assertEquals(anotherVC.getGenotype("three").isFiltered(), false);

        Assert.assertEquals(anotherVC.getGenotype("four").isFiltered(), false);

        Assert.assertEquals(anotherVC.getGenotype("five").isFiltered(), false);

        Assert.assertEquals(anotherVC.getGenotype("six").isFiltered(), false);

        final VariantContextBuilder yetAnotherVCBuilder = VariantFiltration.applyGenotypeFilters(anotherVC, genotypeFilterExps, false, true, false);
        final VariantContext yetAnotherVC = yetAnotherVCBuilder.filters().make();
        Assert.assertEquals(yetAnotherVC.getGenotype("six").isFiltered(), true);
        Assert.assertTrue(yetAnotherVC.getGenotype("six").getFilters().equals(filterName));
    }

    @Test
    public void testApplyVCFilters(){

        final VariantContext vcNoFilters = buildDataForFilters().make(); // assumes this vc doesn't hold any filters yet

        String filterName = "LowDP";
        String filterExpr = "DP < 23";
        List<VariantContextUtils.JexlVCMatchExp> vcFilterExps = VariantContextUtils.initializeMatchExps(Arrays.asList(filterName), Arrays.asList(filterExpr));

        final Set<String> filters = VariantFiltration.buildVCfilters(vcNoFilters, vcFilterExps, false, false);
        Assert.assertFalse(vcNoFilters.isFiltered());
        Assert.assertEquals(filters.size(), 1);
        Assert.assertTrue(filters.contains(filterName));

        filterName = "ID";
        filterExpr = "ID = rs123";
        vcFilterExps = VariantContextUtils.initializeMatchExps(Arrays.asList(filterName), Arrays.asList(filterExpr));
        Set<String> filterWhenFailMissing = VariantFiltration.buildVCfilters(vcNoFilters, vcFilterExps, false, true);
//        Assert.assertEquals(filterWhenFailMissing.size(), 1);
//        Assert.assertTrue(filterWhenFailMissing.contains(filterName));
        filterWhenFailMissing = VariantFiltration.buildVCfilters(vcNoFilters, vcFilterExps, false, false);
        Assert.assertTrue(filterWhenFailMissing.isEmpty());
    }


    private static VariantContextBuilder buildDataForFilters() {
        /**
         * Uses (part of) the following (semi fake) data for testing (data was modified from real data so expect some minor inconsistencies in annotations)
         * 1    1234567    .    T    C    152.03  .
         *          AC=6;AF=1.00;AN=6;DP=22;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;SOR=0.693;set=variant3-variant4-variant6
         *          GT:AD:DP:GQ:PGT:PID:PL:RGQ
         *          1/1:0,2:9:6:1|1:15870493_CT_C:90,6,0
         *          1/1:0,4:10:12:1|1:15870493_CT_C:180,12,0
         *          ./.:0:3:.:.:.:.:0
         *          ./.:0:0:.:.:.:.:0
         *          ./.:0:0:.:.:.:.:0
         *          1/1:0,0:.:6:1|1:15870493_CT_C:90,6,0
         */

        final Allele refT = Allele.create("T", true);
        final Allele altC = Allele.create("C", false);
        final Allele nocall = Allele.NO_CALL;

        final VariantContextBuilder vcBuilder = new VariantContextBuilder("", "1", 1234567, 1234567, Arrays.asList(refT, altC));

        vcBuilder.noID();
        vcBuilder.attribute("AC", 6);
        vcBuilder.attribute("AF", 1.00);
        vcBuilder.attribute("AN", 6);
        vcBuilder.attribute("DP", 22);
        vcBuilder.attribute("ExcessHet", 3.0103);
        vcBuilder.attribute("FS", 0.000);
        vcBuilder.attribute("MLEAC", 2);
        vcBuilder.attribute("MLEAF", 1.00);
        vcBuilder.attribute("SOR", 0.693);

        GenotypeBuilder gtBuilder = new GenotypeBuilder("one", Arrays.asList(altC,altC));
        final String ATTRIBUTE_NOT_IN_VC = "Z";
        final Genotype firstSample = gtBuilder.attribute(VCFConstants.GENOTYPE_KEY, GenotypeType.HOM_VAR)
                                                .attribute(ATTRIBUTE_NOT_IN_VC, 9) // edge case not passing "Z < 10"
                                                .make();

        gtBuilder = new GenotypeBuilder("two", Arrays.asList(altC,altC));
        final Genotype secondSample = gtBuilder.attribute(VCFConstants.GENOTYPE_KEY, GenotypeType.HOM_VAR)
                                                .attribute(ATTRIBUTE_NOT_IN_VC, 10) // edge case passing    "Z = 10"
                                                .make();

        gtBuilder = new GenotypeBuilder("three", Arrays.asList(nocall,nocall));
        final Genotype thirdSample = gtBuilder.attribute(VCFConstants.GENOTYPE_KEY, GenotypeType.NO_CALL)
                                                .attribute(ATTRIBUTE_NOT_IN_VC, 3)
                                                .make();

        gtBuilder = new GenotypeBuilder("four", Arrays.asList(nocall,nocall));
        final Genotype fourthSample = gtBuilder.attribute(VCFConstants.GENOTYPE_KEY, GenotypeType.NO_CALL)
                                                .attribute(ATTRIBUTE_NOT_IN_VC, (0))
                                                .make();

        gtBuilder = new GenotypeBuilder("five", Arrays.asList(nocall,nocall));
        final Genotype fifthSample = gtBuilder.attribute(VCFConstants.GENOTYPE_KEY, GenotypeType.NO_CALL)
                                                .attribute(ATTRIBUTE_NOT_IN_VC, 0)
                                                .make();

        gtBuilder = new GenotypeBuilder("six", Arrays.asList(altC,altC));
        final Genotype sixthSample = gtBuilder.attribute(VCFConstants.GENOTYPE_KEY, GenotypeType.HOM_VAR)
                                                //no Z
                                                .make();

        vcBuilder.genotypes(firstSample, secondSample, thirdSample, fourthSample, fifthSample, sixthSample);
        return vcBuilder;
    }
}
