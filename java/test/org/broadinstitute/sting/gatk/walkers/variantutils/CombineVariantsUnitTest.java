package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.vcf.VCFCodec;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.testng.Assert;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeaderUnitTest;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * test out pieces of the combine variants code
 */
public class CombineVariantsUnitTest {

    // this header is a small subset of the header in VCFHeaderUnitTest: VCF4headerStrings
    public static String[] VCF4headerStringsSmallSubset = {
                "##fileformat=VCFv4.0",
                "##filedate=2010-06-21",
                "##reference=NCBI36",
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">",
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">",
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">",
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">",
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">",
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">",
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">",
                };

    // altered info field
    public static String[] VCF4headerStringsBrokenInfo = {
                "##fileformat=VCFv4.0",
                "##filedate=2010-06-21",
                "##reference=NCBI36",
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">",
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">",
                "##INFO=<ID=AF, Number=1, Type=Integer, Description=\"Dindel estimated population allele frequency\">", // float to integer
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">",
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">",
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">",
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">",
                };

    // altered format field
    public static String[] VCF4headerStringsBrokenFormat = {
                "##fileformat=VCFv4.0",
                "##filedate=2010-06-21",
                "##reference=NCBI36",
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">",
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">",
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">",
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">",
                "##FORMAT=<ID=GT, Number=6, Type=String, Description=\"Genotype\">", // changed 1 to 6 here
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">",
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">",
                };

    private VCFHeader createHeader(String[] headerStr) {
        VCFCodec codec = new VCFCodec();
        List<String> headerFields = new ArrayList<String>();
        for (String str : headerStr)
            headerFields.add(str);
        VCFHeader head = (VCFHeader)codec.createHeader(headerFields,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        Assert.assertEquals(head.getMetaData().size(), headerStr.length /* for the # line */);
        return head;
    }

    @Test
    public void testHeadersWhereOneIsAStrictSubsetOfTheOther() {
        VCFHeader one = createHeader(VCFHeaderUnitTest.VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsSmallSubset);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, null);
        Assert.assertEquals(lines.size(), VCFHeaderUnitTest.VCF4headerStrings.length);
    }

    @Test(expectedExceptions=IllegalStateException.class)
    public void testHeadersInfoDifferentValues() {
        VCFHeader one = createHeader(VCFHeaderUnitTest.VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsBrokenInfo);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, null);
        Assert.assertEquals(lines.size(), VCFHeaderUnitTest.VCF4headerStrings.length);
    }

    @Test
    public void testHeadersFormatDifferentValues() {
        VCFHeader one = createHeader(VCFHeaderUnitTest.VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsBrokenFormat);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, null);
        Assert.assertEquals(lines.size(), VCFHeaderUnitTest.VCF4headerStrings.length);
    }
}
