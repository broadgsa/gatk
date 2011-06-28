package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.testng.Assert;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeaderUnitTest;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;

import org.testng.annotations.Test;

import java.io.StringBufferInputStream;
import java.util.ArrayList;
import java.util.Set;

/**
 * test out pieces of the combine variants code
 */
public class CombineVariantsUnitTest {

    // this header is a small subset of the header in VCFHeaderUnitTest: VCF4headerStrings
    public static String VCF4headerStringsSmallSubset =
                "##fileformat=VCFv4.0\n" +
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    // altered info field
    public static String VCF4headerStringsBrokenInfo =
                "##fileformat=VCFv4.0\n"+
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=String, Description=\"Dindel estimated population allele frequency\">\n"+ // string to integer
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">\n"+
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    // altered format field
    public static String VCF4headerStringsBrokenFormat =
                "##fileformat=VCFv4.0\n"+
                "##filedate=2010-06-21\n"+
                "##reference=NCBI36\n"+
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">\n"+
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">\n"+
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">\n"+
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">\n"+
                "##FORMAT=<ID=GT, Number=6, Type=String, Description=\"Genotype\">\n"+ // changed 1 to 6 here
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">\n"+
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">\n"+
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    private VCFHeader createHeader(String headerStr) {
        VCFCodec codec = new VCFCodec();
        VCFHeader head = (VCFHeader)codec.readHeader(new AsciiLineReader(new StringBufferInputStream(headerStr)));
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
        Assert.assertEquals(lines.size(), VCFHeaderUnitTest.VCF4headerStringCount);
    }

    @Test(expectedExceptions=IllegalStateException.class)
    public void testHeadersInfoDifferentValues() {
        VCFHeader one = createHeader(VCFHeaderUnitTest.VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsBrokenInfo);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, null);
        Assert.assertEquals(lines.size(), VCFHeaderUnitTest.VCF4headerStringCount);
    }

    @Test
    public void testHeadersFormatDifferentValues() {
        VCFHeader one = createHeader(VCFHeaderUnitTest.VCF4headerStrings);
        VCFHeader two = createHeader(VCF4headerStringsBrokenFormat);
        ArrayList<VCFHeader> headers = new ArrayList<VCFHeader>();
        headers.add(one);
        headers.add(two);
        Set<VCFHeaderLine> lines = VCFUtils.smartMergeHeaders(headers, null);
        Assert.assertEquals(lines.size(), VCFHeaderUnitTest.VCF4headerStringCount);
    }
}
