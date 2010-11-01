package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.vcf.*;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;

import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Jun 30, 2010
 * Time: 3:32:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class VCFHeaderUnitTest extends BaseTest {

    private VCFHeader createHeader(String[] headerStr) {
        VCFCodec codec = new VCFCodec();
        List<String> headerFields = new ArrayList<String>();
        for (String str : headerStr)
            headerFields.add(str);
        VCFHeader header = (VCFHeader)codec.createHeader(headerFields,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        Assert.assertEquals(header.getMetaData().size(), headerStr.length /* for the # line */);
        return header;
    }

    @Test
    public void testVCF4ToVCF4() {
        VCFHeader header = createHeader(VCF4headerStrings);
        checkMD5ofHeaderFile(header, "4648aa1169257e0a8a9d30131adb5f35");
    }

    @Test
    public void testVCF4ToVCF4_alternate() {
        VCFHeader header = createHeader(VCF4headerStrings_with_negitiveOne);
        checkMD5ofHeaderFile(header, "ad8c4cf85e868b0261ab49ee2c613088");
    }

    private void checkMD5ofHeaderFile(VCFHeader header, String md5sum) {
        File myTempFile = null;
        PrintWriter pw = null;
        try {
            myTempFile = File.createTempFile("VCFHeader","vcf");
            myTempFile.deleteOnExit();
            pw = new PrintWriter(myTempFile);
        } catch (IOException e) {
            Assert.fail("Unable to make a temp file!");
        }
        for (VCFHeaderLine line : header.getMetaData())
            pw.println(line);
        pw.close();
        Assert.assertTrue(md5sum.equals(md5SumFile(myTempFile)));
    }

    public static String[] VCF4headerStrings = {
                "##fileformat=VCFv4.0",
                "##filedate=2010-06-21",
                "##reference=NCBI36",
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">",
                "##INFO=<ID=DP, Number=1, Type=Integer, Description=\"Total number of reads in haplotype window\">",
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">",
                "##INFO=<ID=CA, Number=1, Type=String, Description=\"Pilot 1 callability mask\">",
                "##INFO=<ID=HP, Number=1, Type=Integer, Description=\"Reference homopolymer tract length\">",
                "##INFO=<ID=NS, Number=1, Type=Integer, Description=\"Number of samples with data\">",
                "##INFO=<ID=DB, Number=0, Type=Flag, Description=\"dbSNP membership build 129 - type match and indel sequence length match within 25 bp\">",
                "##INFO=<ID=NR, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on reverse strand\">",
                "##INFO=<ID=NF, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on forward strand\">",
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">",
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">",
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">",
                "##FORMAT=<ID=GQ, Number=1, Type=Integer, Description=\"Genotype quality\">",
                };

    public static String[] VCF4headerStrings_with_negitiveOne = {
                "##fileformat=VCFv4.0",
                "##filedate=2010-06-21",
                "##reference=NCBI36",
                "##INFO=<ID=GC, Number=0, Type=Flag, Description=\"Overlap with Gencode CCDS coding sequence\">",
                "##INFO=<ID=YY, Number=., Type=Integer, Description=\"Some weird value that has lots of parameters\">",
                "##INFO=<ID=AF, Number=1, Type=Float, Description=\"Dindel estimated population allele frequency\">",
                "##INFO=<ID=CA, Number=1, Type=String, Description=\"Pilot 1 callability mask\">",
                "##INFO=<ID=HP, Number=1, Type=Integer, Description=\"Reference homopolymer tract length\">",
                "##INFO=<ID=NS, Number=1, Type=Integer, Description=\"Number of samples with data\">",
                "##INFO=<ID=DB, Number=0, Type=Flag, Description=\"dbSNP membership build 129 - type match and indel sequence length match within 25 bp\">",
                "##INFO=<ID=NR, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on reverse strand\">",
                "##INFO=<ID=NF, Number=1, Type=Integer, Description=\"Number of reads covering non-ref variant on forward strand\">",
                "##FILTER=<ID=NoQCALL, Description=\"Variant called by Dindel but not confirmed by QCALL\">",
                "##FORMAT=<ID=GT, Number=1, Type=String, Description=\"Genotype\">",
                "##FORMAT=<ID=HQ, Number=2, Type=Integer, Description=\"Haplotype quality\">",
                "##FORMAT=<ID=TT, Number=., Type=Integer, Description=\"Lots of TTs\">",
                };
}
