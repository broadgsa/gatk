package org.broadinstitute.sting.gatk.refdata.features.vcf4;

import org.broad.tribble.util.AsciiLineReader;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * test out pieces of the VCF 4 codec.
 */
public class VCF4UnitTest extends BaseTest {
    File vcfGenotypeFile = new File("testdata/vcf/vcfWithGenotypes.vcf");
    File vcfNoGenotypeFile = new File("testdata/vcf/vcfWithoutGenotypes.vcf");

    // setup the contig ordering
    @BeforeClass
    public static void setupContig() {
        IndexedFastaSequenceFile seq;
        try {
            seq = new IndexedFastaSequenceFile(new File(oneKGLocation + "reference/human_b36_both.fasta"));
        } catch (FileNotFoundException e) {
            throw new StingException("unable to load the sequence dictionary");
        }
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    @Test
    public void testReadBasicHeader() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);

        int seenRecords = 0;

        // check some field entries of each type
        Set<VCFHeaderLine> lines = testSetup.getHeader().getMetaData();

        for (VCFHeaderLine line : lines) {
            // check the vcf info header lines
            if (line instanceof VCFInfoHeaderLine) {
                VCFInfoHeaderLine ihLIne = (VCFInfoHeaderLine)line;

                // test a normal info line
                if (ihLIne.getmName().equals("NS")) {
                    Assert.assertEquals(VCFInfoHeaderLine.INFO_TYPE.Integer,ihLIne.getmType());
                    Assert.assertEquals(1,ihLIne.getmCount());
                    Assert.assertTrue("Number of Samples With Data".equals(ihLIne.getmDescription()));
                    seenRecords++;
                }
                // test a info line that uses the period to represent an unbounded value
                if (ihLIne.getmName().equals("AF")) {
                    Assert.assertEquals(VCFInfoHeaderLine.INFO_TYPE.Float,ihLIne.getmType());
                    Assert.assertEquals(VCFInfoHeaderLine.UNBOUNDED,ihLIne.getmCount());
                    Assert.assertTrue("Allele Frequency".equals(ihLIne.getmDescription()));
                    seenRecords++;
                }
            }
            // check the vcf filter header lines
            if (line instanceof VCFFilterHeaderLine) {
                VCFFilterHeaderLine fhLIne = (VCFFilterHeaderLine)line;
                if (fhLIne.getmName().equals("q10")) {
                    Assert.assertTrue("Quality below 10".equals(fhLIne.getmDescription()));
                    seenRecords++;
                }
            }

            // check the vcf info header lines
            if (line instanceof VCFFormatHeaderLine) {
                VCFFormatHeaderLine ifLIne = (VCFFormatHeaderLine)line;
                if (ifLIne.getmName().equals("GT")) {
                    Assert.assertEquals(VCFFormatHeaderLine.FORMAT_TYPE.String,ifLIne.getmType());
                    Assert.assertEquals(1,ifLIne.getmCount());
                    Assert.assertTrue("Genotype".equals(ifLIne.getmDescription()));
                    seenRecords++;
                }
            }
        }

        Assert.assertEquals("We expected to see three records (one of each type we check), but didn't.",4,seenRecords);
    }

    @Test
    public void testOutputHeader() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);

        File tempFile = null;
        try {
            tempFile = File.createTempFile("VCF4Test","vcf");
        } catch (IOException e) {
            Assert.fail("Couldn't create a temporary file ");
        }

        // write it to disk
        VCFWriter writer = new VCFWriter(tempFile);
        writer.writeHeader(testSetup.getHeader());
        writer.close();

        // md5 sum the file
        Assert.assertTrue("expecting md5sum of 637f3c92806d993a9c7b7116da398d6, but got " + md5SumFile(tempFile),"637f3c92806d993a9c7b7116da398d6".equals(md5SumFile(tempFile)));
    }

    @Test
    public void testCountVCF4Records() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        AsciiLineReader reader = testSetup.getReader();

        // now parse the lines
        String line = null;
        try {
            line = reader.readLine();
        } catch (IOException e) {
            Assert.fail("Failed to read a line");
        }

        // our record count
        int recordCount = 0;
        while (line != null) {
            try {
                //System.err.println(codec.decode(line).toString());
                recordCount++;
                testSetup.codec.decode(line);
                line = reader.readLine();
            } catch (IOException e) {
                Assert.fail("Failed to read a line");
            }
        }
        Assert.assertEquals(7,recordCount);
    }

    @Test
    public void testCountVCF4RecordsWithoutGenotypes() {
        TestSetup testSetup = new TestSetup().invoke(vcfNoGenotypeFile);
        AsciiLineReader reader = testSetup.getReader();

        // now parse the lines
        String line = null;
        try {
            line = reader.readLine();
        } catch (IOException e) {
            Assert.fail("Failed to read a line");
        }

        // our record count
        int recordCount = 0;
        while (line != null) {
            try {
                recordCount++;
                testSetup.codec.decode(line);
                line = reader.readLine();
            } catch (IOException e) {
                Assert.fail("Failed to read a line");
            }
        }
        Assert.assertEquals(6,recordCount);
    }


    // two constants for testing
    String regularLine = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
    String twoFewInfoLine = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
    String twoManyInfoLine = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2;HG=12\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";

    @Test
    public void testCheckInfoValidation() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        testSetup.codec.decode(regularLine);
    }

    @Test
    public void testCheckTwoFewInfoValidation() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        testSetup.codec.decode(twoFewInfoLine);
    }

    @Test(expected=StingException.class)
    public void testCheckTwoManyInfoValidation() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        testSetup.codec.decode(twoManyInfoLine);
    }

    File largeVCF = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/1000GenomesTable1/dindel-v2/CEU.low_coverage.2010_06.indel.genotypes.vcf");

    //@Test
    public void checkLargeVCF() {
        TestSetup testSetup = new TestSetup().invoke(largeVCF);
        AsciiLineReader reader = testSetup.getReader();
        
        // now parse the lines
        String line = null;
        try {
            line = reader.readLine();
        } catch (IOException e) {
            Assert.fail("Failed to read a line");
        }

        // our record count
        int recordCount = 0;
        int badRecordCount = 0;
        long milliseconds = System.currentTimeMillis();
        while (line != null) {
            try {
                recordCount++;
                 try {
                     testSetup.codec.decode(line);
                 } catch (Exception e) {
                     System.err.println(e.getMessage() + " -> " + line);
                     System.err.println(line);
                     badRecordCount++;
                 }
                line = reader.readLine();
                if (recordCount % 1000 == 0)
                    System.err.println("record count == " + recordCount);
            } catch (IOException e) {
                Assert.fail("Failed to read a line");
            }
        }
        System.err.println("Total time = " + (System.currentTimeMillis() - milliseconds));
        Assert.assertEquals(0,badRecordCount);
        Assert.assertEquals(728075,recordCount);
    }

    //@Test
    public void checkBobsCNVVCF() {
        TestSetup testSetup = new TestSetup().invoke(new File("bobs.vcf"));
        AsciiLineReader reader = testSetup.getReader();
        VCF4Codec codec = testSetup.getCodec();

        // now parse the lines
        String line = null;
        try {
            line = reader.readLine();
        } catch (IOException e) {
            Assert.fail("Failed to read a line");
        }

        // our record count
        int recordCount = 0;
        int badRecordCount = 0;
        long milliseconds = System.currentTimeMillis();
        while (line != null) {
            try {
                recordCount++;
                 try {
                     testSetup.codec.decode(line);
                 } catch (Exception e) {
                     System.err.println(e.getMessage() + " -> " + line);
                     System.err.println(line);
                     badRecordCount++;
                 }
                line = reader.readLine();
                if (recordCount % 1000 == 0)
                    System.err.println("record count == " + recordCount);
            } catch (IOException e) {
                Assert.fail("Failed to read a line");
            }
        }
        System.err.println("Total time = " + (System.currentTimeMillis() - milliseconds));
        Assert.assertEquals(0,badRecordCount);
        Assert.assertEquals(15947,recordCount);
    }

    /**
     * test out the clipping of alleles (removing extra context provided by VCF implementation).
     */
    @Test
    public void testClippingOfAllelesDeletionAndInsertion() {
        String ref = "GGTT";
        ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(ref,true));
        alleles.add(Allele.create("GGAATT",false));
        alleles.add(Allele.create("GT",false));

        Pair<GenomeLoc, List<Allele>> locAndList = VCF4Codec.clipAlleles("1",1,ref,alleles);
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",2,3)));

        // we know the ordering
        //System.err.println(locAndList.second.get(0).toString());
        //System.err.println(locAndList.second.get(1).toString());
        //System.err.println(locAndList.second.get(2).toString());
        Assert.assertTrue(locAndList.second.get(0).toString().equals("GT*"));
        Assert.assertTrue(locAndList.second.get(0).isReference());
        Assert.assertTrue(locAndList.second.get(1).toString().equals("GAAT"));
        Assert.assertTrue(locAndList.second.get(2).toString().equals("-"));
    }

    /**
     * test out the clipping of alleles (removing extra context provided by VCF implementation).
     */
    @Test
    public void testClippingManyPotentialFrontClippedBases() {
        String ref = "GGGGTT";
        ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(ref,true));
        alleles.add(Allele.create("GGGGAATT",false));
        alleles.add(Allele.create("GGGT",false));

        Pair<GenomeLoc, List<Allele>> locAndList = VCF4Codec.clipAlleles("1",1,ref,alleles);
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",2,5)));

        // we know the ordering
        //System.err.println(locAndList.second.get(0).toString());
        //System.err.println(locAndList.second.get(1).toString());
        //System.err.println(locAndList.second.get(2).toString());
        Assert.assertTrue(locAndList.second.get(0).toString().equals("GGGT*"));
        Assert.assertTrue(locAndList.second.get(0).isReference());
        Assert.assertTrue(locAndList.second.get(1).toString().equals("GGGAAT"));
        Assert.assertTrue(locAndList.second.get(2).toString().equals("GG"));
    }

    /**
     * test out the clipping of alleles (removing extra context provided by VCF implementation).
     */
    @Test
    public void testClippingOfAllelesLongRefRepeat() {
        String ref = "GGGG";
        ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(ref,true));
        alleles.add(Allele.create("G",false));
        alleles.add(Allele.create("GG",false));

        Pair<GenomeLoc, List<Allele>> locAndList = VCF4Codec.clipAlleles("1",1,ref,alleles);
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",2,4)));

        // we know the ordering
        Assert.assertTrue(locAndList.second.get(0).toString().equals("GGG*"));
        Assert.assertTrue(locAndList.second.get(0).isReference());
        Assert.assertTrue(locAndList.second.get(1).toString().equals("-"));
        Assert.assertTrue(locAndList.second.get(2).toString().equals("G"));
    }

    /**
     * test out the clipping of alleles (removing extra context provided by VCF implementation).
     * TODO - this is kind of a tricky test... we don't know which way clipped the reads, but the position should be accurate
     */
    @Test
    public void testClippingOfAllelesLongRefRepeatClippable() {
        String ref = "GGGGG";
        ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(ref,true));
        alleles.add(Allele.create("GG",false));
        alleles.add(Allele.create("GGG",false));

        Pair<GenomeLoc, List<Allele>> locAndList = VCF4Codec.clipAlleles("1",1,ref,alleles);
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",2,4)));

        // we know the ordering
        Assert.assertTrue(locAndList.second.get(0).toString().equals("GGG*"));
        Assert.assertTrue(locAndList.second.get(0).isReference());
        Assert.assertTrue(locAndList.second.get(1).toString().equals("-"));
        Assert.assertTrue(locAndList.second.get(2).toString().equals("G"));
    }

    /**
     * test out the clipping of alleles (removing extra context provided by VCF implementation).
     */
    @Test
    public void testClippingOfAllelesPlainPolyMorph() {
        String ref = "C";
        ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(ref,true));
        alleles.add(Allele.create("T",false));
        alleles.add(Allele.create("G",false));

        Pair<GenomeLoc, List<Allele>> locAndList = VCF4Codec.clipAlleles("1",1,ref,alleles);
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",1,1)));

        // we know the ordering
        Assert.assertTrue(locAndList.second.get(0).toString().equals("C*"));
        Assert.assertTrue(locAndList.second.get(0).isReference());
        Assert.assertTrue(locAndList.second.get(1).toString().equals("T"));
        Assert.assertTrue(locAndList.second.get(2).toString().equals("G"));
    }

    /**
     * test out the clipping of alleles (removing extra context provided by VCF implementation).
     */
    @Test
    public void testClippingOfAllelesInsertions() {
        String ref = "C";
        ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create(ref,true));
        alleles.add(Allele.create("CTTTTT",false));
        alleles.add(Allele.create("GGGGGG",false));

        Pair<GenomeLoc, List<Allele>> locAndList = VCF4Codec.clipAlleles("1",1,ref,alleles);
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",1,1)));

        // we know the ordering
        Assert.assertTrue(locAndList.second.get(0).toString().equals("C*"));
        Assert.assertTrue(locAndList.second.get(0).isReference());
        Assert.assertTrue(locAndList.second.get(1).toString().equals("CTTTTT"));
        Assert.assertTrue(locAndList.second.get(2).toString().equals("GGGGGG"));
    }

    /**
     * a test setup for the VCF 4 codec
     */
    private class TestSetup {
        private AsciiLineReader reader;
        private VCF4Codec codec;
        private VCFHeader header;

        public AsciiLineReader getReader() {
            return reader;
        }

        public VCF4Codec getCodec() {
            return codec;
        }

        public VCFHeader getHeader() {
            return header;
        }

        public TestSetup invoke(File vcfFile) {
            reader = null;
            try {
                reader = new AsciiLineReader(new FileInputStream(vcfFile));
            } catch (FileNotFoundException e) {
                Assert.fail("Unable to parse out VCF file " + vcfFile);
            }
            codec = new VCF4Codec();
            codec.readHeader(reader);
            header = codec.getHeader(VCFHeader.class);
            return this;
        }
    }
}
