package org.broadinstitute.sting.gatk.refdata.features.vcf4;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.collections.Pair;
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
import java.util.Map;
import java.util.Set;

import net.sf.picard.reference.IndexedFastaSequenceFile;

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
        seq = new IndexedFastaSequenceFile(new File(oneKGLocation + "reference/human_b36_both.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq.getSequenceDictionary());
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
                if (ihLIne.getName().equals("NS")) {
                    Assert.assertEquals(VCFHeaderLineType.Integer,ihLIne.getType());
                    Assert.assertEquals(1,ihLIne.getCount());
                    Assert.assertTrue("Number of Samples With Data".equals(ihLIne.getDescription()));
                    seenRecords++;
                }
                // test a info line that uses the period to represent an unbounded value
                if (ihLIne.getName().equals("AF")) {
                    Assert.assertEquals(VCFHeaderLineType.Float,ihLIne.getType());
                    Assert.assertEquals(VCFInfoHeaderLine.UNBOUNDED,ihLIne.getCount());
                    Assert.assertTrue("Allele Frequency".equals(ihLIne.getDescription()));
                    seenRecords++;
                }
            }
            // check the vcf filter header lines
            if (line instanceof VCFFilterHeaderLine) {
                VCFFilterHeaderLine fhLIne = (VCFFilterHeaderLine)line;
                if (fhLIne.getName().equals("q10")) {
                    Assert.assertTrue("Quality below 10".equals(fhLIne.getDescription()));
                    seenRecords++;
                }
            }

            // check the vcf info header lines
            if (line instanceof VCFFormatHeaderLine) {
                VCFFormatHeaderLine ifLIne = (VCFFormatHeaderLine)line;
                if (ifLIne.getName().equals("GT")) {
                    Assert.assertEquals(VCFHeaderLineType.String,ifLIne.getType());
                    Assert.assertEquals(1,ifLIne.getCount());
                    Assert.assertTrue("Genotype".equals(ifLIne.getDescription()));
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
            tempFile.deleteOnExit();
        } catch (IOException e) {
            Assert.fail("Couldn't create a temporary file ");
        }
        // write it to disk
        VCFWriter writer = new VCFWriter(tempFile);
        writer.writeHeader(testSetup.getHeader());
        writer.close();

        // md5 sum the file
        // TODO -- uncomment this when we have a better solution than using md5s in a unit test
        //Assert.assertTrue("expecting md5sum of e376c7cb1831d3cbdca670f360b7f022, but got " + md5SumFile(tempFile),"e376c7cb1831d3cbdca670f360b7f022".equals(md5SumFile(tempFile)));
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


    // test too many info fields - NOT a valid test with validation turned off in the VCF4 reader
    String twoManyInfoLine = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2;HH\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:0,0";
    //@Test(expected=StingException.class)
    public void testCheckTooManyInfoFields() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        testSetup.codec.decode(twoManyInfoLine);
    }
    // test a regular line
    String regularLine = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:0,0";
    @Test
    public void testCheckInfoValidation() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        testSetup.codec.decode(regularLine);
    }
    // test too few info lines, we don't provide the DP in this line
    // test GT field in the incorrect position (!= 0)
    String GTFieldInTheWrongPosition = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;AF=0.5;DB;H2\tGQ:DP:HQ:GT\t48:1:51,51:0|0\t48:8:51,51:0|0\t43:5:0,0:0|0";
    @Test(expected=RuntimeException.class)
    public void testCheckGTFieldOrdering() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        testSetup.codec.decode(GTFieldInTheWrongPosition);
    }
    
    // test too few info lines, we don't provide the DP in this line
    String twoFewInfoLine = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;AF=0.5;DB\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t0|0:48:1:51,51\t0|0:48:1:51,51";
    @Test
    public void testCheckTwoFewInfoValidation() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        testSetup.codec.decode(twoFewInfoLine);
    }


    
    // test that we're getting the right genotype for a multi-base polymorphism
    String MNPLine = "20\t14370\trs6054257\tGG\tAT\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
    @Test
    public void testMNPValidation() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        VariantContext vc = (VariantContext)testSetup.codec.decode(MNPLine);
        Map<String, Genotype> genotypes = vc.getGenotypes();
        Assert.assertTrue(genotypes.containsKey("NA00003"));
        Genotype g = genotypes.get("NA00003");
        Assert.assertTrue("Expected AT genotype, saw = " + g.getAllele(0),"AT".equals(g.getAllele(0).toString()));
        Assert.assertTrue(vc.getType()== VariantContext.Type.MNP);
    }

    // test that we're getting the right genotype for what appears to be a multi-base polymorphism, but is really just a SNP
    String MNPLine2 = "20\t14370\trs6054257\tGT\tAT\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
    @Test
    public void testMNPWannabeButReallyASNPValidation() {
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        VariantContext vc = (VariantContext)testSetup.codec.decode(MNPLine2);
        Map<String, Genotype> genotypes = vc.getGenotypes();
        Assert.assertTrue(genotypes.containsKey("NA00003"));
        Genotype g = genotypes.get("NA00003");
        Assert.assertTrue("Expected A genotype, saw = " + g.getAllele(0),"A".equals(g.getAllele(0).toString()));
        Assert.assertTrue(vc.getType()== VariantContext.Type.SNP);
    }

    File largeVCF = new File("yri.vcf"); // change to whatever file you'd like to test in the following test

    // @Test uncomment to re-enable testing
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
                     //System.err.println(e.getMessage() + " -> " + line);
                     //System.err.println(line);
                     Assert.fail("Bad record from line " + line + " message = " + e.getMessage());
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
                     Assert.fail("Bad record from line " + line + " message = " + e.getMessage());
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
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",1,3)));

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
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",1,5)));

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
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",1,4)));

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
        Assert.assertTrue(locAndList.first.equals(GenomeLocParser.createGenomeLoc("1",1,4)));

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

    @Test
    public void testGenotypeConversionPhasing() {
        String[] parts = {"GT:GD:DP", "0|0", "0|1", "1\\1"};
        List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(Allele.create("A", true));
        alleles.add(Allele.create("G", false));
        Pair<GenomeLoc, List<Allele>> locAndAlleles = new Pair<GenomeLoc, List<Allele>>(GenomeLocParser.createGenomeLoc("1",1),alleles);
        TestSetup testSetup = new TestSetup().invoke(vcfGenotypeFile);
        Map<String, Genotype> genotypes = testSetup.getCodec().createGenotypeMap(parts, locAndAlleles,0);
        // assert the first genotype is phased, and the third is not
        Assert.assertTrue(genotypes.get("NA00001").genotypesArePhased());
        Assert.assertTrue(!genotypes.get("NA00003").genotypesArePhased());
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
            header = (VCFHeader)codec.readHeader(reader);            
            return this;
        }
    }
}
