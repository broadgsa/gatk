package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.List;
import java.util.Map;

/** test the VCFReader class test */
public class VCFReaderTest extends BaseTest {

    private static File vcfFile = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample.vcf");
    private static File multiSampleVCF = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/MultiSample.vcf");

    @Test
    public void testVCFInput() {
        VCFReader reader = new VCFReader(vcfFile);
        int counter = 0;
        while (reader.hasNext()) {
            counter++;
            reader.next();
        }
        Assert.assertEquals(5, counter);
    }

    @Test
    public void testMultiSampleVCFInput() {
        VCFReader reader = new VCFReader(multiSampleVCF);
        int counter = 0;
        while (reader.hasNext()) {
            counter++;
            reader.next();
        }
        Assert.assertEquals(99, counter);
    }

    @Test
    public void testNoCallSites() {
        VCFReader reader = new VCFReader(multiSampleVCF);
        if (!reader.hasNext()) Assert.fail("The reader should have a record");
        VCFRecord rec = reader.next();
        final int numberOfNoCallsTruth = 9;
        int noCalls = 0;
        for (VCFGenotypeRecord record : rec.getVCFGenotypeRecords()) {
            List<VCFGenotypeEncoding> encodings = record.getAlleles();
            if (encodings.get(0).getType() == VCFGenotypeEncoding.TYPE.UNCALLED &&
                    encodings.get(1).getType() == VCFGenotypeEncoding.TYPE.UNCALLED)
                noCalls++;
        }
        Assert.assertEquals(numberOfNoCallsTruth, noCalls);
    }


    @Test
    public void testKnownCallSites() {
        VCFReader reader = new VCFReader(multiSampleVCF);
        if (!reader.hasNext()) Assert.fail("The reader should have a record");
        VCFRecord rec = reader.next();
        boolean seenNA11992 = false;
        boolean seenNA12287  = false;
        for (VCFGenotypeRecord record : rec.getVCFGenotypeRecords()) {
            if (record.getSampleName().equals("NA11992")) {
                List<VCFGenotypeEncoding> encodings = record.getAlleles();
                if (!encodings.get(0).getBases().equals("A") ||
                    !encodings.get(1).getBases().equals("A")) {
                    Assert.fail("Sample NA11992 at site 1:10000005 should be AA");
                }
                seenNA11992 = true;
            }
            if (record.getSampleName().equals("NA12287")) {
                List<VCFGenotypeEncoding> encodings = record.getAlleles();
                if (!encodings.get(0).getBases().equals("A") ||
                    !encodings.get(1).getBases().equals("T")) {
                    Assert.fail("Sample NA11992 at site 1:10000005 should be AA");
                }
                seenNA12287 = true;
            }
        }
        Assert.assertTrue(seenNA11992 && seenNA12287);
    }

    /** test the basic parsing */
    @Test
    public void testBasicParsing() {
        String formatString = "GT:B:C:D";
        String genotypeString = "0|1:2:3:4";
        String altAlleles[] = {"A", "G", "T"};
        char referenceBase = 'C';
        VCFGenotypeRecord rec = VCFReader.getVCFGenotype("test", formatString, genotypeString, altAlleles, referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.PHASE.PHASED, rec.getPhaseType());
        Assert.assertEquals("C", rec.getAlleles().get(0).toString());
        Assert.assertEquals("A", rec.getAlleles().get(1).toString());
        Map<String, String> values = rec.getFields();
        Assert.assertEquals(3, values.size());
        Assert.assertTrue(values.get("B").equals("2"));
        Assert.assertTrue(values.get("C").equals("3"));
        Assert.assertTrue(values.get("D").equals("4"));
    }


    /** test the parsing of a genotype field with missing parameters */
    @Test
    public void testMissingFieldParsing() {
        String formatString = "GT:B:C:D";
        String genotypeString = "0|1:::4";
        String altAlleles[] = {"A", "G", "T"};
        char referenceBase = 'C';
        VCFGenotypeRecord rec = VCFReader.getVCFGenotype("test", formatString, genotypeString, altAlleles, referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.PHASE.PHASED, rec.getPhaseType());
        Assert.assertEquals("C", rec.getAlleles().get(0).toString());
        Assert.assertEquals("A", rec.getAlleles().get(1).toString());
        Map<String, String> values = rec.getFields();
        Assert.assertEquals(3, values.size());
        Assert.assertTrue(values.get("B").equals(""));
        Assert.assertTrue(values.get("C").equals(""));
        Assert.assertTrue(values.get("D").equals("4"));
    }

    /** test the parsing of a genotype field with different missing parameters */
    @Test
    public void testMissingAllFields() {
        String formatString = "GT:B:C:D";
        String genotypeString = "0|1:::";
        String altAlleles[] = {"A", "G", "T"};
        char referenceBase = 'C';
        VCFGenotypeRecord rec = VCFReader.getVCFGenotype("test", formatString, genotypeString, altAlleles, referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.PHASE.PHASED, rec.getPhaseType());
        Assert.assertEquals("C", rec.getAlleles().get(0).toString());
        Assert.assertEquals("A", rec.getAlleles().get(1).toString());
        Map<String, String> values = rec.getFields();
        Assert.assertEquals(3, values.size());
        Assert.assertTrue(values.get("B").equals(""));
        Assert.assertTrue(values.get("C").equals(""));
        Assert.assertTrue(values.get("D").equals(""));
    }
}
