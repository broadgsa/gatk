package org.broadinstitute.sting.utils.genotype.vcf;

import org.junit.Test;
import org.junit.Assert;
import org.broadinstitute.sting.BaseTest;

import java.io.File;
import java.util.Map;

/**
 * test the VCFReader class test
 */
public class VCFReaderTest extends BaseTest {

    private static File vcfFile = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample.vcf");

    @Test
    public void testVCFInput() {
        VCFReader reader = new VCFReader(vcfFile);
        int counter = 0;
        while (reader.hasNext()) {
            counter++;
            reader.next();
        }
        Assert.assertEquals(5,counter);
    }

     /**
     * test the basic parsing
     */
    @Test
    public void testBasicParsing() {
        String formatString = "GT:B:C:D";
        String genotypeString = "0|1:2:3:4";
        String altAlleles[] = {"A","C","G","T"};
        char referenceBase = 'N';
        VCFGenotypeRecord rec = VCFReader.getVCFGenotype("test",formatString,genotypeString,altAlleles,referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.PHASE.PHASED,rec.getPhaseType());
        Assert.assertEquals(referenceBase,rec.getReference());
        Assert.assertEquals("N",rec.getAllele().get(0));
        Assert.assertEquals("A",rec.getAllele().get(1));
        Map<String,String> values = rec.getFields();
        Assert.assertEquals(4,values.size());
        Assert.assertTrue(values.get("B").equals("2"));
        Assert.assertTrue(values.get("C").equals("3"));
        Assert.assertTrue(values.get("D").equals("4"));
    }


    /**
     * test the parsing of a genotype field with missing parameters
     */
    @Test
    public void testMissingFieldParsing() {
        String formatString = "GT:B:C:D";
        String genotypeString = "0|1:::4";
        String altAlleles[] = {"A","C","G","T"};
        char referenceBase = 'N';
        VCFGenotypeRecord rec = VCFReader.getVCFGenotype("test",formatString,genotypeString,altAlleles,referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.PHASE.PHASED,rec.getPhaseType());
        Assert.assertEquals(referenceBase,rec.getReference());
        Assert.assertEquals("N",rec.getAllele().get(0));
        Assert.assertEquals("A",rec.getAllele().get(1));
        Map<String,String> values = rec.getFields();
        Assert.assertEquals(4,values.size());
        Assert.assertTrue(values.get("B").equals(""));
        Assert.assertTrue(values.get("C").equals(""));
        Assert.assertTrue(values.get("D").equals("4"));
    }

    /**
     * test the parsing of a genotype field with different missing parameters
     */
    @Test
    public void testMissingAllFields() {
        String formatString = "GT:B:C:D";
        String genotypeString = "0|1:::";
        String altAlleles[] = {"A","C","G","T"};
        char referenceBase = 'N';
        VCFGenotypeRecord rec = VCFReader.getVCFGenotype("test",formatString,genotypeString,altAlleles,referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.PHASE.PHASED,rec.getPhaseType());
        Assert.assertEquals(referenceBase,rec.getReference());
        Assert.assertEquals("N",rec.getAllele().get(0));
        Assert.assertEquals("A",rec.getAllele().get(1));
        Map<String,String> values = rec.getFields();
        Assert.assertEquals(4,values.size());
        Assert.assertTrue(values.get("B").equals(""));
        Assert.assertTrue(values.get("C").equals(""));
        Assert.assertTrue(values.get("D").equals(""));
    }
}
