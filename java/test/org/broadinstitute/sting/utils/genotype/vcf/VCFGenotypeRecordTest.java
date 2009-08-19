package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;

import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class VCFGenotypeRecordTest
 *
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class VCFGenotypeRecordTest extends BaseTest {

    /**
     * test the basic parsing
     */
    @Test
    public void testBasicParsing() {
        String formatString = "GT:B:C:D";
        String genotypeString = "0|1:2:3:4";
        String altAlleles[] = {"A","C","G","T"};
        char referenceBase = 'N';
        VCFGenotypeRecord rec = new VCFGenotypeRecord(formatString,genotypeString,altAlleles,referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.GT_GENOTYPE.PHASED,rec.getPhaseType());
        Assert.assertEquals(referenceBase,rec.getReference());
        Assert.assertEquals("N",rec.getAllele().get(0));
        Assert.assertEquals("A",rec.getAllele().get(1));
        Map<String,String> values = rec.getFields();
        Assert.assertEquals(3,values.size());
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
        VCFGenotypeRecord rec = new VCFGenotypeRecord(formatString,genotypeString,altAlleles,referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.GT_GENOTYPE.PHASED,rec.getPhaseType());
        Assert.assertEquals(referenceBase,rec.getReference());
        Assert.assertEquals("N",rec.getAllele().get(0));
        Assert.assertEquals("A",rec.getAllele().get(1));
        Map<String,String> values = rec.getFields();
        Assert.assertEquals(3,values.size());
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
        VCFGenotypeRecord rec = new VCFGenotypeRecord(formatString,genotypeString,altAlleles,referenceBase);
        Assert.assertEquals(VCFGenotypeRecord.GT_GENOTYPE.PHASED,rec.getPhaseType());
        Assert.assertEquals(referenceBase,rec.getReference());
        Assert.assertEquals("N",rec.getAllele().get(0));
        Assert.assertEquals("A",rec.getAllele().get(1));
        Map<String,String> values = rec.getFields();
        Assert.assertEquals(3,values.size());
        Assert.assertTrue(values.get("B").equals(""));
        Assert.assertTrue(values.get("C").equals(""));
        Assert.assertTrue(values.get("D").equals(""));
    }   
}
