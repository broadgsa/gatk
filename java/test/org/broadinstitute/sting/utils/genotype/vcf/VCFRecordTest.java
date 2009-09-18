package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * 
 * @author aaron 
 * 
 * Class VCFRecordTest
 *
 * test the basic functionality of the vcf record
 */
public class VCFRecordTest extends BaseTest {

    private VCFRecord makeFakeVCFRecord() {
        List<String> altBases = new ArrayList<String>();
        altBases.add("C");
        altBases.add("D1");
        Map<String,String> infoFields = new HashMap<String,String>();
        infoFields.put("DP","50");
        List<VCFGenotypeRecord> genotypeObjects = new ArrayList<VCFGenotypeRecord>();
        Map<String, String> keyValues = new HashMap<String,String>();
        keyValues.put("AA","2");
        List<String> Alleles = new ArrayList<String>();
        Alleles.add("A");
        genotypeObjects.add(new VCFGenotypeRecord("SampleName", Alleles, VCFGenotypeRecord.PHASE.PHASED, keyValues));
        return new VCFRecord('A',"chr1",1,"RANDOM",altBases,0,".",infoFields, "GT:AA",genotypeObjects);
    }


    @Test
    public void testGetGenotypes() {
        VCFRecord rec = makeFakeVCFRecord();
        List<VCFGenotypeRecord> genotypeObjects = rec.getVCFGenotypeRecords();
        Assert.assertEquals(1,genotypeObjects.size());
        Assert.assertTrue(genotypeObjects.get(0).getSampleName().equals("SampleName"));
    }
}
