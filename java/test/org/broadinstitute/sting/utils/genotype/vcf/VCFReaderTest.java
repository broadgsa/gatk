package org.broadinstitute.sting.utils.genotype.vcf;

import org.junit.Test;
import org.junit.Assert;
import org.broadinstitute.sting.BaseTest;

import java.io.File;

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


}
