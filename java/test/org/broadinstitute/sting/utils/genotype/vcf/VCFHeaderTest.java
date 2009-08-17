package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.StingException;
import org.junit.Test;
import org.junit.Assert;

import java.io.*;
import java.util.ArrayList;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class VCFHeaderTest
 *
 * Test the VCF Header class
 */
public class VCFHeaderTest extends BaseTest {

    @Test
    public void test1() {
        File in = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample.vcf");
        if (!in.exists()) throw new StingException("vfc doesn't exist");
        List<String> array = new ArrayList<String>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader("vcfexample.vcf"));
            String line = reader.readLine();
            while (line.startsWith("#")) {
                array.add(line);
                line = reader.readLine();
            }
            VCFHeader header = new VCFHeader(array);
        } catch (FileNotFoundException e) {
            Assert.fail("File not found exception in VCFHeaderTest");
        } catch (IOException e) {
            Assert.fail("IO exception in VCFHeaderTest");
        }
    }
    
}
