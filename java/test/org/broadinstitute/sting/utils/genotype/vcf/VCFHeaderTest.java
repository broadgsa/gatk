package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;


/**
 * 
 * @author aaron 
 * 
 * Class VCFHeaderTest
 *
 * Test the VCF Header class
 */
public class VCFHeaderTest extends BaseTest {

    private Set<VCFHeader.HEADER_FIELDS> headerFields = new LinkedHashSet<VCFHeader.HEADER_FIELDS>();
    private Set<String> metaData = new HashSet();
    private Set<String> additionalColumns = new HashSet<String>();

    /**
     * give it fake data, and make sure we get back the right fake data
     */
    @Test
    public void testHeaderConstructor() {
        metaData.add(VCFHeader.FULL_FORMAT_LINE); // required
        metaData.add("two=2");
        additionalColumns.add("extra1");
        additionalColumns.add("extra2");
        // this should create a header that is valid

        VCFHeader header = new VCFHeader(metaData, additionalColumns);

        // check the fields
        int index = 0;
        for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) {
            Assert.assertEquals(VCFHeader.HEADER_FIELDS.values()[index],field);
            index++;
        }
        Assert.assertEquals(metaData.size(), header.getMetaData().size());
        index = 0;
        for (String key: header.getGenotypeSamples()) {
            Assert.assertTrue(additionalColumns.contains(key));
            index++;
        }
        Assert.assertEquals(additionalColumns.size(),index);
    }

    
}
