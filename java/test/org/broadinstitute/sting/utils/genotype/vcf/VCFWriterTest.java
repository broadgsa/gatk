package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class VCFWriterTest
 *         <p/>
 *         This class tests out the ability of the VCF writer to correctly write VCF files
 */
public class VCFWriterTest extends BaseTest {
    private Set<VCFHeader.HEADER_FIELDS> headerFields = new LinkedHashSet<VCFHeader.HEADER_FIELDS>();
    private Map<String, String> metaData = new HashMap();
    private List<String> additionalColumns = new ArrayList<String>();
    private File fakeVCFFile = new File("FAKEVCFFILEFORTESTING.vcf");

    /** test, using the writer and reader, that we can output and input a VCF file without problems */
    @Test
    public void testBasicWriteAndRead() {
        VCFHeader header = createFakeHeader();
        VCFWriter writer = new VCFWriter(header,fakeVCFFile);
        writer.addRecord(createVCFRecord(header));
        writer.addRecord(createVCFRecord(header));
        writer.close();
        VCFReader reader = new VCFReader(fakeVCFFile);
        int counter = 0;
        // validate what we're reading in
        validateHeader(reader.getHeader());
        for(VCFRecord rec :reader) {
            counter++;
        }
        Assert.assertEquals(2,counter);        
        reader.close();
        fakeVCFFile.delete();
    }

    /**
     * create a fake header of known quantity
     * @return a fake VCF header
     */
    private VCFHeader createFakeHeader() {
        for (VCFHeader.HEADER_FIELDS field : VCFHeader.HEADER_FIELDS.values()) {
            headerFields.add(field);
        }
        metaData.put("format", "VCFv3.2"); // required 
        metaData.put("two", "2");
        additionalColumns.add("FORMAT");
        additionalColumns.add("extra1");
        additionalColumns.add("extra2");
        // this should create a header that is valid

        return new VCFHeader(headerFields, metaData, additionalColumns);
    }

    private VCFRecord createVCFRecord(VCFHeader header) {
        int totalVals = header.getColumnCount();
        List<String> array = new ArrayList<String>();
        for (int x = 0; x < totalVals; x++) {
            array.add(String.valueOf(x));    
        }
        return new VCFRecord(header,array);
    }


    /**
     * validate a VCF header
     * @param header the header to validate
     */
    public void validateHeader(VCFHeader header) {
        // check the fields
        int index = 0;
        for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) {
            Assert.assertEquals(VCFHeader.HEADER_FIELDS.values()[index], field);
            index++;
        }
        index = 0;
        for (String key : header.getMetaData().keySet()) {
            Assert.assertEquals(header.getMetaData().get(key), metaData.get(key));
            index++;
        }
        Assert.assertEquals(metaData.size(), index);
        index = 0;
        for (String key : header.getGenotypeSamples()) {
            Assert.assertTrue(additionalColumns.contains(key));
            index++;
        }
        Assert.assertEquals(additionalColumns.size(), index+1 /* for the header field we don't see */);
    }
}
