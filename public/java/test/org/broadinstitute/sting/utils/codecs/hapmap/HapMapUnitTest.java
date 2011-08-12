/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.hapmap;

import org.broad.tribble.annotation.Strand;
import org.broad.tribble.readers.AsciiLineReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Unit tests for the HapMap codec
 */
public class HapMapUnitTest {
    // our sample hapmap file
    private final static File hapMapFile = new File("public/testdata/genotypes_chr1_ASW_phase3.3_first500.hapmap");
    private final static String knownLine = "rs2185539 C/T chr1 556738 + ncbi_b36 bbs urn:lsid:bbs.hapmap.org:Protocol:Phase3.r3:1 urn:lsid:bbs.hapmap.org:Assay:Phase3.r3_r" +
            "s2185539:1 urn:lsid:dcc.hapmap.org:Panel:US_African-30-trios:4 QC+ CC TC TT CT CC CC CC CC CC CC CC CC CC";
    /**
     * test reading the header off of the file.  We take in the file, read off the first line,
     * close the reader, and then ask the HapMap decoder for the header with a new reader.  These should
     * be equal (i.e. they return the same object).
     */
    @Test
    public void testReadHeader() {
        RawHapMapCodec codec = new RawHapMapCodec();
        AsciiLineReader reader = getReader();
        try {
            String header = reader.readLine();
            reader.close();
            Assert.assertTrue(header.equals(codec.readHeader(getReader())));
        } catch (IOException e) {
            Assert.fail("Unable to read from file " + hapMapFile);
        }
        reader.close();
    }

    @Test
    public void testKnownRecordConversion() {
        RawHapMapCodec codec = new RawHapMapCodec();
        RawHapMapFeature feature = (RawHapMapFeature)codec.decode(knownLine);


        // check that the alleles are right
        Assert.assertEquals(feature.getAlleles().length,2);
        Assert.assertTrue("C".equals(feature.getAlleles()[0]));
        Assert.assertTrue("T".equals(feature.getAlleles()[1]));

        // check the name
        Assert.assertTrue("rs2185539".equals(feature.getName()));

        // check the position
        Assert.assertEquals(feature.getStart(),556738);
        Assert.assertEquals(feature.getEnd(),556738);

        // check the contig
        Assert.assertTrue("chr1".equals(feature.getChr()));
                
        // check the assembly, center, protLSID, assayLSID, panelLSID, and qccode
        Assert.assertTrue("ncbi_b36".equals(feature.getAssembly()));
        Assert.assertTrue("bbs".equals(feature.getCenter()));
        Assert.assertTrue("urn:lsid:bbs.hapmap.org:Protocol:Phase3.r3:1".equals(feature.getProtLSID()));
        Assert.assertTrue("urn:lsid:bbs.hapmap.org:Assay:Phase3.r3_rs2185539:1".equals(feature.getAssayLSID()));
        Assert.assertTrue("urn:lsid:dcc.hapmap.org:Panel:US_African-30-trios:4".equals(feature.getPanelLSID()));
        Assert.assertTrue("QC+".equals(feature.getQCCode()));

        // check the strand
        Assert.assertEquals(feature.getStrand(),Strand.POSITIVE);

        // check the genotypes
        int x = 0;
        for (; x < feature.getGenotypes().length; x++) {
            switch (x) {
                case 1: Assert.assertTrue("TC".equals(feature.getGenotypes()[x])); break;
                case 2: Assert.assertTrue("TT".equals(feature.getGenotypes()[x])); break;
                case 3: Assert.assertTrue("CT".equals(feature.getGenotypes()[x])); break;
                default: Assert.assertTrue("CC".equals(feature.getGenotypes()[x])); break;
            }
        }
        // assert that we found the correct number of records
        Assert.assertEquals(x,13);
    }

    @Test
    public void testReadCorrectNumberOfRecords() {
        // setup the record for reading our 500 line file (499 records, 1 header line)
        RawHapMapCodec codec = new RawHapMapCodec();
        AsciiLineReader reader = getReader();

        String line;
        int count = 0;
        try {
            codec.readHeader(reader);
            line = reader.readLine();
            while (line != null) {
                codec.decode(line);
                line = reader.readLine();
                ++count;
            }
        } catch (IOException e) {
            Assert.fail("IOException " + e.getMessage());
        }
        reader.close();
        Assert.assertEquals(count,499);
    }

    @Test
    public void testGetSampleNames() {
        // setup the record for reading our 500 line file (499 records, 1 header line)
        RawHapMapCodec codec = new RawHapMapCodec();
        AsciiLineReader reader = getReader();

        String line;
        try {
            codec.readHeader(reader);
            line = reader.readLine();
            RawHapMapFeature feature = (RawHapMapFeature) codec.decode(line);
            Assert.assertEquals(feature.getSampleIDs().length,87);

        } catch (IOException e) {
            Assert.fail("IOException " + e.getMessage());
        }
        reader.close();
    }


    public AsciiLineReader getReader() {
        try {
            return new AsciiLineReader(new FileInputStream(hapMapFile));            
        } catch (FileNotFoundException e) {
            Assert.fail("Unable to open hapmap file : " + hapMapFile);
        }
        return null; // for intellij, it doesn't know that assert.fail is fatal
    }

}
