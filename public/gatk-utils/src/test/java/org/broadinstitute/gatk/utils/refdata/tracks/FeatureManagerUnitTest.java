/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.refdata.tracks;


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.codecs.table.BedTableCodec;
import org.broadinstitute.gatk.utils.codecs.table.TableFeature;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;


/**
 * @author depristo
 *
 * UnitTests for RMD FeatureManager
 */
public class FeatureManagerUnitTest extends BaseTest {
    private static final File RANDOM_FILE = new File(publicTestDir+ "exampleGATKReport.eval");
    private static final File VCF3_FILE = new File(privateTestDir + "vcf3.vcf");
    private static final File VCF4_FILE = new File(privateTestDir + "HiSeq.10000.vcf");
    private static final File VCF4_FILE_GZ = new File(privateTestDir + "HiSeq.10000.vcf.gz");
    private static final File VCF4_FILE_BGZIP = new File(privateTestDir + "HiSeq.10000.bgzip.vcf.gz");

    private FeatureManager manager;
    private GenomeLocParser genomeLocParser;

    @BeforeMethod
    public void setup() {
        File referenceFile = new File(b36KGReference);
        try {
            IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(referenceFile);
            genomeLocParser = new GenomeLocParser(seq);
            manager = new FeatureManager();
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
    }

    @Test
    public void testManagerCreation() {
        Assert.assertTrue(manager.getFeatureDescriptors().size() > 0);
    }

    private class FMTest extends BaseTest.TestDataProvider {
        public Class codec;
        public Class<? extends Feature> feature;
        public String name;
        public File associatedFile;

        private FMTest(final Class feature, final Class codec, final String name, final File file) {
            super(FMTest.class);
            this.codec = codec;
            this.feature = feature;
            this.name = name;
            this.associatedFile = file;
        }

        public void assertExpected(FeatureManager.FeatureDescriptor featureDescriptor) {
            Assert.assertEquals(featureDescriptor.getCodecClass(), codec);
            Assert.assertEquals(featureDescriptor.getFeatureClass(), feature);
            Assert.assertEquals(featureDescriptor.getName().toLowerCase(), name.toLowerCase());
        }

        public String toString() {
            return String.format("FMTest name=%s codec=%s feature=%s file=%s",
                    name, codec.getSimpleName(), feature.getSimpleName(), associatedFile);
        }
    }

    @DataProvider(name = "tests")
    public Object[][] createTests() {
        new FMTest(VariantContext.class, VCF3Codec.class, "VCF3", VCF3_FILE);
        new FMTest(VariantContext.class, VCFCodec.class, "VCF", VCF4_FILE);
        new FMTest(VariantContext.class, VCFCodec.class, "VCF", VCF4_FILE_GZ);
        new FMTest(VariantContext.class, VCFCodec.class, "VCF", VCF4_FILE_BGZIP);
        new FMTest(TableFeature.class, BedTableCodec.class, "bedtable", null);
        return FMTest.getTests(FMTest.class);
    }

    @Test(dataProvider = "tests")
    public void testGetByFile(FMTest params) {
        if ( params.associatedFile != null ) {
            FeatureManager.FeatureDescriptor byFile = manager.getByFiletype(params.associatedFile);
            Assert.assertNotNull(byFile, "Couldn't find any type associated with file " + params.associatedFile);
            params.assertExpected(byFile);
        }
    }

    @Test
    public void testGetByFileNoMatch() {
        FeatureManager.FeatureDescriptor byFile = manager.getByFiletype(RANDOM_FILE);
        Assert.assertNull(byFile, "Found type " + byFile + " associated with RANDOM, non-Tribble file " + RANDOM_FILE);
    }

    @Test(dataProvider = "tests")
    public void testGetters(FMTest params) {
        params.assertExpected(manager.getByCodec(params.codec));
        params.assertExpected(manager.getByName(params.name));
        params.assertExpected(manager.getByName(params.name.toLowerCase()));
        params.assertExpected(manager.getByName(params.name.toUpperCase()));

        Collection<FeatureManager.FeatureDescriptor> descriptors = manager.getByFeature(params.feature);
        Assert.assertTrue(descriptors.size() > 0, "Look up by FeatureClass failed");
    }

    @Test
    public void testUserFriendlyList() {
        Assert.assertTrue(manager.userFriendlyListOfAvailableFeatures().length() > 0, "Expected at least one codec to be listed");
        Assert.assertTrue(manager.userFriendlyListOfAvailableFeatures().split(",").length > 0, "Expected at least two codecs, but only saw one");
    }

    @Test
    public void testCodecCreation() {
        FeatureManager.FeatureDescriptor descriptor = manager.getByName("vcf");
        Assert.assertNotNull(descriptor, "Couldn't find VCF feature descriptor!");

        FeatureCodec c = manager.createCodec(descriptor, "foo", genomeLocParser, null);
        Assert.assertNotNull(c, "Couldn't create codec");
        Assert.assertEquals(c.getClass(), descriptor.getCodecClass());
        Assert.assertEquals(c.getFeatureType(), descriptor.getFeatureClass());
    }

}

