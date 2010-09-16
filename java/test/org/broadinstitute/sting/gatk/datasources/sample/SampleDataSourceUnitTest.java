package org.broadinstitute.sting.gatk.datasources.sample;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Sep 9, 2010
 * Time: 8:21:00 AM
 */
public class SampleDataSourceUnitTest extends BaseTest {

    // this empty header used to instantiate sampledatasource objects
    private static SAMFileHeader header = new SAMFileHeader();

    // all the test sample files are located here
    private String sampleFilesDir = validationDataLocation +  "samples/";

    // make sure samples are created from the SAM file correctly
    @Test()
    public void loadSAMSamplesTest() {
        SampleDataSource s = new SampleDataSource(header, null);
    }

    // tests that a basic sample with relationships loads correctly
    // Note that this is the only test for family relationships - we may want to expand this
    @Test()
    public void basicLoadSampleFileTest() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFile.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
        Assert.assertTrue(s.sampleCount() == 4);
        Sample sampleA = s.getSampleById("sampleA");
        Sample sampleB = s.getSampleById("sampleB");
        Assert.assertTrue(sampleB.getMother() == sampleA);
        Assert.assertTrue(s.getChildren(sampleA).contains(sampleB));
        Set<Sample> family = s.getFamily("family1");
        Assert.assertTrue(family.size() == 2);
        Assert.assertTrue(family.contains(sampleA));
        Assert.assertTrue(family.contains(sampleB));
    }

    // but that file should fail if it has an extra character in it...
    @Test(expected = StingException.class)
    public void loadInvalidSampleExtraCharText() {
        File sampleFile = new File(sampleFilesDir + "invalidSyntaxExtraChar.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
    }

    // ...or a typo...
    @Test(expected = StingException.class)
    public void loadInvalidSampleTypoText() {
        File sampleFile = new File(sampleFilesDir + "invalidSyntaxTypo.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));

    }

    // ...or an extra unrecognized array
    @Test(expected = StingException.class)
    public void loadInvalidSampleExtraArrayText() {
        File sampleFile = new File(sampleFilesDir + "invalidSyntaxExtraArray.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
    }

    // make sure aliases work
    @Test(expected = StingException.class)
    public void sampleAliasText() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFileWithAlias.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
        // this file has two samples, but one has an alias. let's make sure that checks out...
        Assert.assertTrue(s.sampleCount() == 2);
        Assert.assertTrue(s.getSampleById("sampleA") == s.getSampleById("sampleC"));
    }

    // error is thrown if property is included that's not in properties array
    @Test(expected = StingException.class)
    public void unallowedPropertySampleTest() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFileUnallowedProperty.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
    }

    // same as above, with relationship
    @Test(expected = StingException.class)
    public void unallowedRelationshipSampleTest() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFileUnallowedRelationship.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
    }

    // two sample files
    @Test()
    public void twoSampleFilesTest() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFile.yaml");
        File secondFile = new File(sampleFilesDir + "basicSampleFileExt.yaml");
        ArrayList<File> files = new ArrayList<File>();
        files.add(sampleFile);
        files.add(secondFile);
        SampleDataSource s = new SampleDataSource(header, files);
        Assert.assertTrue(s.getSampleById("sampleA").getProperty("propC").equals("valC"));
        Assert.assertTrue(s.getSampleById("sampleA").getProperty("propA").equals("valA"));
    }

    // two sample files, with contradictory properties
    @Test(expected = StingException.class)
    public void twoContradictorySampleFilesTest() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFile.yaml");
        File secondFile = new File(sampleFilesDir + "basicSampleFileInvalidExt.yaml");
        ArrayList<File> files = new ArrayList<File>();
        files.add(sampleFile);
        files.add(secondFile);
        SampleDataSource s = new SampleDataSource(header, files);
    }

    // three sample files
    @Test()
    public void threeSamplesTest() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFile.yaml");
        ArrayList<File> files = new ArrayList<File>();
        files.add(sampleFile);
        files.add(new File(sampleFilesDir + "basicSampleFileExt.yaml"));
        files.add(new File(sampleFilesDir + "basicSampleFileExt2.yaml"));
        SampleDataSource s = new SampleDataSource(header, files);
        Assert.assertTrue(s.sampleCount() == 5);
        Assert.assertTrue(s.getSampleById("sampleE").getProperty("propC").equals("valC"));
        Assert.assertTrue(s.getSampleById("sampleA").getProperty("propA").equals("valA"));
    }                                                                                  

    // make sure we can import data types other than Strings
    @Test()
    public void sampleTestPropertyType() {
        File sampleFile = new File(sampleFilesDir + "sampleFileOtherTypes.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
        Sample sample = s.getSampleById("sampleA");
        Assert.assertTrue(sample.getProperty("a").getClass() == Integer.class);
        Assert.assertTrue(sample.getProperty("b").getClass() == String.class);
        Assert.assertTrue(sample.getProperty("c").getClass() == Double.class);
        Assert.assertTrue(sample.getProperty("b").getClass() == String.class);
    }
    

    // we create lots of single item lists...
    private ArrayList<File> makeFileList(File file) {
        ArrayList<File> a = new ArrayList<File>();
        a.add(file);
        return a;
    }
}
