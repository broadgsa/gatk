package org.broadinstitute.sting.gatk.datasources.sample;

import net.sf.samtools.SAMFileHeader;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.StingException;

import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

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
        Assert.assertTrue(s.sampleCount() == 5);
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
    @Test(expectedExceptions=StingException.class)
    public void loadInvalidSampleExtraCharText() {
        File sampleFile = new File(sampleFilesDir + "invalidSyntaxExtraChar.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
    }

    // ...or a typo...
    @Test(expectedExceptions=StingException.class)
    public void loadInvalidSampleTypoText() {
        File sampleFile = new File(sampleFilesDir + "invalidSyntaxTypo.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));

    }

    // ...or an extra unrecognized array
    @Test(expectedExceptions=StingException.class)
    public void loadInvalidSampleExtraArrayText() {
        File sampleFile = new File(sampleFilesDir + "invalidSyntaxExtraArray.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
    }

    // make sure aliases work
    @Test(expectedExceptions=StingException.class)
    public void sampleAliasText() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFileWithAlias.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
        // this file has two samples, but one has an alias. let's make sure that checks out...
        Assert.assertTrue(s.sampleCount() == 3);
        Assert.assertTrue(s.getSampleById("sampleA") == s.getSampleById("sampleC"));
    }

    // error is thrown if property is included that's not in properties array
    @Test(expectedExceptions=StingException.class)
    public void unallowedPropertySampleTest() {
        File sampleFile = new File(sampleFilesDir + "basicSampleFileUnallowedProperty.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
    }

    // same as above, with relationship
    @Test(expectedExceptions=StingException.class)
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
    @Test(expectedExceptions=StingException.class)
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
        Assert.assertTrue(s.sampleCount() == 6);
        Assert.assertTrue(s.getSampleById("sampleE").getProperty("propC").equals("valC"));
        Assert.assertTrue(s.getSampleById("sampleA").getProperty("propA").equals("valA"));
    }

    /**
     * testing getSamplesWithProperty
     * in this file there are 5 samples - 2 with population "CEU", 1 with population "ABC", 1 with no population,
     * and then the default null sample
     */
    @Test()
    public void getSamplesWithPropertyTest() {
        File sampleFile = new File(sampleFilesDir + "sampleFileWithProperties.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
        Assert.assertTrue(s.sampleCount() == 5);
        Set<Sample> ceuSamples = s.getSamplesWithProperty("population", "CEU");
        Assert.assertTrue(ceuSamples.size() == 2);

        Iterator<Sample> i = ceuSamples.iterator();
        ArrayList<String> sampleNames = new ArrayList<String>();
        sampleNames.add(i.next().getId());
        sampleNames.add(i.next().getId());
        Assert.assertTrue(sampleNames.contains("sampleA"));
        Assert.assertTrue(sampleNames.contains("sampleB"));
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

    /**
     * check that getSamplesFromVariantContext works
     * create a variant context with two sample names, and make sure the right samples are there
     */
    @Test()
    public void variantContextTest() {
        SampleDataSource s = new SampleDataSource(header, null);
        List<Allele> alleleCollection = new ArrayList<Allele>();
        Allele a1 = Allele.create("A", true);
        alleleCollection.add(a1);

        Set<Genotype> genotypeCollection = new HashSet<Genotype>();
        genotypeCollection.add(new Genotype("NA123", alleleCollection));
        genotypeCollection.add(new Genotype("NA456", alleleCollection));

        VariantContext v = new VariantContext("contextName", "chr1", 1, 1, alleleCollection, genotypeCollection);

        // make sure the set that's returned is the right size
        HashSet<Sample> set = (HashSet) s.getSamplesByVariantContext(v);
        Assert.assertTrue(set.size() == 2);

        // make sure both samples are included
        Iterator<Sample> i = set.iterator();
        ArrayList<String> sampleNames = new ArrayList<String>();
        sampleNames.add(i.next().getId());
        sampleNames.add(i.next().getId());
        Assert.assertTrue(sampleNames.contains("NA123"));
        Assert.assertTrue(sampleNames.contains("NA456"));
    }

    /**
     * checking subContextFromSampleProperty
     */

    /**
     * check that subContextFromSampleProperty works
     * create a variant context with four sample names, make sure that it filters correctly to 2
     */
    @Test()
    public void subContextFromSamplePropertyTest() {

        File sampleFile = new File(sampleFilesDir + "sampleFileWithProperties.yaml");
        SampleDataSource s = new SampleDataSource(header, makeFileList(sampleFile));
        Assert.assertTrue(s.sampleCount() == 5);

        List<Allele> alleleCollection = new ArrayList<Allele>();
        Allele a1 = Allele.create("A", true);
        alleleCollection.add(a1);

        Set<Genotype> genotypeCollection = new HashSet<Genotype>();
        genotypeCollection.add(new Genotype("NA123", alleleCollection));
        genotypeCollection.add(new Genotype("sampleA", alleleCollection));
        genotypeCollection.add(new Genotype("sampleB", alleleCollection));
        genotypeCollection.add(new Genotype("sampleC", alleleCollection));

        VariantContext v = new VariantContext("contextName", "chr1", 1, 1, alleleCollection, genotypeCollection);
        VariantContext subContext = s.subContextFromSampleProperty(v, "population", "CEU");

        Assert.assertTrue(subContext.getSampleNames().contains("sampleA"));
        Assert.assertTrue(subContext.getSampleNames().contains("sampleA"));
        Assert.assertTrue(subContext.getSampleNames().size() == 2);

    }
    

    // we create lots of single item lists...
    private ArrayList<File> makeFileList(File file) {
        ArrayList<File> a = new ArrayList<File>();
        a.add(file);
        return a;
    }
}
