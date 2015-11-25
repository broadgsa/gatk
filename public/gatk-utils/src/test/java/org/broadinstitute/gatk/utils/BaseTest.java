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

package org.broadinstitute.gatk.utils;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.log4j.AppenderSkeleton;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.apache.log4j.spi.LoggingEvent;
import org.broadinstitute.gatk.utils.variant.VCIterable;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.commandline.CommandLineUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.SkipException;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 *
 * User: aaron
 * Date: Apr 14, 2009
 * Time: 10:24:30 AM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 * <p/>
 * Class BaseTest
 * <p/>
 * This is the base test class for all of our test cases.  All test cases should extend from this
 * class; it sets up the logger, and resolves the location of directories that we rely on.
 */
@SuppressWarnings("unchecked")
public abstract class BaseTest {
    /** our log, which we want to capture anything from org.broadinstitute.sting */
    public static final Logger logger = CommandLineUtils.getStingLogger();

    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");
    public static final String gatkDirectory = System.getProperty("gatkdir", CURRENT_DIRECTORY) + "/";
    public static final String baseDirectory = System.getProperty("basedir", CURRENT_DIRECTORY) + "/";
    public static final String testType = System.getProperty("testType"); // May be null
    public static final String testTypeSubDirectory = testType == null ? "" : ("/" + testType); // May be empty

    public static final String hg18Reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    public static final String hg19Reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
    public static final String b36KGReference = "/humgen/1kg/reference/human_b36_both.fasta";
    public static final String b37KGReference = "/humgen/1kg/reference/human_g1k_v37.fasta";
    public static final String b37KGReferenceWithDecoy = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37_decoy.fasta";
    public static final String hg19ReferenceWithChrPrefixInChromosomeNames = "/humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta";
    public static final String GATKDataLocation = "/humgen/gsa-hpprojects/GATK/data/";
    public static final String validationDataLocation = GATKDataLocation + "Validation_Data/";
    public static final String evaluationDataLocation = GATKDataLocation + "Evaluation_Data/";
    public static final String comparisonDataLocation = GATKDataLocation + "Comparisons/";
    public static final String annotationDataLocation = GATKDataLocation + "Annotations/";

    public static final String b37GoodBAM = validationDataLocation + "/CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    public static final String b37GoodNA12878BAM = validationDataLocation + "/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam";
    public static final String b37_NA12878_OMNI = validationDataLocation + "/NA12878.omni.vcf";

    public static final String dbsnpDataLocation = GATKDataLocation;
    public static final String b36dbSNP129 = dbsnpDataLocation + "dbsnp_129_b36.vcf";
    public static final String b37dbSNP129 = dbsnpDataLocation + "dbsnp_129_b37.vcf";
    public static final String b37dbSNP132 = dbsnpDataLocation + "dbsnp_132_b37.vcf";
    public static final String b37dbSNP138 = "/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_138.b37.vcf";
    public static final String hg18dbSNP132 = dbsnpDataLocation + "dbsnp_132.hg18.vcf";

    public static final String hapmapDataLocation = comparisonDataLocation + "Validated/HapMap/3.3/";
    public static final String b37hapmapGenotypes = hapmapDataLocation + "genotypes_r27_nr.b37_fwd.vcf";

    public static final String intervalsLocation = "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/";
    public static final String hg19Intervals = intervalsLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list";
    public static final String hg19Chr20Intervals = GATKDataLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.chr20.interval_list";

    public static final boolean REQUIRE_NETWORK_CONNECTION = false;
    private static final String networkTempDirRoot = "/broad/hptmp/";
    private static final boolean networkTempDirRootExists = new File(networkTempDirRoot).exists();
    private static final File networkTempDirFile;

    private static final String privateTestDirRelative = "private/gatk-tools-private/src/test/resources/";
    public static final String privateTestDir = new File(gatkDirectory, privateTestDirRelative).getAbsolutePath() + "/";
    protected static final String privateTestDirRoot = privateTestDir.replace(privateTestDirRelative, "");

    private static final String publicTestDirRelative = "public/gatk-utils/src/test/resources/";
    public static final String publicTestDir = new File(gatkDirectory, publicTestDirRelative).getAbsolutePath() + "/";
    protected static final String publicTestDirRoot = publicTestDir.replace(publicTestDirRelative, "");

    public static final String keysDataLocation = validationDataLocation + "keys/";

    public static final String exampleFASTA = publicTestDir + "exampleFASTA.fasta";

    public final static String NA12878_PCRFREE = privateTestDir + "PCRFree.2x250.Illumina.20_10_11.bam";
    public final static String NA12878_WEx = privateTestDir + "CEUTrio.HiSeq.WEx.b37_decoy.NA12878.20_10_11mb.bam";

    public static final boolean queueTestRunModeIsSet = System.getProperty("queuetest.run", "").equals("true");

    /** before the class starts up */
    static {
        // setup a basic log configuration
        CommandLineUtils.configureConsoleLogging();

        // setup our log layout
        PatternLayout layout = new PatternLayout();
        layout.setConversionPattern("TEST %C{1}.%M - %d{HH:mm:ss,SSS} - %m%n");

        // now set the layout of all the loggers to our layout
        CommandLineUtils.setLayout(logger, layout);

        // Set the Root logger to only output warnings.
        logger.setLevel(Level.WARN);

        if (networkTempDirRootExists) {
            networkTempDirFile = IOUtils.tempDir("temp.", ".dir", new File(networkTempDirRoot + System.getProperty("user.name")));
            networkTempDirFile.deleteOnExit();
        } else {
            networkTempDirFile = null;
        }


        if ( REQUIRE_NETWORK_CONNECTION ) {
            // find our file sources
            if (!fileExist(hg18Reference) || !fileExist(hg19Reference) || !fileExist(b36KGReference)) {
                logger.fatal("We can't locate the reference directories.  Aborting!");
                throw new RuntimeException("BaseTest setup failed: unable to locate the reference directories");
            }
        }
    }

    /**
     * Simple generic utility class to creating TestNG data providers:
     *
     * 1: inherit this class, as in
     *
     *      private class SummarizeDifferenceTest extends TestDataProvider {
     *         public SummarizeDifferenceTest() {
     *           super(SummarizeDifferenceTest.class);
     *         }
     *         ...
     *      }
     *
     * Provide a reference to your class to the TestDataProvider constructor.
     *
     * 2: Create instances of your subclass.  Return from it the call to getTests, providing
     * the class type of your test
     *
     * <code>
     * {@literal @}DataProvider(name = "summaries")
     * public Object[][] createSummaries() {
     *   new SummarizeDifferenceTest().addDiff("A", "A").addSummary("A:2");
     *   new SummarizeDifferenceTest().addDiff("A", "B").addSummary("A:1", "B:1");
     *   return SummarizeDifferenceTest.getTests(SummarizeDifferenceTest.class);
     * }
     * </code>
     *
     * This class magically tracks created objects of this
     */
    public static class TestDataProvider {
        private static final Map<Class, List<Object>> tests = new HashMap<>();
        protected String name;

        /**
         * Create a new TestDataProvider instance bound to the class variable C
         */
        public TestDataProvider(Class c, String name) {
            if ( ! tests.containsKey(c) )
                tests.put(c, new ArrayList<>());
            tests.get(c).add(this);
            this.name = name;
        }

        public TestDataProvider(Class c) {
            this(c, "");
        }

        public void setName(final String name) {
            this.name = name;
        }

        /**
         * Return all of the data providers in the form expected by TestNG of type class C
         * @param c
         * @return
         */
        public static Object[][] getTests(Class c) {
            List<Object[]> params2 = new ArrayList<Object[]>();
            for ( Object x : tests.get(c) ) params2.add(new Object[]{x});
            return params2.toArray(new Object[][]{});
        }

        @Override
        public String toString() {
            return "TestDataProvider("+name+")";
        }
    }

    /**
     * test if the file exists
     *
     * @param file name as a string
     * @return true if it exists
     */
    public static boolean fileExist(String file) {
        File temp = new File(file);
        return temp.exists();
    }
    
    /**
     * this appender looks for a specific message in the log4j stream.
     * It can be used to verify that a specific message was generated to the logging system.
     */
    public static class ValidationAppender extends AppenderSkeleton {

        private boolean foundString = false;
        private String targetString = "";

        public ValidationAppender(String target) {
            targetString = target;
        }

        @Override
        protected void append(LoggingEvent loggingEvent) {
            if (loggingEvent.getMessage().equals(targetString))
                foundString = true;
        }

        public void close() {
            // do nothing
        }

        public boolean requiresLayout() {
            return false;
        }

        public boolean foundString() {
            return foundString;
        }
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(final String name, final String extension) {
        try {
            final File file = File.createTempFile(name, extension);
            file.deleteOnExit();

            // Mark corresponding indices for deletion on exit as well just in case an index is created for the temp file:
            new File(file.getAbsolutePath() + Tribble.STANDARD_INDEX_EXTENSION).deleteOnExit();
            new File(file.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION).deleteOnExit();
            new File(file.getAbsolutePath() + ".bai").deleteOnExit();
            new File(file.getAbsolutePath().replaceAll(extension + "$", ".bai")).deleteOnExit();

            return file;
        } catch (IOException ex) {
            throw new ReviewedGATKException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

    /**
     * Creates a temp list file that will be deleted on exit after tests are complete.
     * @param tempFilePrefix Prefix of the file.
     * @param lines lines to write to the file.
     * @return A list file in the temporary directory starting with tempFilePrefix, which will be deleted after the program exits.
     */
    public static File createTempListFile(final String tempFilePrefix, final String... lines) {
        try {
            final File tempListFile = createTempFile(tempFilePrefix, ".list");

            final PrintWriter out = new PrintWriter(tempListFile);
            for (final String line : lines) {
                out.println(line);
            }
            out.close();

            return tempListFile;
        } catch (IOException ex) {
            throw new ReviewedGATKException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     * @param name Name of the file.
     * @return A file in the network temporary directory with name, which will be deleted after the program exits.
     * @throws SkipException when the network is not available.
     */
    public static File tryCreateNetworkTempFile(String name) {
        if (!networkTempDirRootExists)
            throw new SkipException("Network temporary directory does not exist: " + networkTempDirRoot);
        File file = new File(networkTempDirFile, name);
        file.deleteOnExit();
        return file;
    }

    /**
     * Log this message so that it shows up inline during output as well as in html reports
     *
     * @param message
     */
    public static void log(final String message) {
        Reporter.log(message, true);
    }

    private static final double DEFAULT_FLOAT_TOLERANCE = 1e-1;

    public static final void assertEqualsDoubleSmart(final Object actual, final Double expected) {
        Assert.assertTrue(actual instanceof Double, "Not a double");
        assertEqualsDoubleSmart((double)(Double)actual, (double)expected);
    }

    public static final void assertEqualsDoubleSmart(final Object actual, final Double expected, final double tolerance) {
        Assert.assertTrue(actual instanceof Double, "Not a double");
        assertEqualsDoubleSmart((double)(Double)actual, (double)expected, tolerance);
    }

    public static final void assertEqualsDoubleSmart(final double actual, final double expected) {
        assertEqualsDoubleSmart(actual, expected, DEFAULT_FLOAT_TOLERANCE);
    }

    public static final <T> void assertEqualsSet(final Set<T> actual, final Set<T> expected, final String info) {
        final Set<T> actualSet = new HashSet<T>(actual);
        final Set<T> expectedSet = new HashSet<T>(expected);
        Assert.assertTrue(actualSet.equals(expectedSet), info); // note this is necessary due to testng bug for set comps
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance) {
        assertEqualsDoubleSmart(actual, expected, tolerance, null);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance, final String message) {
        if ( Double.isNaN(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isNaN(actual), "expected is nan, actual is not");
        else if ( Double.isInfinite(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isInfinite(actual), "expected is infinite, actual is not");
        else {
            final double delta = Math.abs(actual - expected);
            final double ratio = Math.abs(actual / expected - 1.0);
            Assert.assertTrue(delta < tolerance || ratio < tolerance, "expected = " + expected + " actual = " + actual
                    + " not within tolerance " + tolerance
                    + (message == null ? "" : "message: " + message));
        }
    }

    public static void assertVariantContextsAreEqual( final VariantContext actual, final VariantContext expected ) {
        Assert.assertNotNull(actual, "VariantContext expected not null");
        Assert.assertEquals(actual.getChr(), expected.getChr(), "chr");
        Assert.assertEquals(actual.getStart(), expected.getStart(), "start");
        Assert.assertEquals(actual.getEnd(), expected.getEnd(), "end");
        Assert.assertEquals(actual.getID(), expected.getID(), "id");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "alleles for " + expected + " vs " + actual);

        assertAttributesEquals(actual.getAttributes(), expected.getAttributes());
        Assert.assertEquals(actual.filtersWereApplied(), expected.filtersWereApplied(), "filtersWereApplied");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "isFiltered");
        assertEqualsSet(actual.getFilters(), expected.getFilters(), "filters");
        assertEqualsDoubleSmart(actual.getPhredScaledQual(), expected.getPhredScaledQual());

        Assert.assertEquals(actual.hasGenotypes(), expected.hasGenotypes(), "hasGenotypes");
        if ( expected.hasGenotypes() ) {
            assertEqualsSet(actual.getSampleNames(), expected.getSampleNames(), "sample names set");
            Assert.assertEquals(actual.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName(), "sample names");
            final Set<String> samples = expected.getSampleNames();
            for ( final String sample : samples ) {
                assertGenotypesAreEqual(actual.getGenotype(sample), expected.getGenotype(sample));
            }
        }
    }

    public static void assertVariantContextStreamsAreEqual(final Iterable<VariantContext> actual, final Iterable<VariantContext> expected) {
        final Iterator<VariantContext> actualIT = actual.iterator();
        final Iterator<VariantContext> expectedIT = expected.iterator();

        while ( expectedIT.hasNext() ) {
            final VariantContext expectedVC = expectedIT.next();
            if ( expectedVC == null )
                continue;

            VariantContext actualVC;
            do {
                Assert.assertTrue(actualIT.hasNext(), "Too few records found in actual");
                actualVC = actualIT.next();
            } while ( actualIT.hasNext() && actualVC == null );

            if ( actualVC == null )
                Assert.fail("Too few records in actual");

            assertVariantContextsAreEqual(actualVC, expectedVC);
        }
        Assert.assertTrue(! actualIT.hasNext(), "Too many records found in actual");
    }


    public static void assertGenotypesAreEqual(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype names");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "Genotype alleles");
        Assert.assertEquals(actual.getGenotypeString(), expected.getGenotypeString(), "Genotype string");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");

        // filters are the same
        Assert.assertEquals(actual.getFilters(), expected.getFilters(), "Genotype fields");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "Genotype isFiltered");

        // inline attributes
        Assert.assertEquals(actual.getDP(), expected.getDP(), "Genotype dp");
        Assert.assertTrue(Arrays.equals(actual.getAD(), expected.getAD()));
        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "Genotype gq");
        Assert.assertEquals(actual.hasPL(), expected.hasPL(), "Genotype hasPL");
        Assert.assertEquals(actual.hasAD(), expected.hasAD(), "Genotype hasAD");
        Assert.assertEquals(actual.hasGQ(), expected.hasGQ(), "Genotype hasGQ");
        Assert.assertEquals(actual.hasDP(), expected.hasDP(), "Genotype hasDP");

        Assert.assertEquals(actual.hasLikelihoods(), expected.hasLikelihoods(), "Genotype haslikelihoods");
        Assert.assertEquals(actual.getLikelihoodsString(), expected.getLikelihoodsString(), "Genotype getlikelihoodsString");
        Assert.assertEquals(actual.getLikelihoods(), expected.getLikelihoods(), "Genotype getLikelihoods");
        Assert.assertTrue(Arrays.equals(actual.getPL(), expected.getPL()));

        Assert.assertEquals(actual.getPhredScaledQual(), expected.getPhredScaledQual(), "Genotype phredScaledQual");
        assertAttributesEquals(actual.getExtendedAttributes(), expected.getExtendedAttributes());
        Assert.assertEquals(actual.isPhased(), expected.isPhased(), "Genotype isPhased");
        Assert.assertEquals(actual.getPloidy(), expected.getPloidy(), "Genotype getPloidy");
    }

    public static void assertVCFHeadersAreEqual(final VCFHeader actual, final VCFHeader expected) {
        Assert.assertEquals(actual.getMetaDataInSortedOrder().size(), expected.getMetaDataInSortedOrder().size(), "No VCF header lines");

        // for some reason set.equals() is returning false but all paired elements are .equals().  Perhaps compare to is busted?
        //Assert.assertEquals(actual.getMetaDataInInputOrder(), expected.getMetaDataInInputOrder());
        final List<VCFHeaderLine> actualLines = new ArrayList<VCFHeaderLine>(actual.getMetaDataInSortedOrder());
        final List<VCFHeaderLine> expectedLines = new ArrayList<VCFHeaderLine>(expected.getMetaDataInSortedOrder());
        for ( int i = 0; i < actualLines.size(); i++ ) {
            Assert.assertEquals(actualLines.get(i), expectedLines.get(i), "VCF header lines");
        }
    }

    public static void assertVCFandBCFFilesAreTheSame(final File vcfFile, final File bcfFile) throws IOException {
        final Pair<VCFHeader, VCIterable<LineIterator>> vcfData = VCIterable.readAllVCs(vcfFile, new VCFCodec());
        final Pair<VCFHeader, VCIterable<PositionalBufferedStream>> bcfData = VCIterable.readAllVCs(bcfFile, new BCF2Codec());
        assertVCFHeadersAreEqual(bcfData.getFirst(), vcfData.getFirst());
        assertVariantContextStreamsAreEqual(bcfData.getSecond(), vcfData.getSecond());
    }

    private static void assertAttributeEquals(final String key, final Object actual, final Object expected) {
        if ( expected instanceof Double ) {
            // must be very tolerant because doubles are being rounded to 2 sig figs
            assertEqualsDoubleSmart(actual, (Double) expected, 1e-2);
        } else
            Assert.assertEquals(actual, expected, "Attribute " + key);
    }

    private static void assertAttributesEquals(final Map<String, Object> actual, Map<String, Object> expected) {
        final Set<String> expectedKeys = new HashSet<String>(expected.keySet());

        for ( final Map.Entry<String, Object> act : actual.entrySet() ) {
            final Object actualValue = act.getValue();
            if ( expected.containsKey(act.getKey()) && expected.get(act.getKey()) != null ) {
                final Object expectedValue = expected.get(act.getKey());
                if ( expectedValue instanceof List ) {
                    final List<Object> expectedList = (List<Object>)expectedValue;
                    Assert.assertTrue(actualValue instanceof List, act.getKey() + " should be a list but isn't");
                    final List<Object> actualList = (List<Object>)actualValue;
                    Assert.assertEquals(actualList.size(), expectedList.size(), act.getKey() + " size");
                    for ( int i = 0; i < expectedList.size(); i++ )
                        assertAttributeEquals(act.getKey(), actualList.get(i), expectedList.get(i));
                } else
                    assertAttributeEquals(act.getKey(), actualValue, expectedValue);
            } else {
                // it's ok to have a binding in x -> null that's absent in y
                Assert.assertNull(actualValue, act.getKey() + " present in one but not in the other");
            }
            expectedKeys.remove(act.getKey());
        }

        // now expectedKeys contains only the keys found in expected but not in actual,
        // and they must all be null
        for ( final String missingExpected : expectedKeys ) {
            final Object value = expected.get(missingExpected);
            Assert.assertTrue(isMissing(value), "Attribute " + missingExpected + " missing in one but not in other" );
        }
    }

    private static final boolean isMissing(final Object value) {
        if ( value == null ) return true;
        else if ( value.equals(VCFConstants.MISSING_VALUE_v4) ) return true;
        else if ( value instanceof List ) {
            // handles the case where all elements are null or the list is empty
            for ( final Object elt : (List)value)
                if ( elt != null )
                    return false;
            return true;
        } else
            return false;
    }

    /**
     * Checks whether two double array contain the same values or not.
     * @param actual actual produced array.
     * @param expected expected array.
     * @param tolerance maximum difference between double value to be consider equivalent.
     */
    protected static void assertEqualsDoubleArray(final double[] actual, final double[] expected, final double tolerance) {
        if (expected == null)
            Assert.assertNull(actual);
        else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual.length,expected.length,"array length");
        }
        for (int i = 0; i < actual.length; i++)
            Assert.assertEquals(actual[i],expected[i],tolerance,"array position " + i);
    }
}
