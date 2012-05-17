/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.utils.variantcontext;

import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variantcontext.writer.Options;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.testng.Assert;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

/**
 * Routines for generating all sorts of VCs for testing
 *
 * @author Your Name
 * @since Date created
 */
public class VariantContextTestProvider {
    final private static boolean ENABLE_VARARRAY_TESTS = false;
    final private static boolean ENABLE_PLOIDY_TESTS = false;
    final private static boolean ENABLE_PL_TESTS = true;
    private static VCFHeader syntheticHeader;
    final static List<VariantContextTestData> TEST_DATAs = new ArrayList<VariantContextTestData>();
    private static VariantContext ROOT;

    private final static List<File> testSourceVCFs = Arrays.asList(
            new File(BaseTest.testDir + "ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf"),
            new File(BaseTest.testDir + "dbsnp_135.b37.1000.vcf")
            );

    public abstract static class VariantContextIOTest {
        public String toString() {
            return "VariantContextIOTest:" + getExtension();
        }
        public abstract String getExtension();
        public abstract FeatureCodec<VariantContext> makeCodec();
        public abstract VariantContextWriter makeWriter(final File outputFile, final EnumSet<Options> baseOptions);

        public List<VariantContext> preprocess(final VCFHeader header, List<VariantContext> vcsBeforeIO) {
            return vcsBeforeIO;
        }

        public List<VariantContext> postprocess(final VCFHeader header, List<VariantContext> vcsAfterIO) {
            return vcsAfterIO;
        }
    }

    public static class VariantContextTestData {
        public final VCFHeader header;
        public List<VariantContext> vcs;

        public VariantContextTestData(final VCFHeader header, final VariantContextBuilder builder) {
            this(header, Collections.singletonList(builder.make()));
        }

        public VariantContextTestData(final VCFHeader header, final List<VariantContext> vcs) {
            final Set<String> samples = new HashSet<String>();
            for ( final VariantContext vc : vcs )
                if ( vc.hasGenotypes() )
                    samples.addAll(vc.getSampleNames());
            this.header = samples.isEmpty() ? header : new VCFHeader(header.getMetaData(), samples);
            this.vcs = vcs;
        }

        public boolean hasGenotypes() {
            return vcs.get(0).hasGenotypes();
        }

        public String toString() {
            StringBuilder b = new StringBuilder();
            b.append("VariantContextTestData: [");
            final VariantContext vc = vcs.get(0);
            final VariantContextBuilder builder = new VariantContextBuilder(vc);
            builder.noGenotypes();
            b.append(builder.make().toString()).append(" nGenotypes = ").append(vc.getNSamples());
            if ( vcs.size() > 1 ) b.append(" ----- with another ").append(vcs.size() - 1).append(" VariantContext records");
            b.append("]");
            return b.toString();
        }
    }

    private final static VariantContextBuilder builder() {
        return new VariantContextBuilder(ROOT);
    }

    private final static void add(VariantContextBuilder builder) {
        TEST_DATAs.add(new VariantContextTestData(syntheticHeader, builder));
    }

    public static void initializeTests() throws IOException {
        createSyntheticHeader();
        makeSyntheticTests();
        makeEmpiricalTests();
    }

    private static void makeEmpiricalTests() throws IOException {
        for ( final File file : testSourceVCFs ) {
            VCFCodec codec = new VCFCodec();
            Pair<VCFHeader, List<VariantContext>> x = readAllVCs( file, codec );
            List<VariantContext> fullyDecoded = new ArrayList<VariantContext>(x.getSecond().size());
            for ( final VariantContext raw : x.getSecond() )
                fullyDecoded.add(raw.fullyDecode(x.getFirst()));
            TEST_DATAs.add(new VariantContextTestData(x.getFirst(), x.getSecond()));
        }
    }

    private static void createSyntheticHeader() {
        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();

        metaData.add(new VCFInfoHeaderLine("STRING1", 1, VCFHeaderLineType.String, "x"));
        metaData.add(new VCFInfoHeaderLine("STRING3", 3, VCFHeaderLineType.String, "x"));
        metaData.add(new VCFInfoHeaderLine("STRING20", 20, VCFHeaderLineType.String, "x"));

        metaData.add(new VCFInfoHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        metaData.add(new VCFInfoHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "Genotype Quality"));
        metaData.add(new VCFInfoHeaderLine("PL", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"));
        // prep the header
        metaData.add(new VCFContigHeaderLine(VCFHeader.CONTIG_KEY, Collections.singletonMap("ID", "1"), 0));

        metaData.add(new VCFFilterHeaderLine("FILTER1"));
        metaData.add(new VCFFilterHeaderLine("FILTER2"));

        metaData.add(new VCFInfoHeaderLine("INT1", 1, VCFHeaderLineType.Integer, "x"));
        metaData.add(new VCFInfoHeaderLine("INT3", 3, VCFHeaderLineType.Integer, "x"));
        metaData.add(new VCFInfoHeaderLine("INT20", 20, VCFHeaderLineType.Integer, "x"));
        metaData.add(new VCFInfoHeaderLine("INT.VAR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "x"));
        metaData.add(new VCFInfoHeaderLine("FLOAT1", 1, VCFHeaderLineType.Float, "x"));
        metaData.add(new VCFInfoHeaderLine("FLOAT3", 3, VCFHeaderLineType.Float, "x"));
        metaData.add(new VCFInfoHeaderLine("FLAG", 1, VCFHeaderLineType.Flag, "x"));

        syntheticHeader = new VCFHeader(metaData);
    }


    private static void makeSyntheticTests() {
        VariantContextBuilder rootBuilder = new VariantContextBuilder();
        rootBuilder.source("test");
        rootBuilder.loc("1", 10, 10);
        rootBuilder.alleles("A", "C");
        rootBuilder.unfiltered();
        ROOT = rootBuilder.make();

        add(builder());
        add(builder().alleles("A"));
        add(builder().alleles("A", "C", "T"));
        add(builder().alleles("-", "C").referenceBaseForIndel("A"));
        add(builder().alleles("-", "CAGT").referenceBaseForIndel("A"));
        add(builder().loc("1", 10, 11).alleles("C", "-").referenceBaseForIndel("A"));
        add(builder().loc("1", 10, 13).alleles("CGT", "-").referenceBaseForIndel("A"));

        // make sure filters work
        add(builder().unfiltered());
        add(builder().passFilters());
        add(builder().filters("FILTER1"));
        add(builder().filters("FILTER1", "FILTER2"));

        add(builder().log10PError(VariantContext.NO_LOG10_PERROR));
        add(builder().log10PError(-1));
        add(builder().log10PError(-1.234e6));

        add(builder().noID());
        add(builder().id("rsID12345"));


        add(builder().attribute("INT1", 1));
        add(builder().attribute("INT1", 100));
        add(builder().attribute("INT1", 1000));
        add(builder().attribute("INT1", 100000));
        add(builder().attribute("INT1", null));
        add(builder().attribute("INT3", Arrays.asList(1, 2, 3)));
        add(builder().attribute("INT3", Arrays.asList(1000, 2000, 3000)));
        add(builder().attribute("INT3", Arrays.asList(100000, 200000, 300000)));
        add(builder().attribute("INT3", null));
        add(builder().attribute("INT20", Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)));

        add(builder().attribute("FLOAT1", 1.0));
        add(builder().attribute("FLOAT1", 100.0));
        add(builder().attribute("FLOAT1", 1000.0));
        add(builder().attribute("FLOAT1", 100000.0));
        add(builder().attribute("FLOAT1", null));
        add(builder().attribute("FLOAT3", Arrays.asList(1.0, 2.0, 3.0)));
        add(builder().attribute("FLOAT3", Arrays.asList(1000.0, 2000.0, 3000.0)));
        add(builder().attribute("FLOAT3", Arrays.asList(100000.0, 200000.0, 300000.0)));
        add(builder().attribute("FLOAT3", null));

        add(builder().attribute("FLAG", true));
        add(builder().attribute("FLAG", false));

        add(builder().attribute("STRING1", "s1"));
        add(builder().attribute("STRING1", null));
        add(builder().attribute("STRING3", Arrays.asList("s1", "s2", "s3")));
        add(builder().attribute("STRING3", null));
        add(builder().attribute("STRING20", Arrays.asList("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20")));

        addGenotypesToTestData();
    }

    private static void addGenotypesToTestData() {
        final ArrayList<VariantContext> sites = new ArrayList<VariantContext>();

        sites.add(builder().alleles("A").make());
        sites.add(builder().alleles("A", "C", "T").make());
        sites.add(builder().alleles("-", "C").referenceBaseForIndel("A").make());
        sites.add(builder().alleles("-", "CAGT").referenceBaseForIndel("A").make());

        for ( VariantContext site : sites ) {
            addGenotypes(site);
        }
    }

    private static void addGenotypeTests( final VariantContext site, Genotype ... genotypes ) {
        // for each sites VC, we are going to add create two root genotypes.
        // The first is the primary, and will be added to each new test
        // The second is variable.  In some tests it's absent (testing 1 genotype), in others it is duplicated
        // 1 once, 10, 100, or 1000 times to test scaling

        final VariantContextBuilder builder = new VariantContextBuilder(site);

        // add a single context
        builder.genotypes(genotypes[0]);
        add(builder);

        if ( genotypes.length > 1 ) {
            // add all
            add(builder.genotypes(Arrays.asList(genotypes)));

            // add all with the last replicated 10x, 100x, 1000x times
            for ( int nCopiesOfLast : Arrays.asList(10, 100, 1000) ) {
                final GenotypesContext gc = new GenotypesContext();
                final Genotype last = genotypes[genotypes.length-1];
                for ( int i = 0; i < genotypes.length - 1; i++ )
                    gc.add(genotypes[i]);
                for ( int i = 0; i < nCopiesOfLast; i++ )
                    gc.add(new Genotype("copy" + i, last));
                add(builder.genotypes(gc));
            }
        }
    }


    private static void addGenotypes( final VariantContext site) {
        final GenotypesContext gc = new GenotypesContext();

        // test ref/ref
        final Allele ref = site.getReference();
        final Allele alt1 = site.getNAlleles() > 1 ? site.getAlternateAllele(0) : null;
        final Genotype homRef = new Genotype("homRef", Arrays.asList(ref, ref));
        addGenotypeTests(site, homRef);

        if ( alt1 != null ) {
            final Genotype het = new Genotype("het", Arrays.asList(ref, alt1));
            final Genotype homVar = new Genotype("homVar", Arrays.asList(alt1, alt1));
            addGenotypeTests(site, homRef, het);
            addGenotypeTests(site, homRef, het, homVar);

            // ploidy
            if ( ENABLE_PLOIDY_TESTS ) {
                addGenotypeTests(site,
                        new Genotype("dip", Arrays.asList(ref, alt1)),
                        new Genotype("hap", Arrays.asList(ref)));

                addGenotypeTests(site,
                        new Genotype("dip", Arrays.asList(ref, alt1)),
                        new Genotype("tet", Arrays.asList(ref, alt1, alt1)));
            }
        }

        if ( ENABLE_PL_TESTS ) {
            // testing PLs
            addGenotypeTests(site,
                    new Genotype("g1", Arrays.asList(ref, ref), -1, new double[]{0, -1, -2}),
                    new Genotype("g2", Arrays.asList(ref, ref), -1, new double[]{0, -2, -3}));

            addGenotypeTests(site,
                    new Genotype("g1", Arrays.asList(ref, ref), -1, new double[]{-1, 0, -2}),
                    new Genotype("g2", Arrays.asList(ref, ref), -1, new double[]{0, -2, -3}));

            addGenotypeTests(site,
                    new Genotype("g1", Arrays.asList(ref, ref), -1, new double[]{-1, 0, -2}),
                    new Genotype("g2", Arrays.asList(ref, ref), -1, new double[]{0, -2000, -1000}));

            addGenotypeTests(site, // missing PLs
                    new Genotype("g1", Arrays.asList(ref, ref), -1, new double[]{-1, 0, -2}),
                    new Genotype("g2", Arrays.asList(ref, ref), -1));
        }

        // test attributes
        addGenotypeTests(site,
                attr("g1", ref, "INT1", 1),
                attr("g2", ref, "INT1", 2));
        addGenotypeTests(site,
                attr("g1", ref, "INT1", 1),
                attr("g2", ref, "INT1"));
        addGenotypeTests(site,
                attr("g1", ref, "INT3", 1, 2, 3),
                attr("g2", ref, "INT3", 4, 5, 6));
        addGenotypeTests(site,
                attr("g1", ref, "INT3", 1, 2, 3),
                attr("g2", ref, "INT3"));

        if (ENABLE_VARARRAY_TESTS) {
            addGenotypeTests(site,
                    attr("g1", ref, "INT.VAR", 1, 2, 3),
                    attr("g2", ref, "INT.VAR", 4, 5),
                    attr("g3", ref, "INT.VAR", 6));
            addGenotypeTests(site,
                    attr("g1", ref, "INT.VAR", 1, 2, 3),
                    attr("g2", ref, "INT.VAR"),
                    attr("g3", ref, "INT.VAR", 5));
        }

        addGenotypeTests(site,
                attr("g1", ref, "FLOAT1", 1.0),
                attr("g2", ref, "FLOAT1", 2.0));
        addGenotypeTests(site,
                attr("g1", ref, "FLOAT1", 1.0),
                attr("g2", ref, "FLOAT1"));
        addGenotypeTests(site,
                attr("g1", ref, "FLOAT3", 1.0, 2.0, 3.0),
                attr("g2", ref, "FLOAT3", 4.0, 5.0, 6.0));
        addGenotypeTests(site,
                attr("g1", ref, "FLOAT3", 1.0, 2.0, 3.0),
                attr("g2", ref, "FLOAT3"));

        // test test Integer, Float, Flag, String atomic, vector, and missing types of different lengths per sample
    }

    private static Genotype attr(final String name, final Allele ref, final String key, final Object ... value) {
        if ( value.length == 0 )
            return new Genotype(name, Arrays.asList(ref, ref), -1);
        else {
            final Object toAdd = value.length == 1 ? value[0] : Arrays.asList(value);
            Map<String, Object> attr = Collections.singletonMap(key, toAdd);
            return new Genotype(name, Arrays.asList(ref, ref), -1, null, attr, false);
        }
    }

    public static List<VariantContextTestData> generateSiteTests() {
        return TEST_DATAs;
    }

    public static void testReaderWriter(final VariantContextIOTest tester, final VariantContextTestData data) throws IOException {
        final File tmpFile = File.createTempFile("testReaderWriter", tester.getExtension());
        tmpFile.deleteOnExit();

        // todo -- test all options

        // write
        final EnumSet<Options> options = EnumSet.of(Options.INDEX_ON_THE_FLY);
        final VariantContextWriter writer = tester.makeWriter(tmpFile, options);
        writer.writeHeader(data.header);
        final List<VariantContext> expected = data.vcs;
        for ( VariantContext vc : expected )
            writer.add(vc);
        writer.close();

        final List<VariantContext> actual = readAllVCs(tmpFile, tester.makeCodec()).getSecond();
        Assert.assertEquals(actual.size(), expected.size());

        for ( int i = 0; i < expected.size(); i++ )
            VariantContextTestProvider.assertEquals(actual.get(i), expected.get(i));
    }

    private final static Pair<VCFHeader, List<VariantContext>> readAllVCs( final File source, final FeatureCodec<VariantContext> codec ) throws IOException {
        // read in the features
        PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(source));
        FeatureCodecHeader header = codec.readHeader(pbs);
        pbs.close();

        pbs = new PositionalBufferedStream(new FileInputStream(source));
        pbs.skip(header.getHeaderEnd());

        final List<VariantContext> read = new ArrayList<VariantContext>();
        while ( ! pbs.isDone() ) {
            final VariantContext vc = codec.decode(pbs);
            if ( vc != null ) read.add(vc);
        };

        return new Pair<VCFHeader, List<VariantContext>>((VCFHeader)header.getHeaderValue(), read);
    }

    public static void assertEquals( final VariantContext actual, final VariantContext expected ) {
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getChr(), expected.getChr());
        Assert.assertEquals(actual.getStart(), expected.getStart());
        Assert.assertEquals(actual.getEnd(), expected.getEnd());
        // TODO -- expand me
    }
}