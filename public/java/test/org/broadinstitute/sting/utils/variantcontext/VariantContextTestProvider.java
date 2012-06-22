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

import org.apache.log4j.Logger;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Codec;
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
    final protected static Logger logger = Logger.getLogger(VariantContextTestProvider.class);

    final private static boolean ENABLE_VARARRAY_TESTS = true;
    final private static boolean ENABLE_PLOIDY_TESTS = true;
    final private static boolean ENABLE_PL_TESTS = true;
    final private static boolean ENABLE_SOURCE_VCF_TESTS = true;
    final private static boolean ENABLE_VARIABLE_LENGTH_GENOTYPE_STRING_TESTS = true;

    private static VCFHeader syntheticHeader;
    final static List<VariantContextTestData> TEST_DATAs = new ArrayList<VariantContextTestData>();
    private static VariantContext ROOT;

    private final static List<File> testSourceVCFs = Arrays.asList(
            new File(BaseTest.privateTestDir + "ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf"),
            new File(BaseTest.privateTestDir + "ex2.vcf"),
            new File(BaseTest.privateTestDir + "dbsnp_135.b37.1000.vcf")
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
            this(header, Collections.singletonList(builder.fullyDecoded(true).make()));
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
            b.append(builder.make().toString());
            if ( vc.getNSamples() < 5 ) {
                for ( final Genotype g : vc.getGenotypes() )
                    b.append(g.toString());
            } else {
                b.append(" nGenotypes = ").append(vc.getNSamples());
            }

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
        if ( ENABLE_SOURCE_VCF_TESTS ) {
            for ( final File file : testSourceVCFs ) {
                VCFCodec codec = new VCFCodec();
                Pair<VCFHeader, Iterable<VariantContext>> x = readAllVCs( file, codec );
                List<VariantContext> fullyDecoded = new ArrayList<VariantContext>();

                logger.warn("Reading records from " + file);
                for ( final VariantContext raw : x.getSecond() ) {
                    if ( raw != null )
                        fullyDecoded.add(raw.fullyDecode(x.getFirst()));
                }
                logger.warn("Done reading " + file);

                TEST_DATAs.add(new VariantContextTestData(x.getFirst(), fullyDecoded));
            }
        }
    }

    private final static void addHeaderLine(final Set<VCFHeaderLine> metaData, final String id, final int count, final VCFHeaderLineType type) {
        metaData.add(new VCFInfoHeaderLine(id, count, type, "x"));
        if ( type != VCFHeaderLineType.Flag )
            metaData.add(new VCFFormatHeaderLine(id, count, type, "x"));
    }

    private final static void addHeaderLine(final Set<VCFHeaderLine> metaData, final String id, final VCFHeaderLineCount count, final VCFHeaderLineType type) {
        metaData.add(new VCFInfoHeaderLine(id, count, type, "x"));
        if ( type != VCFHeaderLineType.Flag )
            metaData.add(new VCFFormatHeaderLine(id, count, type, "x"));
    }

    private static void createSyntheticHeader() {
        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();

        addHeaderLine(metaData, "STRING1", 1, VCFHeaderLineType.String);
        addHeaderLine(metaData, "STRING3", 3, VCFHeaderLineType.String);
        addHeaderLine(metaData, "STRING20", 20, VCFHeaderLineType.String);
        addHeaderLine(metaData, "VAR.INFO.STRING", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String);

        addHeaderLine(metaData, "GT", 1, VCFHeaderLineType.String);
        addHeaderLine(metaData, "GQ", 1, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "PL", VCFHeaderLineCount.G, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "GS", 2, VCFHeaderLineType.String);
        addHeaderLine(metaData, "GV", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String);
        addHeaderLine(metaData, "FT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String);

        // prep the header
        metaData.add(new VCFContigHeaderLine(VCFHeader.CONTIG_KEY, Collections.singletonMap("ID", "1"), 0));

        metaData.add(new VCFFilterHeaderLine("FILTER1"));
        metaData.add(new VCFFilterHeaderLine("FILTER2"));

        addHeaderLine(metaData, "INT1", 1, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "INT3", 3, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "INT20", 20, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "INT.VAR", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer);
        addHeaderLine(metaData, "FLOAT1", 1, VCFHeaderLineType.Float);
        addHeaderLine(metaData, "FLOAT3", 3, VCFHeaderLineType.Float);
        addHeaderLine(metaData, "FLAG", 0, VCFHeaderLineType.Flag);

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
        //add(builder().attribute("FLAG", false)); // NOTE -- VCF doesn't allow false flags

        add(builder().attribute("STRING1", "s1"));
        add(builder().attribute("STRING1", null));
        add(builder().attribute("STRING3", Arrays.asList("s1", "s2", "s3")));
        add(builder().attribute("STRING3", null));
        add(builder().attribute("STRING20", Arrays.asList("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20")));

        add(builder().attribute("VAR.INFO.STRING", "s1"));
        add(builder().attribute("VAR.INFO.STRING", Arrays.asList("s1", "s2")));
        add(builder().attribute("VAR.INFO.STRING", Arrays.asList("s1", "s2", "s3")));
        add(builder().attribute("VAR.INFO.STRING", null));

        addGenotypesToTestData();

        addComplexGenotypesTest();
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

            // add all with the last replicated 10x and 100x times
            for ( int nCopiesOfLast : Arrays.asList(10, 100, 1000) ) {
                final GenotypesContext gc = new GenotypesContext();
                final Genotype last = genotypes[genotypes.length-1];
                for ( int i = 0; i < genotypes.length - 1; i++ )
                    gc.add(genotypes[i]);
                for ( int i = 0; i < nCopiesOfLast; i++ )
                    gc.add(new GenotypeBuilder(last).name("copy" + i).make());
                add(builder.genotypes(gc));
            }
        }
    }


    private static void addGenotypes( final VariantContext site) {
        // test ref/ref
        final Allele ref = site.getReference();
        final Allele alt1 = site.getNAlleles() > 1 ? site.getAlternateAllele(0) : null;
        final Genotype homRef = GenotypeBuilder.create("homRef", Arrays.asList(ref, ref));
        addGenotypeTests(site, homRef);

        if ( alt1 != null ) {
            final Genotype het = GenotypeBuilder.create("het", Arrays.asList(ref, alt1));
            final Genotype homVar = GenotypeBuilder.create("homVar", Arrays.asList(alt1, alt1));
            addGenotypeTests(site, homRef, het);
            addGenotypeTests(site, homRef, het, homVar);

            // test no GT at all
            addGenotypeTests(site, new GenotypeBuilder("noGT", new ArrayList<Allele>(0)).attribute("INT1", 10).make());

            final List<Allele> noCall = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

            // ploidy
            if ( ENABLE_PLOIDY_TESTS ) {
                addGenotypeTests(site,
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("hap", Arrays.asList(ref)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("hap", Arrays.asList(ref)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall",  noCall),
                        GenotypeBuilder.create("noCall2", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("hap", Arrays.asList(ref)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("noCall", noCall),
                        GenotypeBuilder.create("noCall2", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));

                addGenotypeTests(site,
                        GenotypeBuilder.create("nocall", noCall),
                        GenotypeBuilder.create("dip", Arrays.asList(ref, alt1)),
                        GenotypeBuilder.create("tet", Arrays.asList(ref, alt1, alt1)));
            }


            //
            //
            // TESTING PHASE
            //
            //
            final Genotype gUnphased = new GenotypeBuilder("gUnphased", Arrays.asList(ref, alt1)).make();
            final Genotype gPhased = new GenotypeBuilder("gPhased", Arrays.asList(ref, alt1)).phased(true).make();
            final Genotype gPhased2 = new GenotypeBuilder("gPhased2", Arrays.asList(alt1, alt1)).phased(true).make();
            final Genotype gPhased3 = new GenotypeBuilder("gPhased3", Arrays.asList(ref, ref)).phased(true).make();
            final Genotype haploidNoPhase = new GenotypeBuilder("haploidNoPhase", Arrays.asList(ref)).make();
            addGenotypeTests(site, gUnphased, gPhased);
            addGenotypeTests(site, gUnphased, gPhased2);
            addGenotypeTests(site, gUnphased, gPhased3);
            addGenotypeTests(site, gPhased, gPhased2);
            addGenotypeTests(site, gPhased, gPhased3);
            addGenotypeTests(site, gPhased2, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased);
            addGenotypeTests(site, haploidNoPhase, gPhased2);
            addGenotypeTests(site, haploidNoPhase, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased, gPhased2);
            addGenotypeTests(site, haploidNoPhase, gPhased, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased2, gPhased3);
            addGenotypeTests(site, haploidNoPhase, gPhased, gPhased2, gPhased3);

            final Genotype gUnphasedTet = new GenotypeBuilder("gUnphasedTet", Arrays.asList(ref, alt1, ref, alt1)).make();
            final Genotype gPhasedTet = new GenotypeBuilder("gPhasedTet", Arrays.asList(ref, alt1, alt1, alt1)).phased(true).make();
            addGenotypeTests(site, gUnphasedTet, gPhasedTet);
        }

        if ( ENABLE_PL_TESTS ) {
            if ( site.getNAlleles() == 2 ) {
                // testing PLs
                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{0, -1, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2, -3}));

                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{-1, 0, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2, -3}));

                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{-1, 0, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2000, -1000}));

                addGenotypeTests(site, // missing PLs
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{-1, 0, -2}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref)));
            }
            else if ( site.getNAlleles() == 3 ) {
                // testing PLs
                addGenotypeTests(site,
                        GenotypeBuilder.create("g1", Arrays.asList(ref, ref), new double[]{0, -1, -2, -3, -4, -5}),
                        GenotypeBuilder.create("g2", Arrays.asList(ref, ref), new double[]{0, -2, -3, -4, -5, -6}));
            }
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

        if (ENABLE_VARIABLE_LENGTH_GENOTYPE_STRING_TESTS) {
            //
            //
            // TESTING MULTIPLE SIZED LISTS IN THE GENOTYPE FIELD
            //
            //
            addGenotypeTests(site,
                    attr("g1", ref, "GS", Arrays.asList("S1", "S2")),
                    attr("g2", ref, "GS", Arrays.asList("S3", "S4")));

            addGenotypeTests(site, // g1 is missing the string, and g2 is missing FLOAT1
                    attr("g1", ref, "FLOAT1", 1.0),
                    attr("g2", ref, "GS", Arrays.asList("S3", "S4")));

            // variable sized lists
            addGenotypeTests(site,
                    attr("g1", ref, "GV", "S1"),
                    attr("g2", ref, "GV", Arrays.asList("S3", "S4")));

            addGenotypeTests(site,
                    attr("g1", ref, "GV", Arrays.asList("S1", "S2")),
                    attr("g2", ref, "GV", Arrays.asList("S3", "S4", "S5")));

            addGenotypeTests(site, // missing value in varlist of string
                    attr("g1", ref, "FLOAT1", 1.0),
                    attr("g2", ref, "GV", Arrays.asList("S3", "S4", "S5")));


            //
            //
            // TESTING GENOTYPE FILTERS
            //
            //
            addGenotypeTests(site,
                    new GenotypeBuilder("g1-x", Arrays.asList(ref, ref)).filters("X").make(),
                    new GenotypeBuilder("g2-x", Arrays.asList(ref, ref)).filters("X").make());
            addGenotypeTests(site,
                    new GenotypeBuilder("g1-unft", Arrays.asList(ref, ref)).unfiltered().make(),
                    new GenotypeBuilder("g2-x", Arrays.asList(ref, ref)).filters("X").make());
            addGenotypeTests(site,
                    new GenotypeBuilder("g1-unft", Arrays.asList(ref, ref)).unfiltered().make(),
                    new GenotypeBuilder("g2-xy", Arrays.asList(ref, ref)).filters("X", "Y").make());
            addGenotypeTests(site,
                    new GenotypeBuilder("g1-unft", Arrays.asList(ref, ref)).unfiltered().make(),
                    new GenotypeBuilder("g2-x", Arrays.asList(ref, ref)).filters("X").make(),
                    new GenotypeBuilder("g3-xy", Arrays.asList(ref, ref)).filters("X", "Y").make());
        }
    }

    private static Genotype attr(final String name, final Allele ref, final String key, final Object ... value) {
        if ( value.length == 0 )
            return GenotypeBuilder.create(name, Arrays.asList(ref, ref));
        else {
            final Object toAdd = value.length == 1 ? value[0] : Arrays.asList(value);
            return new GenotypeBuilder(name, Arrays.asList(ref, ref)).attribute(key, toAdd).make();
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

        final Iterable<VariantContext> actual = readAllVCs(tmpFile, tester.makeCodec()).getSecond();
        assertEquals(actual, expected);

    }

    /**
     * Utility class to read all of the VC records from a file
     *
     * @param source
     * @param codec
     * @return
     * @throws IOException
     */
    private final static Pair<VCFHeader, Iterable<VariantContext>> readAllVCs( final File source, final FeatureCodec<VariantContext> codec ) throws IOException {
        // read in the features
        PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(source));
        FeatureCodecHeader header = codec.readHeader(pbs);
        pbs.close();

        pbs = new PositionalBufferedStream(new FileInputStream(source));
        pbs.skip(header.getHeaderEnd());

        final VCFHeader vcfHeader = (VCFHeader)header.getHeaderValue();
        return new Pair<VCFHeader, Iterable<VariantContext>>(vcfHeader, new VCIterable(pbs, codec, vcfHeader));
    }

    private static class VCIterable implements Iterable<VariantContext>, Iterator<VariantContext> {
        final PositionalBufferedStream pbs;
        final FeatureCodec<VariantContext> codec;
        final VCFHeader header;

        private VCIterable(final PositionalBufferedStream pbs, final FeatureCodec<VariantContext> codec, final VCFHeader header) {
            this.pbs = pbs;
            this.codec = codec;
            this.header = header;
        }

        @Override
        public Iterator<VariantContext> iterator() {
            return this;
        }

        @Override
        public boolean hasNext() {
            try {
                return ! pbs.isDone();
            } catch ( IOException e ) {
                throw new RuntimeException(e);
            }
        }

        @Override
        public VariantContext next() {
            try {
                final VariantContext vc = codec.decode(pbs);
                return vc == null ? null : vc.fullyDecode(header);
            } catch ( IOException e ) {
                throw new RuntimeException(e);
            }
        }

        @Override
        public void remove() {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    public static void assertVCFandBCFFilesAreTheSame(final File vcfFile, final File bcfFile) throws IOException {
        final Pair<VCFHeader, Iterable<VariantContext>> vcfData = readAllVCs(vcfFile, new VCFCodec());
        final Pair<VCFHeader, Iterable<VariantContext>> bcfData = readAllVCs(bcfFile, new BCF2Codec());
        assertEquals(bcfData.getFirst(), vcfData.getFirst());
        assertEquals(bcfData.getSecond(), vcfData.getSecond());
    }

    public static void assertEquals(final Iterable<VariantContext> actual, final Iterable<VariantContext> expected) {
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

            assertEquals(actualVC, expectedVC);
        }
        Assert.assertTrue(! actualIT.hasNext(), "Too many records found in actual");
    }

    /**
     * Assert that two variant contexts are actually equal
     * @param actual
     * @param expected
     */
    public static void assertEquals( final VariantContext actual, final VariantContext expected ) {
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getChr(), expected.getChr());
        Assert.assertEquals(actual.getStart(), expected.getStart());
        Assert.assertEquals(actual.getEnd(), expected.getEnd());
        Assert.assertEquals(actual.getID(), expected.getID());
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles());

        assertAttributesEquals(actual.getAttributes(), expected.getAttributes());
        Assert.assertEquals(actual.getFilters(), expected.getFilters());
        BaseTest.assertEqualsDoubleSmart(actual.getPhredScaledQual(), expected.getPhredScaledQual());

        Assert.assertEquals(actual.hasGenotypes(), expected.hasGenotypes());
        if ( expected.hasGenotypes() ) {
            Assert.assertEquals(actual.getSampleNames(), expected.getSampleNames());
            final Set<String> samples = expected.getSampleNames();
            for ( final String sample : samples ) {
                assertEquals(actual.getGenotype(sample), expected.getGenotype(sample));
            }
        }
    }

    public static void assertEquals(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName());
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles());
        Assert.assertEquals(actual.getGenotypeString(), expected.getGenotypeString());
        Assert.assertEquals(actual.getType(), expected.getType());

        // filters are the same
        Assert.assertEquals(actual.getFilters(), expected.getFilters());
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered());

        // inline attributes
        Assert.assertEquals(actual.getDP(), expected.getDP());
        Assert.assertEquals(actual.getAD(), expected.getAD());
        Assert.assertEquals(actual.getGQ(), expected.getGQ());
        Assert.assertEquals(actual.hasPL(), expected.hasPL());
        Assert.assertEquals(actual.hasAD(), expected.hasAD());
        Assert.assertEquals(actual.hasGQ(), expected.hasGQ());
        Assert.assertEquals(actual.hasDP(), expected.hasDP());

        Assert.assertEquals(actual.hasLikelihoods(), expected.hasLikelihoods());
        Assert.assertEquals(actual.getLikelihoodsString(), expected.getLikelihoodsString());
        Assert.assertEquals(actual.getLikelihoods(), expected.getLikelihoods());
        Assert.assertEquals(actual.getPL(), expected.getPL());

        Assert.assertEquals(actual.getPhredScaledQual(), expected.getPhredScaledQual());
        assertAttributesEquals(actual.getExtendedAttributes(), expected.getExtendedAttributes());
        Assert.assertEquals(actual.isPhased(), expected.isPhased());
        Assert.assertEquals(actual.getPloidy(), expected.getPloidy());
    }

    private static void assertAttributesEquals(final Map<String, Object> actual, Map<String, Object> expected) {
        final Set<String> expectedKeys = new HashSet<String>(expected.keySet());

        for ( final Map.Entry<String, Object> act : actual.entrySet() ) {
            final Object actualValue = act.getValue();
            if ( expected.containsKey(act.getKey()) && expected.get(act.getKey()) != null ) {
                final Object expectedValue = expected.get(act.getKey());
                if ( expectedValue instanceof List ) {
                    final List<Object> expectedList = (List<Object>)expectedValue;
                    Assert.assertTrue(actualValue instanceof List);
                    final List<Object> actualList = (List<Object>)actualValue;
                    Assert.assertEquals(actualList.size(), expectedList.size());
                    for ( int i = 0; i < expectedList.size(); i++ )
                        assertAttributesEquals(actualList.get(i), expectedList.get(i));
                } else
                    assertAttributesEquals(actualValue, expectedValue);
            } else {
                // it's ok to have a binding in x -> null that's absent in y
                Assert.assertNull(actualValue);
            }
            expectedKeys.remove(act.getKey());
        }

        // now expectedKeys contains only the keys found in expected but not in actual,
        // and they must all be null
        for ( final String missingExpected : expectedKeys ) {
            final Object value = expected.get(missingExpected);
            Assert.assertTrue(isMissing(value));
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

    private static void assertAttributesEquals(final Object actual, final Object expected) {
        if ( expected instanceof Double ) {
            // must be very tolerant because doubles are being rounded to 2 sig figs
            BaseTest.assertEqualsDoubleSmart(actual, (Double)expected, 1e-2);
        } else
            Assert.assertEquals(actual, expected);
    }

    public static void addComplexGenotypesTest() {
        final List<Allele> allAlleles = Arrays.asList(
                Allele.create("A", true),
                Allele.create("C", false),
                Allele.create("G", false));

        for ( int nAlleles : Arrays.asList(2, 3) ) {
            for ( int highestPloidy : Arrays.asList(1, 2, 3) ) {
                // site alleles
                final List<Allele> siteAlleles = allAlleles.subList(0, nAlleles);

                // possible alleles for genotypes
                final List<Allele> possibleGenotypeAlleles = new ArrayList<Allele>(siteAlleles);
                possibleGenotypeAlleles.add(Allele.NO_CALL);

                // there are n^ploidy possible genotypes
                final List<List<Allele>> possibleGenotypes = makeAllGenotypes(possibleGenotypeAlleles, highestPloidy);
                final int nPossibleGenotypes = possibleGenotypes.size();

                VariantContextBuilder vb = new VariantContextBuilder("unittest", "1", 1, 1, siteAlleles);

                // first test -- create n copies of each genotype
                for ( int i = 0; i < nPossibleGenotypes; i++ ) {
                    final List<Genotype> samples = new ArrayList<Genotype>(3);
                    samples.add(GenotypeBuilder.create("sample" + i, possibleGenotypes.get(i)));
                    add(vb.genotypes(samples));
                }

                // second test -- create one sample with each genotype
                {
                    final List<Genotype> samples = new ArrayList<Genotype>(nPossibleGenotypes);
                    for ( int i = 0; i < nPossibleGenotypes; i++ ) {
                        samples.add(GenotypeBuilder.create("sample" + i, possibleGenotypes.get(i)));
                    }
                    add(vb.genotypes(samples));
                }

                // test mixed ploidy
                for ( int i = 0; i < nPossibleGenotypes; i++ ) {
                    for ( int ploidy = 1; ploidy < highestPloidy; ploidy++ ) {
                        final List<Genotype> samples = new ArrayList<Genotype>(highestPloidy);
                        final List<Allele> genotype = possibleGenotypes.get(i).subList(0, ploidy);
                        samples.add(GenotypeBuilder.create("sample" + i, genotype));
                        add(vb.genotypes(samples));
                    }
                }
            }
        }
    }

    private static final List<List<Allele>> makeAllGenotypes(final List<Allele> alleles, final int highestPloidy) {
        final List<List<Allele>> combinations = new ArrayList<List<Allele>>();
        if ( highestPloidy == 1 ) {
            for ( final Allele a : alleles )
                combinations.add(Collections.singletonList(a));
        } else {
            final List<List<Allele>> sub = makeAllGenotypes(alleles, highestPloidy - 1);
            for ( List<Allele> subI : sub ) {
                for ( final Allele a : alleles ) {
                    combinations.add(Utils.cons(a, subI));
                }
            }
        }
        return combinations;
    }

    public static void assertEquals(final VCFHeader actual, final VCFHeader expected) {
        Assert.assertEquals(actual.getMetaData().size(), expected.getMetaData().size());

        // for some reason set.equals() is returning false but all paired elements are .equals().  Perhaps compare to is busted?
        //Assert.assertEquals(actual.getMetaData(), expected.getMetaData());
        final List<VCFHeaderLine> actualLines = new ArrayList<VCFHeaderLine>(actual.getMetaData());
        final List<VCFHeaderLine> expectedLines = new ArrayList<VCFHeaderLine>(expected.getMetaData());
        for ( int i = 0; i < actualLines.size(); i++ ) {
            Assert.assertEquals(actualLines.get(i), expectedLines.get(i));
        }
    }

    public static void main( String argv[] ) {
        final File variants1 = new File(argv[0]);
        final File variants2 = new File(argv[1]);
        try {
            VariantContextTestProvider.assertVCFandBCFFilesAreTheSame(variants1, variants2);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }
    }
}