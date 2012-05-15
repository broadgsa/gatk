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
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.writer.Options;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;

/**
 * Routines for generating all sorts of VCs for testing
 *
 * @author Your Name
 * @since Date created
 */
public class VariantContextTestProvider {
    final static VCFHeader header;
    final static List<VariantContextTestData> TEST_DATAs = new ArrayList<VariantContextTestData>();
    final static VariantContext ROOT;

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
        public List<VariantContext> vcs;

        public VariantContextTestData(final VariantContextBuilder builder) {
            this(Collections.singletonList(builder.make()));
        }

        public VariantContextTestData(final VariantContext vc) {
            this(Collections.singletonList(vc));
        }

        public VariantContextTestData(final List<VariantContext> vcs) {
            this.vcs = vcs;
        }

        public boolean hasGenotypes() {
            return vcs.get(0).hasGenotypes();
        }
    }

    private final static VariantContextBuilder builder() {
        return new VariantContextBuilder(ROOT);
    }

    private final static void add(VariantContextBuilder builder) {
        TEST_DATAs.add(new VariantContextTestData(builder));
    }

    static {
        Set<VCFHeaderLine> metaData = new TreeSet<VCFHeaderLine>();
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
        metaData.add(new VCFFilterHeaderLine("FILTER1"));
        metaData.add(new VCFFilterHeaderLine("FILTER2"));

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
        metaData.add(new VCFInfoHeaderLine("INT1", 1, VCFHeaderLineType.Integer, "x"));
        metaData.add(new VCFInfoHeaderLine("INT3", 3, VCFHeaderLineType.Integer, "x"));
        metaData.add(new VCFInfoHeaderLine("INT20", 20, VCFHeaderLineType.Integer, "x"));

        add(builder().attribute("FLOAT1", 1.0));
        add(builder().attribute("FLOAT1", 100.0));
        add(builder().attribute("FLOAT1", 1000.0));
        add(builder().attribute("FLOAT1", 100000.0));
        add(builder().attribute("FLOAT1", null));
        add(builder().attribute("FLOAT3", Arrays.asList(1.0, 2.0, 3.0)));
        add(builder().attribute("FLOAT3", Arrays.asList(1000.0, 2000.0, 3000.0)));
        add(builder().attribute("FLOAT3", Arrays.asList(100000.0, 200000.0, 300000.0)));
        add(builder().attribute("FLOAT3", null));
        metaData.add(new VCFInfoHeaderLine("FLOAT1", 1, VCFHeaderLineType.Float, "x"));
        metaData.add(new VCFInfoHeaderLine("FLOAT3", 3, VCFHeaderLineType.Float, "x"));

        add(builder().attribute("FLAG", true));
        add(builder().attribute("FLAG", false));
        metaData.add(new VCFInfoHeaderLine("FLAG", 1, VCFHeaderLineType.Flag, "x"));

        add(builder().attribute("STRING1", "s1"));
        add(builder().attribute("STRING1", null));
        // TODO - renable when BCF2 spec is fixed
//        add(builder().attribute("STRING3", Arrays.asList("s1", "s2", "s3")));
//        add(builder().attribute("STRING3", null));
//        add(builder().attribute("STRING20", Arrays.asList("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20")));
        metaData.add(new VCFInfoHeaderLine("STRING1", 1, VCFHeaderLineType.String, "x"));
        metaData.add(new VCFInfoHeaderLine("STRING3", 3, VCFHeaderLineType.String, "x"));
        metaData.add(new VCFInfoHeaderLine("STRING20", 20, VCFHeaderLineType.String, "x"));

        addGenotypesData(new ArrayList<VariantContextTestData>(TEST_DATAs), metaData);

        // prep the header
        metaData.add(new VCFContigHeaderLine(VCFHeader.CONTIG_KEY, Collections.singletonMap("ID", "1"), 0));

        header = new VCFHeader(metaData);
    }

    private static void addGenotypesData(final ArrayList<VariantContextTestData> sites, Set<VCFHeaderLine> metaData) {
        // TODO
        // for each sites VC, we are going to add create two root genotypes.
        // The first is the primary, and will be added to each new test
        // The second is variable.  In some tests it's absent (testing 1 genotype), in others it is duplicated
        // 1 once, 10, 100, or 1000 times to test scaling
        // Also, create a "missing" genotype (corresponding to a . sample) in the VCF for inclusion as well.

        // test GT

        // test GQ

        // test test Integer, Float, Flag, String atomic, vector, and missing types of different lengths per sample
    }


    public static VCFHeader getHeader() {
        return header;
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
        writer.writeHeader(VariantContextTestProvider.getHeader());
        final List<VariantContext> expected = data.vcs;
        for ( VariantContext vc : expected )
            writer.add(vc);
        writer.close();

        // read in the features
        FeatureCodec<VariantContext> codec = tester.makeCodec();
        PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(tmpFile));
        FeatureCodecHeader header = codec.readHeader(pbs);
        pbs.close();
        // TODO -- test header quality

        pbs = new PositionalBufferedStream(new FileInputStream(tmpFile));
        pbs.skip(header.getHeaderEnd());

        final List<VariantContext> actual = new ArrayList<VariantContext>(expected.size());
        while ( ! pbs.isDone() ) { actual.add(codec.decode(pbs)); };

        Assert.assertEquals(actual.size(), expected.size());

        for ( int i = 0; i < expected.size(); i++ )
            VariantContextTestProvider.assertEquals(actual.get(i), expected.get(i));
    }

    public static void assertEquals( final VariantContext actual, final VariantContext expected ) {
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.getChr(), expected.getChr());
        Assert.assertEquals(actual.getStart(), expected.getStart());
        Assert.assertEquals(actual.getEnd(), expected.getEnd());
        // TODO -- expand me
    }
}