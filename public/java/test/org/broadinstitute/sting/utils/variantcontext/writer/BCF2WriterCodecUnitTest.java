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

// our package
package org.broadinstitute.sting.utils.variantcontext.writer;


// the imports for unit testing.


import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.codecs.bcf2.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContextTestProvider;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;


public class BCF2WriterCodecUnitTest extends BaseTest {
    private static File tmpFile;
    private SAMSequenceDictionary dictionary;

//    private final static String START_VCF41_LINES = "##fileformat=VCFv4.1\n" +
//            "##reference=file://" + BaseTest.b37KGReference + "\n" +
//            "##contig=<ID=1,length=249250621,assembly=b37>\n" +
//            "##contig=<ID=2,length=243199373,assembly=b37>\n";
//
////            ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
////            ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
////            ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
////            ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
////            ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
////            ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
////            ##FILTER=<ID=q10,Description="Quality below 10">
////            ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
////            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
////            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
////            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
////            ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
//
//            private final static String SITES_HEADER_LINE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    @BeforeSuite
    public void before() throws IOException {
        tmpFile = File.createTempFile("BCF2WriterCodecUnitTest", ".bcf");
        tmpFile.delete();
        IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        dictionary = seq.getSequenceDictionary();
    }

    @BeforeMethod
    public void beforeMethod() throws IOException {
        tmpFile.delete(); // cleanup the test file
    }

    // --------------------------------------------------------------------------------
    //
    // Provider of VariantContexts for testing
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "SiteVCs")
    public Object[][] SiteVCsTest() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( VariantContextTestProvider.VariantContextsTest test : VariantContextTestProvider.generateSiteTests() )
            tests.add(new Object[]{test.vcs});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SiteVCs")
    public void testBCF2WriterReader(final List<VariantContext> contexts) throws IOException {
        // todo -- test all options

        // write
        final VariantContextWriter writer = VariantContextWriterFactory.create(tmpFile, dictionary);
        writer.writeHeader(VariantContextTestProvider.getHeader());
        for ( VariantContext vc : contexts )
            writer.add(vc);
        writer.close();

        // read in the features
        BCF2Codec codec = new BCF2Codec();
        PositionalBufferedStream pbs = new PositionalBufferedStream(new FileInputStream(tmpFile));
        FeatureCodecHeader header = codec.readHeader(pbs);
        pbs.close();
        pbs = new PositionalBufferedStream(new FileInputStream(tmpFile));
        pbs.skip(header.getHeaderEnd());

        Iterator<VariantContext> it = contexts.iterator();
        while ( ! pbs.isDone() ) {
            VariantContext vc  = it.next();
            VariantContext bcf = codec.decode(pbs);
            VariantContextTestProvider.assertEquals(vc, bcf);
        }
    }
}