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

package org.broadinstitute.gatk.utils.refdata.utils;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class FeatureToGATKFeatureIteratorUnitTest extends BaseTest {
    @Test
    @SuppressWarnings("unchecked")
    public void testCloseFilePointers() throws IOException {
        final String chr = "20";
        IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(BaseTest.hg19Reference));
        GenomeLocParser parser = new GenomeLocParser(seq);
        File file = new File(privateTestDir + "NA12878.hg19.example1.vcf");
        VCFCodec codec = new VCFCodec();
        TestFeatureReader reader = new TestFeatureReader(file.getAbsolutePath(), codec);
        CheckableCloseableTribbleIterator<Feature> tribbleIterator = reader.query(chr, 1, 100000);
        FeatureToGATKFeatureIterator gatkIterator = new FeatureToGATKFeatureIterator(parser, tribbleIterator, "test");
        Assert.assertTrue(gatkIterator.hasNext(), "GATK feature iterator does not have a next value.");
        GenomeLoc gatkLocation = gatkIterator.next().getLocation();
        Assert.assertEquals(gatkLocation.getContig(), chr, "Instead of chr 20 rod iterator was at location " + gatkLocation);
        Assert.assertFalse(tribbleIterator.isClosed(), "Tribble iterator is closed but should be still open.");
        gatkIterator.close();
        Assert.assertTrue(tribbleIterator.isClosed(), "Tribble iterator is open but should be now closed.");
        reader.close();
    }
}
