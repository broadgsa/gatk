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


import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.List;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class MRUCachingSAMSequencingDictionaryUnitTest extends BaseTest {
    private static ReferenceSequenceFile seq;
    private static SAMSequenceDictionary dict;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        dict = seq.getSequenceDictionary();
    }

    @Test
    public void testBasic() {
        final MRUCachingSAMSequenceDictionary caching = new MRUCachingSAMSequenceDictionary(dict);

        Assert.assertEquals(caching.getDictionary(), dict, "Dictionary not the one I expected");

        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            Assert.assertFalse(caching.isCached(rec.getSequenceIndex()), "Expected index to not be cached");
            Assert.assertFalse(caching.isCached(rec.getSequenceName()), "Expected contig to not be cached");

            Assert.assertEquals(caching.getSequence(rec.getSequenceName()), rec, "Couldn't query for sequence");
            Assert.assertEquals(caching.getSequence(rec.getSequenceIndex()), rec, "Couldn't query for sequence index");
            Assert.assertEquals(caching.hasContig(rec.getSequenceName()), true, "hasContig query for sequence");
            Assert.assertEquals(caching.hasContigIndex(rec.getSequenceIndex()), true, "hasContigIndex query for sequence");
            Assert.assertEquals(caching.getSequenceIndex(rec.getSequenceName()), rec.getSequenceIndex(), "Couldn't query for sequence");

            Assert.assertEquals(caching.hasContig(rec.getSequenceName() + "asdfadsfa"), false, "hasContig query for unknown sequence");
            Assert.assertEquals(caching.hasContigIndex(dict.getSequences().size()), false, "hasContigIndex query for unknown index");

            Assert.assertTrue(caching.isCached(rec.getSequenceIndex()), "Expected index to be cached");
            Assert.assertTrue(caching.isCached(rec.getSequenceName()), "Expected contig to be cached");
        }
    }

    @Test(expectedExceptions = ReviewedGATKException.class)
    public void testBadGetSequence() {
        final MRUCachingSAMSequenceDictionary caching = new MRUCachingSAMSequenceDictionary(dict);
        caching.getSequence("notInDictionary");
    }

    @Test(expectedExceptions = ReviewedGATKException.class)
    public void testBadGetSequenceIndex() {
        final MRUCachingSAMSequenceDictionary caching = new MRUCachingSAMSequenceDictionary(dict);
        caching.getSequence(dict.getSequences().size());
    }
}