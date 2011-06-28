package org.broadinstitute.sting.gatk.iterators;

import static org.testng.Assert.fail;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.broadinstitute.sting.utils.GenomeLocParser;

import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import org.testng.annotations.BeforeMethod;

import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;



/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 * <p/>
 * Class BoundedReadIteratorUnitTest
 * <p/>
 * tests for the bounded read iterator.
 */
public class BoundedReadIteratorUnitTest extends BaseTest {

    /** the file list and the fasta sequence */
    private List<File> fl;
    private ReferenceSequenceFile seq;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @BeforeMethod
    public void doForEachTest() throws FileNotFoundException {
        fl = new ArrayList<File>();
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testBounding() {
        logger.warn("Executing testBounding");
        // total reads expected
        final int expected = 20;
        // bound by ten reads
        BoundedReadIterator iter = new BoundedReadIterator(new testIterator(), expected);

        int count = 0;
        for (SAMRecord rec: iter) {
            count++;
        }

        Assert.assertEquals(count, expected);
    }
}

class testIterator implements StingSAMIterator {
    SAMFileHeader header;
    testIterator() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1,1,2000);
    }

    public void close() {

    }

    public boolean hasNext() {
        return true;
    }

    public SAMRecord next() {
        return ArtificialSAMUtils.createArtificialRead(header,"blah",0,1,100);
    }

    public void remove() {
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
