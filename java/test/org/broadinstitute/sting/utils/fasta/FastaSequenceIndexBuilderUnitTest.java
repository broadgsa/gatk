/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils.fasta;

import net.sf.picard.PicardException;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceDataSourceProgressListener;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * Test the fasta sequence index reader.
 */
public class FastaSequenceIndexBuilderUnitTest extends BaseTest {

    private FastaSequenceIndexBuilder builder;
    private ReferenceDataSourceProgressListener progress;
    private File fastaFile;
    private FastaSequenceIndex sequenceIndex;

    @Before
    public void doForEachTest() throws FileNotFoundException {
        sequenceIndex = new FastaSequenceIndex();
    }

    /**
     * Tests basic unix file with one contig
     * File is the exampleFASTA.fasta shipped with GATK
     */
    @Test
    public void unixFileTest() {
        logger.warn("Executing unixFileTest");

        fastaFile = new File(validationDataLocation + "exampleFASTA.fasta");
        builder = new FastaSequenceIndexBuilder(fastaFile, progress);
        sequenceIndex.addIndexEntry("chr1", 6, 100000, 60, 61);

        Assert.assertTrue(sequenceIndex.equals(builder.sequenceIndex));
    }


    /**
     * Tests basic windows file with one contig
     * File is a simple fasta file
     */
    @Test
    public void windowsFileTest() {
        logger.warn("Executing windowsFileTest");

        fastaFile = new File(validationDataLocation + "exampleFASTA-windows.fasta");
        builder = new FastaSequenceIndexBuilder(fastaFile, progress);
        sequenceIndex.addIndexEntry("chr2", 7, 29, 7, 9);

        Assert.assertTrue(sequenceIndex.equals(builder.sequenceIndex));
    }

    /**
     * Tests fasta with the two contigs from above combined
     * File is the exampleFASTA.fasta shipped with GATK
     */
    @Test
    public void combinedWindowsUnix() {
        logger.warn("Executing combinedWindowsUnix");

        fastaFile = new File(validationDataLocation + "exampleFASTA-combined.fasta");
        builder = new FastaSequenceIndexBuilder(fastaFile, progress);
        sequenceIndex.addIndexEntry("chr1", 6, 100000, 60, 61);
        sequenceIndex.addIndexEntry("chr2", 101680, 29, 7, 9);

        Assert.assertTrue(sequenceIndex.equals(builder.sequenceIndex));
    }

    /**
     * Tests fasta with the two contigs from above combined
     * File is the exampleFASTA.fasta shipped with GATK
     */
    @Test
    public void threeVariableLengthContigs() {
        logger.warn("Executing threeVariableLengthContigs");

        fastaFile = new File(validationDataLocation + "exampleFASTA-3contigs.fasta");
        builder = new FastaSequenceIndexBuilder(fastaFile, progress);
        sequenceIndex.addIndexEntry("chr1", 6, 17, 5, 6);
        sequenceIndex.addIndexEntry("chr2", 35, 21, 7, 8);
        sequenceIndex.addIndexEntry("chr3", 66, 100, 10, 11);

        Assert.assertTrue(sequenceIndex.equals(builder.sequenceIndex));
    }
}