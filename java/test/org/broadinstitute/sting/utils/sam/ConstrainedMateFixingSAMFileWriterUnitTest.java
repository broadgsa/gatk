package org.broadinstitute.sting.utils.sam;

import com.sun.xml.internal.messaging.saaj.packaging.mime.util.OutputUtil;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.sam.SamFileValidator;
import net.sf.samtools.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;


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
 *         <p/>
 *         Class ArtificialSAMFileWriter
 *         <p/>
 *         Test out the ArtificialSAMFileWriter class
 */
public class ConstrainedMateFixingSAMFileWriterUnitTest extends BaseTest {
    final int MAX_TEMP_FILES = 10;
    IndexedFastaSequenceFile fasta = null;
    SAMFileReader bamIn;

//    File referenceFile = new File("/Users/depristo/Desktop/broadLocal/localData/Homo_sapiens_assembly18.fasta");       // todo -- replace me with network version
//    final File BAM_FILE = new File("/Users/depristo/Desktop/broadLocal/GATK/trunk/HiSeq.test.bam");
//    final File OUTPUT_FILE = new File("/Users/depristo/Desktop/broadLocal/GATK/trunk/HiSeq.1mb.CMF.bam");

    File referenceFile = new File(hg18Reference);
    final File BAM_FILE = new File(validationDataLocation + "HiSeq.1mb.bam");
    final File OUTPUT_FILE = new File("HiSeq.1mb.CMF.bam");

    final int MAX_ISIZE_FOR_MOVE = 1000;

    @BeforeMethod
    public void beforeMethod(Object[] data) {
        if (OUTPUT_FILE.exists()) OUTPUT_FILE.delete();
        bamIn = new SAMFileReader(BAM_FILE);

        try {
            fasta = new IndexedFastaSequenceFile(referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
    }

    private ConstrainedMateFixingSAMFileWriter makeWriter(final int maxInsertSizeForMovingReadPairs) {
        return new ConstrainedMateFixingSAMFileWriter(bamIn.getFileHeader(), OUTPUT_FILE, 5, maxInsertSizeForMovingReadPairs);
    }

    private List<SAMRecord> readBAM(File file) {
        List<SAMRecord> l = new ArrayList<SAMRecord>();
        for ( SAMRecord read : new SAMFileReader(file) ) { l.add(read); }
        return l;
    }

    private void assertBamsAreEqual(File bam1File, File bam2File) {
        List<SAMRecord> reads1 = readBAM(bam1File);
        List<SAMRecord> reads2 = readBAM(bam2File);
        Assert.assertEquals(reads1.size(), reads2.size());

        for ( SAMRecord read1 : reads1 )
            Assert.assertTrue(reads2.contains(read1), "Couldn't find equal read in new BAM " + read1);
    }

    private void writeReads(Collection<SAMRecord> reads) {
        ConstrainedMateFixingSAMFileWriter writer = makeWriter(MAX_ISIZE_FOR_MOVE);
        for ( SAMRecord read : reads ) {
            writer.addAlignment(read);
        }
        writer.close();
        logger.warn("Max reads in memory: " + writer.getMaxReadsInQueue());
    }

    private void validateOutput(final File bamFile) {
        SamFileValidator validator = new SamFileValidator(new PrintWriter(System.err), MAX_TEMP_FILES);
        validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.INVALID_TAG_NM, SAMValidationError.Type.MATE_NOT_FOUND));
        boolean validated = validator.validateSamFileVerbose(new SAMFileReader(bamFile), fasta);
        Assert.assertTrue(validated, "SAM file failed to validate");
    }

    @Test(enabled = true)
    public void unmodifiedWrite() {
        writeReads(readBAM(BAM_FILE));
        validateOutput(OUTPUT_FILE);
        assertBamsAreEqual(BAM_FILE, OUTPUT_FILE);
    }

    @Test(enabled = true)
    public void writeResortingOnTheFlyNoPairs() {
        List<SAMRecord> reads = readBAM(BAM_FILE);
        for ( SAMRecord read : reads ) {
            if ( ! ConstrainedMateFixingSAMFileWriter.iSizeTooBigToMove(read, MAX_ISIZE_FOR_MOVE) )
                read.setAlignmentStart(read.getAlignmentStart() - 10);
        }
        //writeReads(Utils.sortedSAMRecords(reads));
        writeReads(reads);
        validateOutput(OUTPUT_FILE);
    }
}
