/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.commandline.ArgumentException;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.walkers.PrintReadsWalker;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;

/**
 * Tests selected functionality in the GenomeAnalysisEngine class
 */
public class GenomeAnalysisEngineUnitTest extends BaseTest {

    @Test(expectedExceptions=ArgumentException.class)
    public void testDuplicateSamFileHandlingSingleDuplicate() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"), new Tags()));

        testEngine.setSAMFileIDs(samFiles);
        testEngine.checkForDuplicateSamFiles();
    }

    @Test(expectedExceptions=ArgumentException.class)
    public void testDuplicateSamFileHandlingMultipleDuplicates() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleNORG.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"),  new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleNORG.bam"), new Tags()));

        testEngine.setSAMFileIDs(samFiles);
        testEngine.checkForDuplicateSamFiles();
    }

    @Test
    public void testEmptyIntervalSetHandling() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        testEngine.setWalker(new PrintReadsWalker());
        testEngine.setIntervals(new GenomeLocSortedSet(null));

        testEngine.validateSuppliedIntervals();
    }
}
