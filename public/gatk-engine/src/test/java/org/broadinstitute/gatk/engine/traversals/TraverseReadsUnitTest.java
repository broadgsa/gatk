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

package org.broadinstitute.gatk.engine.traversals;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.gatk.engine.walkers.TestCountReadsWalker;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.reads.*;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.testng.Assert.fail;

/**
 *
 * User: aaron
 * Date: Apr 24, 2009
 * Time: 3:42:16 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 24, 2009
 * <p/>
 * Class TraverseReadsUnitTest
 * <p/>
 * test traversing reads
 */
public class TraverseReadsUnitTest extends BaseTest {

    private ReferenceSequenceFile seq;
    private SAMReaderID bam = new SAMReaderID(new File(validationDataLocation + "index_test.bam"),new Tags()); // TCGA-06-0188.aligned.duplicates_marked.bam");
    private File refFile = new File(validationDataLocation + "Homo_sapiens_assembly17.fasta");
    private List<SAMReaderID> bamList;
    private ReadWalker countReadWalker;
    private File output;
    private TraverseReadsNano traversalEngine = null;

    private IndexedFastaSequenceFile ref = null;
    private GenomeLocParser genomeLocParser = null;
    private GenomeAnalysisEngine engine = null;

    @BeforeClass
    public void doOnce() {
        try {
            ref = new CachingIndexedFastaSequenceFile(refFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(refFile,ex);
        }
        genomeLocParser = new GenomeLocParser(ref);

        engine = new GenomeAnalysisEngine();
        engine.setReferenceDataSource(refFile);
        engine.setGenomeLocParser(genomeLocParser);
    }

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @BeforeMethod
    public void doForEachTest() {
        output = new File("testOut.txt");
        FileOutputStream out = null;
        PrintStream ps; // declare a print stream object

        try {
            out = new FileOutputStream(output);
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("Couldn't open the output file");
        }

        bamList = new ArrayList<SAMReaderID>();
        bamList.add(bam);
        countReadWalker = new TestCountReadsWalker();
        
        traversalEngine = new TraverseReadsNano(1);
        traversalEngine.initialize(engine, countReadWalker);
    }

    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testUnmappedReadCount() {
        SAMDataSource dataSource = new SAMDataSource(refFile, bamList,new ThreadAllocation(),null,genomeLocParser);
        Iterable<Shard> shardStrategy = dataSource.createShardIteratorOverAllReads(new ReadShardBalancer());

        countReadWalker.initialize();
        Object accumulator = countReadWalker.reduceInit();

        for(Shard shard: shardStrategy) {
            if (shard == null) {
                fail("Shard == null");
            }

            ReadShardDataProvider dataProvider = new ReadShardDataProvider(shard,genomeLocParser,dataSource.seek(shard),null, Collections.<ReferenceOrderedDataSource>emptyList());
            accumulator = traversalEngine.traverse(countReadWalker, dataProvider, accumulator);
            dataProvider.close();
        }

        countReadWalker.onTraversalDone(accumulator);

        if (!(accumulator instanceof Long)) {
            fail("Count read walker should return a Long.");
        }
        if (!accumulator.equals(new Long(10000))) {
            fail("there should be 10000 mapped reads in the index file, there was " + (accumulator));
        }
    }

}
