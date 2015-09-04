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

package org.broadinstitute.gatk.utils.fasta;


// the imports for unit testing.


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Priority;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Basic unit test for CachingIndexedFastaSequenceFile
 */
public class CachingIndexedFastaSequenceFileUnitTest extends BaseTest {
    private File simpleFasta = new File(publicTestDir + "/exampleFASTA.fasta");
    private static final int STEP_SIZE = 1;
    private final static boolean DEBUG = false;

    //private static final List<Integer> QUERY_SIZES = Arrays.asList(1);
    private static final List<Integer> QUERY_SIZES = Arrays.asList(1, 10, 100);
    private static final List<Integer> CACHE_SIZES = Arrays.asList(-1, 100, 1000);

    @DataProvider(name = "fastas")
    public Object[][] createData1() {
        List<Object[]> params = new ArrayList<Object[]>();
        for ( File fasta : Arrays.asList(simpleFasta) ) {
            for ( int cacheSize : CACHE_SIZES ) {
                for ( int querySize : QUERY_SIZES ) {
                    params.add(new Object[]{fasta, cacheSize, querySize});
                }
            }
        }

        return params.toArray(new Object[][]{});
    }

    private static long getCacheSize(final long cacheSizeRequested) {
        return cacheSizeRequested == -1 ? CachingIndexedFastaSequenceFile.DEFAULT_CACHE_SIZE : cacheSizeRequested;
    }

    @Test(dataProvider = "fastas", enabled = true && ! DEBUG)
    public void testCachingIndexedFastaReaderSequential1(File fasta, int cacheSize, int querySize) throws FileNotFoundException {
        final CachingIndexedFastaSequenceFile caching = new CachingIndexedFastaSequenceFile(fasta, getCacheSize(cacheSize), true, false);

        SAMSequenceRecord contig = caching.getSequenceDictionary().getSequence(0);
        logger.warn(String.format("Checking contig %s length %d with cache size %d and query size %d",
                contig.getSequenceName(), contig.getSequenceLength(), cacheSize, querySize));
        testSequential(caching, fasta, querySize);
    }

    private void testSequential(final CachingIndexedFastaSequenceFile caching, final File fasta, final int querySize) throws FileNotFoundException {
        Assert.assertTrue(caching.isPreservingCase(), "testSequential only works for case preserving CachingIndexedFastaSequenceFile readers");

        final IndexedFastaSequenceFile uncached = new IndexedFastaSequenceFile(fasta);

        SAMSequenceRecord contig = uncached.getSequenceDictionary().getSequence(0);
        for ( int i = 0; i < contig.getSequenceLength(); i += STEP_SIZE ) {
            int start = i;
            int stop = start + querySize;
            if ( stop <= contig.getSequenceLength() ) {
                ReferenceSequence cachedVal = caching.getSubsequenceAt(contig.getSequenceName(), start, stop);
                ReferenceSequence uncachedVal = uncached.getSubsequenceAt(contig.getSequenceName(), start, stop);

                Assert.assertEquals(cachedVal.getName(), uncachedVal.getName());
                Assert.assertEquals(cachedVal.getContigIndex(), uncachedVal.getContigIndex());
                Assert.assertEquals(cachedVal.getBases(), uncachedVal.getBases());
            }
        }

        // asserts for efficiency.  We are going to make contig.length / STEP_SIZE queries
        // at each of range: start -> start + querySize against a cache with size of X.
        // we expect to hit the cache each time range falls within X.  We expect a hit
        // on the cache if range is within X.  Which should happen at least (X - query_size * 2) / STEP_SIZE
        // times.
        final int minExpectedHits = (int)Math.floor((Math.min(caching.getCacheSize(), contig.getSequenceLength()) - querySize * 2.0) / STEP_SIZE);
        caching.printEfficiency(Priority.WARN);
        Assert.assertTrue(caching.getCacheHits() >= minExpectedHits, "Expected at least " + minExpectedHits + " cache hits but only got " + caching.getCacheHits());

    }

    // Tests grabbing sequences around a middle cached value.
    @Test(dataProvider = "fastas", enabled = true && ! DEBUG)
    public void testCachingIndexedFastaReaderTwoStage(File fasta, int cacheSize, int querySize) throws FileNotFoundException {
        final IndexedFastaSequenceFile uncached = new IndexedFastaSequenceFile(fasta);
        final CachingIndexedFastaSequenceFile caching = new CachingIndexedFastaSequenceFile(fasta, getCacheSize(cacheSize), true, false);

        SAMSequenceRecord contig = uncached.getSequenceDictionary().getSequence(0);

        int middleStart = (contig.getSequenceLength() - querySize) / 2;
        int middleStop = middleStart + querySize;

        logger.warn(String.format("Checking contig %s length %d with cache size %d and query size %d with intermediate query",
                contig.getSequenceName(), contig.getSequenceLength(), cacheSize, querySize));

        for ( int i = 0; i < contig.getSequenceLength(); i += 10 ) {
            int start = i;
            int stop = start + querySize;
            if ( stop <= contig.getSequenceLength() ) {
                ReferenceSequence grabMiddle = caching.getSubsequenceAt(contig.getSequenceName(), middleStart, middleStop);
                ReferenceSequence cachedVal = caching.getSubsequenceAt(contig.getSequenceName(), start, stop);
                ReferenceSequence uncachedVal = uncached.getSubsequenceAt(contig.getSequenceName(), start, stop);

                Assert.assertEquals(cachedVal.getName(), uncachedVal.getName());
                Assert.assertEquals(cachedVal.getContigIndex(), uncachedVal.getContigIndex());
                Assert.assertEquals(cachedVal.getBases(), uncachedVal.getBases());
            }
        }
    }

    @DataProvider(name = "ParallelFastaTest")
    public Object[][] createParallelFastaTest() {
        List<Object[]> params = new ArrayList<Object[]>();

        for ( File fasta : Arrays.asList(simpleFasta) ) {
            for ( int cacheSize : CACHE_SIZES ) {
                for ( int querySize : QUERY_SIZES ) {
                    for ( int nt : Arrays.asList(1, 2, 3, 4) ) {
                        params.add(new Object[]{fasta, cacheSize, querySize, nt});
                    }
                }
            }
        }

        return params.toArray(new Object[][]{});
    }


    @Test(dataProvider = "ParallelFastaTest", enabled = true && ! DEBUG, timeOut = 60000)
    public void testCachingIndexedFastaReaderParallel(final File fasta, final int cacheSize, final int querySize, final int nt) throws FileNotFoundException, InterruptedException {
        final CachingIndexedFastaSequenceFile caching = new CachingIndexedFastaSequenceFile(fasta, getCacheSize(cacheSize), true, false);

        logger.warn(String.format("Parallel caching index fasta reader test cacheSize %d querySize %d nt %d", caching.getCacheSize(), querySize, nt));
        for ( int iterations = 0; iterations < 1; iterations++ ) {
            final ExecutorService executor = Executors.newFixedThreadPool(nt);
            final Collection<Callable<Object>> tasks = new ArrayList<Callable<Object>>(nt);
            for ( int i = 0; i < nt; i++ )
                tasks.add(new Callable<Object>() {
                    @Override
                    public Object call() throws Exception {
                        testSequential(caching, fasta, querySize);
                        return null;
                    }
                });
            executor.invokeAll(tasks);
            executor.shutdownNow();
        }
    }

    // make sure some bases are lower case and some are upper case
    @Test(enabled = true)
    public void testMixedCasesInExample() throws FileNotFoundException, InterruptedException {
        final IndexedFastaSequenceFile original = new IndexedFastaSequenceFile(new File(exampleFASTA));
        final CachingIndexedFastaSequenceFile casePreserving = new CachingIndexedFastaSequenceFile(new File(exampleFASTA), true);
        final CachingIndexedFastaSequenceFile allUpper = new CachingIndexedFastaSequenceFile(new File(exampleFASTA));

        int nMixedCase = 0;
        for ( SAMSequenceRecord contig : original.getSequenceDictionary().getSequences() ) {
            nMixedCase += testCases(original, casePreserving, allUpper, contig.getSequenceName(), -1, -1);

            final int step = 100;
            for ( int lastPos = step; lastPos < contig.getSequenceLength(); lastPos += step ) {
                testCases(original, casePreserving, allUpper, contig.getSequenceName(), lastPos - step, lastPos);
            }
        }

        Assert.assertTrue(nMixedCase > 0, "No mixed cases sequences found in file.  Unexpected test state");
    }

    private int testCases(final IndexedFastaSequenceFile original,
                          final IndexedFastaSequenceFile casePreserving,
                          final IndexedFastaSequenceFile allUpper,
                          final String contig, final int start, final int stop ) {
        final String orig = fetchBaseString(original, contig, start, stop);
        final String keptCase = fetchBaseString(casePreserving, contig, start, stop);
        final String upperCase = fetchBaseString(allUpper, contig, start, stop).toUpperCase();

        final String origToUpper = orig.toUpperCase();
        if ( ! orig.equals(origToUpper) ) {
            Assert.assertEquals(keptCase, orig, "Case preserving operation not equal to the original case for contig " + contig);
            Assert.assertEquals(upperCase, origToUpper, "All upper case reader not equal to the uppercase of original case for contig " + contig);
            return 1;
        } else {
            return 0;
        }
    }

    private String fetchBaseString(final IndexedFastaSequenceFile reader, final String contig, final int start, final int stop) {
        if ( start == -1 )
            return new String(reader.getSequence(contig).getBases());
        else
            return new String(reader.getSubsequenceAt(contig, start, stop).getBases());
    }

    @Test(enabled = true)
    public void testIupacChanges() throws FileNotFoundException, InterruptedException {
        final String testFasta = privateTestDir + "iupacFASTA.fasta";
        final CachingIndexedFastaSequenceFile iupacPreserving = new CachingIndexedFastaSequenceFile(new File(testFasta), false, true);
        final CachingIndexedFastaSequenceFile makeNs = new CachingIndexedFastaSequenceFile(new File(testFasta));

        int preservingNs = 0;
        int changingNs = 0;
        for ( SAMSequenceRecord contig : iupacPreserving.getSequenceDictionary().getSequences() ) {
            final String sPreserving = fetchBaseString(iupacPreserving, contig.getSequenceName(), 0, 15000);
            preservingNs += StringUtils.countMatches(sPreserving, "N");

            final String sChanging = fetchBaseString(makeNs, contig.getSequenceName(), 0, 15000);
            changingNs += StringUtils.countMatches(sChanging, "N");
        }

        Assert.assertEquals(changingNs, preservingNs + 4);
    }

    @Test(enabled = true, expectedExceptions = {UserException.class})
    public void testFailOnBadBase() throws FileNotFoundException, InterruptedException {
        final String testFasta = privateTestDir + "problematicFASTA.fasta";
        final CachingIndexedFastaSequenceFile fasta = new CachingIndexedFastaSequenceFile(new File(testFasta));

        for ( SAMSequenceRecord contig : fasta.getSequenceDictionary().getSequences() ) {
            fetchBaseString(fasta, contig.getSequenceName(), -1, -1);
        }
    }
}
