// our package
package org.broadinstitute.sting.utils.fasta;


// the imports for unit testing.


import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.BeforeMethod;
import org.broadinstitute.sting.BaseTest;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Basic unit test for GenomeLoc
 */
public class CachingIndexedFastaSequenceFileUnitTest extends BaseTest {
    private File simpleFasta = new File(testDir + "/exampleFASTA.fasta");
    private static final int STEP_SIZE = 1;

    //private static final List<Integer> QUERY_SIZES = Arrays.asList(1);
    private static final List<Integer> QUERY_SIZES = Arrays.asList(1, 10, 100);
    private static final List<Integer> CACHE_SIZES = Arrays.asList(-1, 1000);

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

    @Test(dataProvider = "fastas", enabled = true)
    public void testCachingIndexedFastaReaderSequential1(File fasta, int cacheSize, int querySize) {
        IndexedFastaSequenceFile caching, uncached;
        try {
            caching = cacheSize == -1 ? new CachingIndexedFastaSequenceFile(fasta) : new CachingIndexedFastaSequenceFile(fasta, cacheSize);
            uncached = new IndexedFastaSequenceFile(fasta);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(fasta,ex);
        }

        SAMSequenceRecord contig = uncached.getSequenceDictionary().getSequence(0);
        logger.warn(String.format("Checking contig %s length %d with cache size %d and query size %d",
                contig.getSequenceName(), contig.getSequenceLength(), cacheSize, querySize));
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
    }

    // Tests grabbing sequences around a middle cached value.
    @Test(dataProvider = "fastas", enabled = true)
    public void testCachingIndexedFastaReaderTwoStage(File fasta, int cacheSize, int querySize) {
        IndexedFastaSequenceFile caching, uncached;
        try {
            uncached = new IndexedFastaSequenceFile(fasta);
            caching = new CachingIndexedFastaSequenceFile(fasta, cacheSize);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(fasta,ex);
        }

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
}
