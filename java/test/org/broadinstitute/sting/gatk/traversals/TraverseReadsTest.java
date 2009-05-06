package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.CountReadsWalker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import static org.junit.Assert.fail;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

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
 * Class TraverseReadsTest
 * <p/>
 * test traversing reads
 */
public class TraverseReadsTest extends BaseTest {

    private FastaSequenceFile2 seq;
    private File bam = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/index_test.bam"); // TCGA-06-0188.aligned.duplicates_marked.bam");
    private File refFile = new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/Homo_sapiens_assembly17.fasta");
    private List<File> bamList;
    private Walker countReadWalker;
    private File output;
    private static long readSize = 100000;
    TraverseReads traversalEngine = null;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
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

        // Connect print stream to the output stream
        ps = new PrintStream(out);
        bamList = new ArrayList<File>();
        bamList.add(bam);
        countReadWalker = new CountReadsWalker();
        try {
            Field f = Walker.class.getDeclaredField("out");
            f.setAccessible(true);
            f.set(countReadWalker, ps);

        } catch (IllegalAccessException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (NoSuchFieldException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("Couldn't set the walkers printstream");
        }
        List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods = new ArrayList<ReferenceOrderedData<? extends ReferenceOrderedDatum>>();

        traversalEngine = new TraverseReads(bamList, refFile, rods);


    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testMappedReadCount() {

        IndexedFastaSequenceFile ref = null;
        try {
            ref = new IndexedFastaSequenceFile(refFile);
        }
        catch (FileNotFoundException ex) {
            throw new RuntimeException("File not found opening fasta file; please do this check before MicroManaging", ex);
        }
        try {
            Thread.sleep(5000);
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        GenomeLoc.setupRefContigOrdering(ref);

        ShardStrategy shardStrategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS,
                ref.getSequenceDictionary(),
                readSize);
        SAMDataSource dataSource = null;
        try {
            dataSource = new SAMDataSource(TraversalEngine.unpackReads(bamList));
            dataSource.viewUnmappedReads(true);
            //dataSource.viewUnmappedReads(false);
        }
        catch (SimpleDataSourceLoadException ex) {
            throw new RuntimeException(ex);
        }
        catch (FileNotFoundException ex) {
            throw new RuntimeException(ex);
        }

        dataSource.viewUnmappedReads(false);

        boolean walkerInitialized = false;
        Object accumulator = null;
        while (shardStrategy.hasNext()) {
            Shard shard = shardStrategy.next();
            BoundedReadIterator readIter = null;
            try {
                readIter = (BoundedReadIterator) dataSource.seek(shard);
            }
            catch (SimpleDataSourceLoadException ex) {
                throw new RuntimeException(ex);
            }

            //LocusContextProvider locusProvider = new LocusContextProvider( readIter );

            // set the sam header of the traversal engine
            traversalEngine.setSAMHeader(readIter.getHeader());

            if (!walkerInitialized) {
                countReadWalker.initialize();
                accumulator = ((ReadWalker<?, ?>) countReadWalker).reduceInit();
                walkerInitialized = true;

            }
            if (shard == null) {
                fail("Shard == null");
            }


            accumulator = traversalEngine.traverse(countReadWalker, shard, readIter, accumulator);
            readIter.close();

        }

        traversalEngine.printOnTraversalDone("loci", accumulator);
        countReadWalker.onTraversalDone(accumulator);

        if (!(accumulator instanceof Integer)) {
            fail("Count read walker should return an interger.");
        }
        if (((Integer) accumulator) != 9721) {
            fail("there should be 9721 mapped reads in the index file");
        }
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testUnmappedReadCount() {

        IndexedFastaSequenceFile ref = null;
        try {
            ref = new IndexedFastaSequenceFile(refFile);
        }
        catch (FileNotFoundException ex) {
            throw new RuntimeException("File not found opening fasta file; please do this check before MicroManaging", ex);
        }
        try {
            Thread.sleep(5000);
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        GenomeLoc.setupRefContigOrdering(ref);

        ShardStrategy shardStrategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS,
                ref.getSequenceDictionary(),
                readSize);
        SAMDataSource dataSource = null;
        try {
            dataSource = new SAMDataSource(TraversalEngine.unpackReads(bamList));
            dataSource.viewUnmappedReads(true);
            //dataSource.viewUnmappedReads(false);
        }
        catch (SimpleDataSourceLoadException ex) {
            throw new RuntimeException(ex);
        }
        catch (FileNotFoundException ex) {
            throw new RuntimeException(ex);
        }

        dataSource.viewUnmappedReads(true);

        boolean walkerInitialized = false;
        Object accumulator = null;
        while (shardStrategy.hasNext()) {
            Shard shard = shardStrategy.next();
            BoundedReadIterator readIter = null;
            try {
                readIter = (BoundedReadIterator) dataSource.seek(shard);
            }
            catch (SimpleDataSourceLoadException ex) {
                throw new RuntimeException(ex);
            }

            //LocusContextProvider locusProvider = new LocusContextProvider( readIter );

            // set the sam header of the traversal engine
            traversalEngine.setSAMHeader(readIter.getHeader());

            if (!walkerInitialized) {
                countReadWalker.initialize();
                accumulator = ((ReadWalker<?, ?>) countReadWalker).reduceInit();
                walkerInitialized = true;

            }
            if (shard == null) {
                fail("Shard == null");
            }


            accumulator = traversalEngine.traverse(countReadWalker, shard, readIter, accumulator);
            readIter.close();

        }

        traversalEngine.printOnTraversalDone("loci", accumulator);
        countReadWalker.onTraversalDone(accumulator);

        if (!(accumulator instanceof Integer)) {
            fail("Count read walker should return an interger.");
        }
        if (((Integer) accumulator) != 10000) {
            fail("there should be 9721 mapped reads in the index file");
        }
    }

}
