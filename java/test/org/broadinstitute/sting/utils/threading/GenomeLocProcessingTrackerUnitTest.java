// our package
package org.broadinstitute.sting.utils.threading;


// the imports for unit testing.


import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;

/**
 * Basic unit test for GenomeLoc
 */
public class GenomeLocProcessingTrackerUnitTest extends BaseTest {
    IndexedFastaSequenceFile fasta = null;
    GenomeLocParser genomeLocParser = null;
    File sharedFile = new File("synchronizationFile.txt");
    static final boolean USE_FILE_LOCK = false;
    String chr1 = null;

    @BeforeTest
    public void before() {
        logger.warn("SharedFile is " + sharedFile.getAbsolutePath());
        File referenceFile = new File(hg18Reference);
        try {
            fasta = new IndexedFastaSequenceFile(referenceFile);
            chr1 = fasta.getSequenceDictionary().getSequence(1).getSequenceName();
            genomeLocParser = new GenomeLocParser(fasta);

        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
    }

    @AfterMethod
    public void cleanup(Object[] data) {
        if ( sharedFile.exists() ) {
            sharedFile.delete();
        }

        ((TestTarget)data[0]).getTracker().close();
    }

    abstract private class TestTarget {
        String name;
        int nShards;
        int shardSize;

        protected TestTarget(String name, int nShards, int shardSize) {
            this.name = name;
            this.nShards = nShards;
            this.shardSize = shardSize;
        }

        public abstract GenomeLocProcessingTracker getTracker();

        public List<GenomeLoc> getShards() {
            List<GenomeLoc> shards = new ArrayList<GenomeLoc>();
            for ( int i = 0; i < nShards; i++ ) {
                int start = shardSize * i;
                int stop = start + shardSize;
                shards.add(genomeLocParser.createGenomeLoc(chr1, start, stop));
            }
            return shards;
        }

        public String toString() {
            return String.format("TestTarget %s: nShards=%d shardSize=%d", name, nShards, shardSize);
        }
    }

    private class SharedFileTest extends TestTarget {
        protected SharedFileTest(int nShards, int shardSize) {
            super("SharedFile", nShards, shardSize);
        }

        GenomeLocProcessingTracker tracker = null;
        public GenomeLocProcessingTracker getTracker() {
            if ( tracker == null )
                //tracker = new SharedMemoryGenomeLocProcessingTracker();
                tracker = new SharedFileGenomeLocProcessingTracker(sharedFile, genomeLocParser, USE_FILE_LOCK);
            return tracker;
        }
    }


    @DataProvider(name = "data")
    public Object[][] createData1() {
        List<TestTarget> params = new ArrayList<TestTarget>();

        for ( int nShard : Arrays.asList(10, 100, 1000) ) {
            for ( int shardSize : Arrays.asList(10) ) {
//        for ( int nShard : Arrays.asList(10, 100, 1000, 10000) ) {
//            for ( int shardSize : Arrays.asList(10, 100) ) {
                // shared mem -- canonical implementation
                params.add(new TestTarget("SharedMem", nShard, shardSize) {
                    SharedMemoryGenomeLocProcessingTracker tracker = new SharedMemoryGenomeLocProcessingTracker();
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                });

//                // shared file -- working implementation
//                params.add(new SharedFileTest(nShard, shardSize));
            }
        }

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( TestTarget x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    private static final String NAME_ONE   = "name1";
    private static final String NAME_TWO   = "name2";

    @Test(dataProvider = "data", enabled = true)
    public void testSingleProcessTracker(TestTarget test) {
        GenomeLocProcessingTracker tracker = test.getTracker();
        List<GenomeLoc> shards = test.getShards();
        logger.warn("testSingleProcessTracker " + test);

        int counter = 0;
        for ( GenomeLoc shard : shards ) {
            counter++;

            Assert.assertNull(tracker.findOwner(shard));
            Assert.assertFalse(tracker.locIsOwned(shard));

            GenomeLocProcessingTracker.ProcessingLoc proc = tracker.claimOwnership(shard,NAME_ONE);
            Assert.assertNotNull(proc);
            Assert.assertNotNull(proc.getLoc());
            Assert.assertNotNull(proc.getOwner());
            Assert.assertEquals(proc.getLoc(), shard);
            Assert.assertEquals(proc.getOwner(), NAME_ONE);
            Assert.assertEquals(tracker.findOwner(shard), proc);
            Assert.assertTrue(tracker.locIsOwned(shard));
            Assert.assertNotNull(tracker.getProcessingLocs());
            Assert.assertEquals(tracker.getProcessingLocs().size(), counter);

            GenomeLocProcessingTracker.ProcessingLoc badClaimAttempt = tracker.claimOwnership(shard,NAME_TWO);
            Assert.assertFalse(badClaimAttempt.getOwner().equals(NAME_TWO));
            Assert.assertEquals(badClaimAttempt.getOwner(), NAME_ONE);
        }
    }

    @Test(dataProvider = "data", enabled = true)
    public void testMarkedProcesses(TestTarget test) {
        GenomeLocProcessingTracker tracker = test.getTracker();
        List<GenomeLoc> shards = test.getShards();
        logger.warn("testMarkedProcesses " + test);

        List<GenomeLoc> markedShards = new ArrayList<GenomeLoc>();

        for ( int i = 0; i < shards.size(); i++ ) {
            if ( i % 2 == 0 ) {
                markedShards.add(shards.get(i));
                tracker.claimOwnership(shards.get(i), NAME_TWO);
            }
        }

        for ( GenomeLoc shard : shards ) {
            GenomeLocProcessingTracker.ProcessingLoc proc = tracker.claimOwnership(shard,NAME_ONE);

            Assert.assertTrue(proc.isOwnedBy(NAME_ONE) || proc.isOwnedBy(NAME_TWO));

            if ( proc.isOwnedBy(NAME_ONE) )
                Assert.assertTrue(! markedShards.contains(shard), "Ran process was already marked!");
            else
                Assert.assertTrue(markedShards.contains(shard), "Unran process wasn't marked");

            if ( ! markedShards.contains(shard) ) {
                Assert.assertEquals(tracker.findOwner(shard), proc);
            }
        }
    }

    public class TestThread implements Callable<Integer> {
        public TestTarget test;
        public String name;
        public List<GenomeLoc> ran, toRun;

        public TestThread(TestTarget test, int count, List<GenomeLoc> toRun) {
            this.test = test;
            this.toRun = toRun;
            this.name = "thread" + count;
            this.ran = new ArrayList<GenomeLoc>();
        }

        public Integer call() {
            logger.warn(String.format("Call() Thread %s", name));
            for ( GenomeLoc shard : toRun ) {
                //System.out.printf("Claiming ownership in %s on %s%n", name, shard);
                GenomeLocProcessingTracker.ProcessingLoc proc = test.getTracker().claimOwnership(shard,name);
                //System.out.printf("  => ownership of %s is %s (I own? %b)%n", shard, proc.getOwner(), proc.isOwnedBy(name));
                if ( proc.isOwnedBy(name) ) {
                    ran.add(proc.getLoc());
                }
                //logger.warn(String.format("Thread %s on %s -> owned by %s", name, shard, proc.getOwner()));
            }

            return 1;
        }
    }

    private static TestThread findOwner(String name, List<TestThread> threads) {
        for ( TestThread thread : threads ) {
            if ( thread.name.equals(name) )
                return thread;
        }
        return null;
    }

    private static final <T> void assertAllThreadsFinished(List<Future<T>> futures) {
        try {
            for ( Future f : futures ) {
                Assert.assertTrue(f.isDone(), "Thread never finished running");
                Assert.assertTrue(f.get() != null, "Finished successfully");
            }
        } catch (InterruptedException e) {
            Assert.fail("Thread failed to run to completion", e);
        } catch (ExecutionException e) {
            Assert.fail("Thread generated an exception", e);
        }
    }

    private static final List<GenomeLoc> subList(List<GenomeLoc> l, int i) {
        List<GenomeLoc> r = new ArrayList<GenomeLoc>();
        for ( int j = 0; j < l.size(); j++ ) {
            if ( j % i == 0 )
                r.add(l.get(j));
        }

        return r;
    }

    @Test(dataProvider = "data", enabled = true)
    public void testThreadedProcesses(TestTarget test) {
        // start up 3 threads
        logger.warn("ThreadedTesting " + test);
        List<TestThread> threads = new ArrayList<TestThread>();
        for ( int i = 0; i < 4; i++) {
            List<GenomeLoc> toRun = subList(test.getShards(), i+1);
            TestThread thread = new TestThread(test, i, toRun);
            threads.add(thread);
        }
        ExecutorService exec = java.util.concurrent.Executors.newFixedThreadPool(threads.size());

        try {
            List<Future<Integer>> results = exec.invokeAll(threads, 300, TimeUnit.SECONDS);
            GenomeLocProcessingTracker tracker = test.getTracker();
            List<GenomeLoc> shards = test.getShards();

            for ( TestThread thread : threads )
                logger.warn(String.format("TestThread %s ran %d jobs of %d to run", thread.name, thread.ran.size(), thread.toRun.size()));

            assertAllThreadsFinished(results);

            // we ran everything
            Assert.assertEquals(tracker.getProcessingLocs().size(), shards.size(), "Not all shards were run");

            for ( GenomeLoc shard : shards ) {
                Assert.assertTrue(tracker.locIsOwned(shard), "Unowned shard");

                GenomeLocProcessingTracker.ProcessingLoc proc = tracker.findOwner(shard);
                Assert.assertNotNull(proc, "Proc was null");

                Assert.assertNotNull(proc.getOwner(), "Owner was null");
                Assert.assertEquals(proc.getLoc(), shard, "Shard loc doesn't make ProcessingLoc");

                TestThread owner = findOwner(proc.getOwner(), threads);
                Assert.assertNotNull(owner, "Couldn't find owner");

                Assert.assertTrue(owner.ran.contains(shard), "Owner doesn't contain ran shard");

                for ( TestThread thread : threads )
                    if ( ! proc.isOwnedBy(thread.name) && thread.ran.contains(shard) )
                        Assert.fail("Shard appears in another run list: proc=" + proc + " shard=" + shard + " also in jobs of " + thread.name + " obj=" + thread.ran.get(thread.ran.indexOf(shard)));

            }
        }  catch (InterruptedException e) {
            Assert.fail("Thread failure", e);
        }
    }
}