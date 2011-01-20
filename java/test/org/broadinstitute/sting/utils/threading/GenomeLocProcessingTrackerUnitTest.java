// our package
package org.broadinstitute.sting.utils.threading;


// the imports for unit testing.


import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.iterators.GenomeLocusIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.*;

/**
 * Basic unit test for GenomeLoc
 */
public class GenomeLocProcessingTrackerUnitTest extends BaseTest {
    IndexedFastaSequenceFile fasta = null;
    GenomeLocParser genomeLocParser = null;
    String chr1 = null;
    private final static String FILE_ROOT = "testdata/GLPTFile";

    @BeforeTest
    public void before() {
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

    @BeforeMethod
    public void beforeMethod(Object[] data) {
        if ( data.length > 0 )
            ((TestTarget)data[0]).init();
    }

    @AfterMethod
    public void afterMethod(Object[] data) {
        if ( data.length > 0 )
            ((TestTarget)data[0]).getTracker().close();
    }

    abstract private class TestTarget {
        String name;
        int nShards;
        int shardSize;

        public void init() {}

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

    @DataProvider(name = "threadData")
    public Object[][] createThreadData() {
        return createData(Arrays.asList(10, 100, 1000, 10000), Arrays.asList(10));
    }

    public Object[][] createData(List<Integer> nShards, List<Integer> shardSizes) {
        List<TestTarget> params = new ArrayList<TestTarget>();

        int counter = 0;
//        for ( int nShard : Arrays.asList(10,100,1000) ) {
//            for ( int shardSize : Arrays.asList(10) ) {
        for ( int nShard : nShards ) {
            for ( int shardSize : shardSizes ) {
                // shared mem -- canonical implementation
                params.add(new TestTarget("ThreadSafeSharedMemory", nShard, shardSize) {
                    SharedMemoryGenomeLocProcessingTracker tracker = GenomeLocProcessingTracker.createSharedMemory();
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                });

                final File file1 = new File(String.format("%s_ThreadSafeFileBacked_%d_%d", FILE_ROOT, counter++, nShard, shardSize));
                params.add(new TestTarget("ThreadSafeFileBacked", nShard, shardSize) {
                    GenomeLocProcessingTracker tracker = GenomeLocProcessingTracker.createFileBackedThreaded(file1, genomeLocParser);
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                    public void init() {
                        if ( file1.exists() )
                            file1.delete();
                    }
                });

                final File file2 = new File(String.format("%s_ThreadSafeFileLockingFileBacked_%d_%d", FILE_ROOT, counter++, nShard, shardSize));
                params.add(new TestTarget("ThreadSafeFileLockingFileBacked", nShard, shardSize) {
                    GenomeLocProcessingTracker tracker = GenomeLocProcessingTracker.createFileBackedDistributed(file2, genomeLocParser);
                    public GenomeLocProcessingTracker getTracker() { return tracker; }
                    public void init() {
                        if ( file2.exists() )
                            file2.delete();
                    }
                });
            }
        }

        List<Object[]> params2 = new ArrayList<Object[]>();
        for ( TestTarget x : params ) params2.add(new Object[]{x});
        return params2.toArray(new Object[][]{});
    }

    @DataProvider(name = "simpleData")
    public Object[][] createSimpleData() {
        return createData(Arrays.asList(1000), Arrays.asList(100));
    }

    private static final String NAME_ONE   = "name1";
    private static final String NAME_TWO   = "name2";

    @Test(enabled = true)
    public void testNoop() {
        GenomeLocProcessingTracker tracker = GenomeLocProcessingTracker.createNoOp();
        for ( int start = 1; start < 100; start++ ) {
            for ( int n = 0; n < 2; n++ ) {
                GenomeLoc loc = genomeLocParser.createGenomeLoc(chr1, start, start +1);
                ProcessingLoc ploc = tracker.claimOwnership(loc, NAME_ONE);
                Assert.assertTrue(ploc.isOwnedBy(NAME_ONE));
                Assert.assertEquals(tracker.getProcessingLocs().size(), 0);
            }
        }
    }

    @Test(dataProvider = "simpleData", enabled = true)
    public void testSingleProcessTracker(TestTarget test) {
        GenomeLocProcessingTracker tracker = test.getTracker();
        List<GenomeLoc> shards = test.getShards();
        logger.warn("testSingleProcessTracker " + test);

        int counter = 0;
        for ( GenomeLoc shard : shards ) {
            counter++;

            Assert.assertNull(tracker.findOwner(shard));
            Assert.assertFalse(tracker.locIsOwned(shard));

            ProcessingLoc proc = tracker.claimOwnership(shard,NAME_ONE);
            Assert.assertNotNull(proc);
            Assert.assertNotNull(proc.getLocation());
            Assert.assertNotNull(proc.getOwner());
            Assert.assertEquals(proc.getLocation(), shard);
            Assert.assertEquals(proc.getOwner(), NAME_ONE);
            Assert.assertEquals(tracker.findOwner(shard), proc);
            Assert.assertTrue(tracker.locIsOwned(shard));
            Assert.assertNotNull(tracker.getProcessingLocs());
            Assert.assertEquals(tracker.getProcessingLocs().size(), counter);

            ProcessingLoc badClaimAttempt = tracker.claimOwnership(shard,NAME_TWO);
            Assert.assertFalse(badClaimAttempt.getOwner().equals(NAME_TWO));
            Assert.assertEquals(badClaimAttempt.getOwner(), NAME_ONE);
        }
    }

    @Test(dataProvider = "simpleData", enabled = true)
    public void testIterator(TestTarget test) {
        GenomeLocProcessingTracker tracker = test.getTracker();
        List<GenomeLoc> shards = test.getShards();
        logger.warn("testIterator " + test);

        List<GenomeLoc> markedShards = new ArrayList<GenomeLoc>();
        List<GenomeLoc> toFind = new ArrayList<GenomeLoc>();

        for ( int i = 0; i < shards.size(); i++ ) {
            if ( ! (i % 10 == 0) ) {
                markedShards.add(shards.get(i));
                tracker.claimOwnership(shards.get(i), NAME_TWO);
            } else {
                toFind.add(shards.get(i));
            }
        }

        int nFound = 0;
        Iterator<GenomeLoc> it = shards.iterator();
        while ( it.hasNext() ) {
            GenomeLoc shard = tracker.claimOwnershipOfNextAvailable(it, NAME_ONE);

            if ( shard == null ) { // everything to get is done
                Assert.assertEquals(nFound, toFind.size(), "Didn't find all of the available shards");
            } else {
                nFound++;
                ProcessingLoc proc = tracker.findOwner(shard);

                Assert.assertTrue(proc.isOwnedBy(NAME_ONE));
                Assert.assertTrue(! markedShards.contains(shard), "Ran process was already marked!");
                Assert.assertTrue(toFind.contains(shard), "Claimed shard wasn't one of the unmarked!");
            }
        }
    }

    @Test(dataProvider = "simpleData", enabled = true)
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
            ProcessingLoc proc = tracker.claimOwnership(shard,NAME_ONE);

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
        boolean useIterator;

        public TestThread(TestTarget test, int count, List<GenomeLoc> toRun, boolean useIterator) {
            this.test = test;
            this.toRun = toRun;
            this.name = "thread" + count;
            this.ran = new ArrayList<GenomeLoc>();
            this.useIterator = useIterator;
        }

        public Integer call() {
            //logger.warn(String.format("Call() Thread %s", name));
            if ( useIterator ) {
                for ( GenomeLoc shard : test.getTracker().onlyOwned(toRun.iterator(), name) ) {
                    if ( shard != null ) { // ignore the unclaimable end of the stream
                        ran.add(shard);
                        // do some work here
                        for ( int sum =0, i = 0; i < 100000; i++) sum += i;
                    }
                }

            } else {
                for ( GenomeLoc shard : toRun ) {
                    //System.out.printf("Claiming ownership in %s on %s%n", name, shard);
                    ProcessingLoc proc = test.getTracker().claimOwnership(shard,name);
                    //System.out.printf("  => ownership of %s is %s (I own? %b)%n", shard, proc.getOwner(), proc.isOwnedBy(name));
                    if ( proc.isOwnedBy(name) ) {
                        ran.add(proc.getLocation());
                        // do some work here
                        for ( int sum =0, i = 0; i < 100000; i++) sum += i;
                    }
                    //logger.warn(String.format("Thread %s on %s -> owned by %s", name, shard, proc.getOwner()));
                }
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

    @Test(dataProvider = "threadData", enabled = true)
    public void testThreadedProcessesLowLevelFunctions(TestTarget test) {
        testThreading(test, false);
    }

    @Test(dataProvider = "threadData", enabled = true)
    public void testThreadedProcessesIterator(TestTarget test) {
        testThreading(test, true);
    }

    private void testThreading(TestTarget test, boolean useIterator) {
        // start up 3 threads
        logger.warn("ThreadedTesting " + test + " using iterator " + useIterator);
        List<TestThread> threads = new ArrayList<TestThread>();
        for ( int i = 0; i < 4; i++) {
            List<GenomeLoc> toRun = subList(test.getShards(), i+1);
            TestThread thread = new TestThread(test, i, toRun, useIterator);
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

                ProcessingLoc proc = tracker.findOwner(shard);
                Assert.assertNotNull(proc, "Proc was null");

                Assert.assertNotNull(proc.getOwner(), "Owner was null");
                Assert.assertEquals(proc.getLocation(), shard, "Shard loc doesn't make ProcessingLoc");

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