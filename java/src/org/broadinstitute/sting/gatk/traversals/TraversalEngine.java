package org.broadinstitute.sting.gatk.traversals;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.util.RuntimeIOException;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;

import java.io.File;
import java.util.*;

public abstract class TraversalEngine {
    // list of reference ordered data objects
    protected List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods = null;

    // Iterator over rods
    List<Pair<String, RODIterator<? extends ReferenceOrderedDatum>>> rodIters;

    // How strict should we be with SAM/BAM parsing?
    protected ValidationStringency strictness = ValidationStringency.STRICT;

    // Time in milliseconds since we initialized this engine
    protected long startTime = -1;
    protected long lastProgressPrintTime = -1;                // When was the last time we printed our progress?

    // How long can we go without printing some progress info?
    protected long MAX_PROGRESS_PRINT_TIME = 10 * 1000;        // 10 seconds in millisecs

    // Maximum number of reads to process before finishing
    protected long maxReads = -1;

    // Name of the reads file, in BAM/SAM format
    protected List<File> readsFiles = null;                          // the name of the reads file
    protected SAMFileReader samReader = null;
    // iterator over the sam records in the readsFile
    protected Iterator<SAMRecord> samReadIter = null;

    // The verifying iterator, it does checking
    protected VerifyingSamIterator verifyingSamReadIter = null;

    // The reference data -- filename, refSeqFile, and iterator
    protected File refFileName = null;                        // the name of the reference file
    //private ReferenceSequenceFile refFile = null;
    protected FastaSequenceFile2 refFile = null;              // todo: merge FastaSequenceFile2 into picard!
    protected ReferenceIterator refIter = null;

    // Progress tracker for the sam file
    protected FileProgressTracker samReadingTracker = null;

    protected boolean DEBUGGING = false;
    protected boolean beSafeP = true;
    protected boolean SORT_ON_FLY = false;
    protected boolean DOWNSAMPLE_BY_FRACTION = false;
    protected boolean DOWNSAMPLE_BY_COVERAGE = false;
    protected boolean FILTER_UNSORTED_READS = false;
    protected boolean walkOverAllSites = false;
    protected int maxOnFlySorts = 100000;
    protected double downsamplingFraction = 1.0;
    protected int downsamplingCoverage = 0;
    protected long N_RECORDS_TO_PRINT = 100000;
    protected boolean THREADED_IO = false;
    protected int THREADED_IO_BUFFER_SIZE = 10000;
    protected boolean filterZeroMappingQualityReads = false;


    // the stored header
    protected SAMFileHeader myHeader = null;

    /**
     * our log, which we want to capture anything from this class
     */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);


    // Locations we are going to process during the traversal
    private List<GenomeLoc> locs = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Setting up the engine
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Creates a new, uninitialized TraversalEngine
     *
     * @param reads SAM/BAM file of reads
     * @param ref   Reference file in FASTA format, assumes a .dict file is also available
     * @param rods  Array of reference ordered data sets
     */
    public TraversalEngine(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        readsFiles = reads;
        refFileName = ref;
        this.rods = rods;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Manipulating the underlying engine parameters
    //
    // --------------------------------------------------------------------------------------------------------------
    //public void setRegion(final String reg) { regionStr = regionStr; }
    //public void setTraversalType(final String type) { traversalType = type; }
    public void setStrictness(final ValidationStringency s) {
        strictness = s;
    }

    public void setMaxReads(final int maxReads) {
        this.maxReads = maxReads;
    }

    public void setFilterZeroMappingQualityReads(final boolean filterZeroMappingQualityReads) {
        this.filterZeroMappingQualityReads = filterZeroMappingQualityReads;
    }

    public void setThreadedIO(final boolean threadedIO) {
        this.THREADED_IO = threadedIO;
    }

    public void setWalkOverAllSites(final boolean walkOverAllSites) {
        this.walkOverAllSites = walkOverAllSites;
    }

    private void setDebugging(final boolean d) {
        DEBUGGING = d;
    }

    private void setSafetyChecking(final boolean beSafeP) {
        if (!beSafeP)
            logger.warn("*** Turning off safety checking, I hope you know what you are doing.  Errors will result in debugging assert failures and other inscrutable messages...");
        this.beSafeP = beSafeP;
    }

    public void setFilterUnsortedReads(final boolean filterUnsorted) {
        if (!filterUnsorted)
            logger.warn("*** Turning on filtering of out of order reads, I *really* hope you know what you are doing, as you are removing data...");
        this.FILTER_UNSORTED_READS = filterUnsorted;
    }

    private void setSortOnFly(final int maxReadsToSort) {
        logger.info("Sorting read file on the fly: max reads allowed is " + maxReadsToSort);
        SORT_ON_FLY = true;
        maxOnFlySorts = maxReadsToSort;
    }

    private void setSortOnFly() { setSortOnFly(100000); }

    private void setDownsampleByFraction(final double fraction) {
        logger.info("Downsampling to approximately " + (fraction * 100.0) + "% of filtered reads");
        DOWNSAMPLE_BY_FRACTION = true;
        downsamplingFraction = fraction;
    }

    private void setDownsampleByCoverage(final int coverage) {
        logger.info("Downsampling to coverage " + coverage);
        DOWNSAMPLE_BY_COVERAGE = true;
        downsamplingCoverage = coverage;
    }

    /**
     * Sets up read filtering performed by the traversal engine.
     * @param reads Read filters to instantiate as the traversal engine walks over the data.
     * @deprecated Should only be used by old-style traversals.
     */
    @Deprecated
    public void setReadFilters( Reads reads) {
        if( reads.getDownsamplingFraction() != null ) setDownsampleByFraction(reads.getDownsamplingFraction());
        if( reads.getDownsampleToCoverage() != null ) setDownsampleByCoverage(reads.getDownsampleToCoverage());
        if( reads.getMaxOnTheFlySorts() != null ) setSortOnFly(reads.getMaxOnTheFlySorts());
        if( reads.getSafetyChecking() != null ) setSafetyChecking(reads.getSafetyChecking());
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing locations (areas of the genome we're traversing over)
    //
    // --------------------------------------------------------------------------------------------------------------


    /**
     * get the associated SAM header for our run
     * @return the header if it's stored, null if not
     */
    public SAMFileHeader getSAMHeader() {
        if ( myHeader == null )
            myHeader = samReader.getFileHeader();
        return myHeader;
    }

    /**
     * set's the SAM header for this traversal, which should
     * be the merged header in the multiple BAM file case.
     *  
     * @param myHeader the passed in header
     */

    public void setSAMHeader(SAMFileHeader myHeader) {
        this.myHeader = myHeader;
    }

    /**
     * Sets the intervals over which the traversal(s) should happen.
     * @param locs
     */
    public void setLocation(final List<GenomeLoc> locs) {
        this.locs = locs;
    }


    public boolean hasLocations() {
        return this.locs != null && !this.locs.isEmpty();
    }

    public List<GenomeLoc> getLocations() {
        return this.locs;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // printing
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * @param curTime (current runtime, in millisecs)
     * @return true if the maximum interval (in millisecs) has passed since the last printing
     */
    protected boolean maxElapsedIntervalForPrinting(final long curTime) {
        return (curTime - this.lastProgressPrintTime) > MAX_PROGRESS_PRINT_TIME;
    }

    /**
     * Forward request to printProgress
     *
     * @param type
     * @param loc
     */
    public void printProgress(final String type, GenomeLoc loc) {
        printProgress(false, type, loc);
    }

    /**
     * Utility routine that prints out process information (including timing) every N records or
     * every M seconds, for N and M set in global variables.
     *
     * @param mustPrint If true, will print out info, regardless of nRecords or time interval
     * @param type      String to print out describing our atomic traversal type ("read", "locus", etc)
     * @param loc       Current location
     */
    public void printProgress(boolean mustPrint, final String type, GenomeLoc loc) {
        final long nRecords = TraversalStatistics.nRecords;
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;
        //System.out.printf("Cur = %d, last print = %d, elapsed=%.2f, nRecords=%d, met=%b%n", curTime, lastProgressPrintTime, elapsed, nRecords, maxElapsedIntervalForPrinting(curTime));

        if (mustPrint || nRecords % N_RECORDS_TO_PRINT == 0 || maxElapsedIntervalForPrinting(curTime)) {
            this.lastProgressPrintTime = curTime;
            final double secsPer1MReads = (elapsed * 1000000.0) / nRecords;
            if (loc != null)
                logger.info(String.format("[PROGRESS] Traversed to %s, processing %,d %s in %.2f secs (%.2f secs per 1M %s)", loc, nRecords, type, elapsed, secsPer1MReads, type));
            else
                logger.info(String.format("[PROGRESS] Traversed %,d %s in %.2f secs (%.2f secs per 1M %s)", nRecords, type, elapsed, secsPer1MReads, type));

            // Currently samReadingTracker will print misleading info if we're not processing the whole file

            // If an index is enabled, file read progress is meaningless because a linear
            // traversal is not being performed.  For now, don't bother printing progress.
            // TODO: Create a sam indexed read tracker that tracks based on percentage through the query.
            if (samReadingTracker != null && this.locs == null)
                logger.info(String.format("[PROGRESS]   -> %s", samReadingTracker.progressMeter()));
        }
    }

    /**
     * A passthrough method so that subclasses can report which types of traversals they're using.
     * TODO: Make this method abstract once all traversals support it.
     * @param sum Result of the computation.
     * @param <T> Type of the computation.
     */
    public <T> void printOnTraversalDone( T sum ) {
        throw new UnsupportedOperationException( "This method is a required override for new traversal engines.  Please port your traversal engine to the new style." );
    }

    /**
     * Called after a traversal to print out information about the traversal process
     *
     * @param type String describing this type of traversal ("loci", "read")
     * @param sum  The reduce result of the traversal
     * @param <T>  ReduceType of the traversal
     */
    protected <T> void printOnTraversalDone(final String type, T sum) {
        printProgress(true, type, null);
        logger.info("Traversal reduce result is " + sum);
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;
        logger.info(String.format("Total runtime %.2f secs, %.2f min, %.2f hours%n", elapsed, elapsed / 60, elapsed / 3600));
        logger.info(String.format("Traversal skipped %d reads out of %d total (%.2f%%)",
                    TraversalStatistics.nSkippedReads,
                    TraversalStatistics.nReads,
                    (TraversalStatistics.nSkippedReads * 100.0) / TraversalStatistics.nReads));
        logger.info(String.format("  -> %d unmapped reads", TraversalStatistics.nUnmappedReads));
        logger.info(String.format("  -> %d duplicate reads", TraversalStatistics.nDuplicates));
        logger.info(String.format("  -> %d non-primary reads", TraversalStatistics.nNotPrimary));
        logger.info(String.format("  -> %d reads with bad alignments", TraversalStatistics.nBadAlignments));
        logger.info(String.format("  -> %d reads with indels", TraversalStatistics.nSkippedIndels));
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Initialization
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Initialize the traversal engine.  After this point traversals can be run over the data
     *
     * @return true on success
     */
    public boolean initialize() {
        lastProgressPrintTime = startTime = System.currentTimeMillis();
        initializeReference();
        // Initial the reference ordered data iterators
        initializeRODs();

        return true;
    }

   
    @Deprecated
    protected StingSAMIterator wrapReadsIterator( final Iterator<SAMRecord> rawIterator, final boolean enableVerification ) {
        // Reads sourceInfo is gone by this point in the traversal engine.  Stub in a null and rely on the iterator to
        // throw an exception if reads info isn't present.
        StingSAMIterator wrappedIterator = StingSAMIteratorAdapter.adapt(null,rawIterator);
        wrappedIterator = applyDecoratingIterators(enableVerification,wrappedIterator);

        if (THREADED_IO) {
            logger.info(String.format("Enabling threaded I/O with buffer of %d reads", THREADED_IO_BUFFER_SIZE));
            wrappedIterator = StingSAMIteratorAdapter.adapt(null,new ThreadedIterator<SAMRecord>(wrappedIterator, THREADED_IO_BUFFER_SIZE));
        }

        return wrappedIterator;
    }

    /**
     * Repackage instance variables and call static method.
     * TODO: This method's days are numbered.
     * @param enableVerification
     * @param wrappedIterator
     * @return
     */
    protected StingSAMIterator applyDecoratingIterators( final boolean enableVerification, final StingSAMIterator wrappedIterator ) {
        return applyDecoratingIterators(enableVerification,
                                        wrappedIterator,
                                        DOWNSAMPLE_BY_FRACTION ? downsamplingFraction : null,
                                        SORT_ON_FLY ? maxOnFlySorts : null,
                                        filterZeroMappingQualityReads,
                                        beSafeP);
    }

    /**
     * WARNING: In TraversalEngine for backward compatibility ONLY.  Reads are not used as the data source, only as parameters
     *          for validation.
     */
    public static StingSAMIterator applyDecoratingIterators(boolean enableVerification,
                                                            StingSAMIterator wrappedIterator,
                                                            Double downsamplingFraction,
                                                            Integer maxOnFlySorts,
                                                            Boolean filterZeroMappingQualityReads, 
                                                            Boolean beSafeP) {
        // NOTE: this (and other filtering) should be done before on-the-fly sorting
        //  as there is no reason to sort something that we will end of throwing away
        if (downsamplingFraction != null)
            wrappedIterator = new DownsampleIterator(wrappedIterator, downsamplingFraction);

        if (maxOnFlySorts != null)
            wrappedIterator = new SortSamIterator(wrappedIterator, maxOnFlySorts);

        if (beSafeP != null && beSafeP && enableVerification)
            wrappedIterator = new VerifyingSamIterator(wrappedIterator);

        if ( filterZeroMappingQualityReads != null && filterZeroMappingQualityReads )
            wrappedIterator = StingSAMIteratorAdapter.adapt(wrappedIterator.getSourceInfo(),
                                                            new FilteringIterator(wrappedIterator, new ZeroMappingQualityReadFilterFunc()));

        return wrappedIterator;
    }

    private static class ZeroMappingQualityReadFilterFunc implements SamRecordFilter {
        public boolean filterOut(SAMRecord rec) {
            if (rec.getMappingQuality() == 0) {
                //System.out.printf("Filtering 0 mapping quality read %s%n", rec.format());
                return true;
            } else {
                return false;
            }
        }
    }

    protected SAMFileReader initializeSAMFile(File samFile) {
        // todo: fixme, this is a hack to try out dynamic merging
        if ( samFile.toString().endsWith(".list") ) {
            return null;
            // todo: omg, this is just scary, just it's just for testing purposes.  fix with the new DataSource system
        } else {
            SAMFileReader samReader = new SAMFileReader(samFile, true);
            samReader.setValidationStringency(strictness);

            final SAMFileHeader header = samReader.getFileHeader();
            logger.debug(String.format("Sort order is: " + header.getSortOrder()));

            // Kludge filename into sam file header.
            if (samReader.getFileHeader().getReadGroups().size() < 1) {
                //logger.warn("Setting header in reader " + f.getName());
                SAMReadGroupRecord rec = new SAMReadGroupRecord(samFile.getName());
                rec.setLibrary(samFile.getName());
                rec.setSample(samFile.getName());

                samReader.getFileHeader().addReadGroup(rec);
            }

            return samReader;
        }
    }

    /**
     * Prepare the reference for stream processing
     */
    protected void initializeReference() {
        if (refFileName != null) {
            //this.refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFileName);
            this.refFile = new FastaSequenceFile2(refFileName);   // todo: replace when FastaSequenceFile2 is in picard
            this.refIter = new ReferenceIterator(this.refFile);
            if (!GenomeLoc.setupRefContigOrdering(this.refFile)) {
              // We couldn't process the reference contig ordering, fail since we need it
                Utils.scareUser(String.format("We couldn't load the contig dictionary associated with %s.  At the current time we require this dictionary file to efficiently access the FASTA file.  Please use /seq/software/picard/current/bin/CreateSequenceDictionary.jar to create a sequence dictionary for your file", refFileName));
            }
        }
    }

    /**
     * Prepare the list of reference ordered data iterators for each of the rods
     *
     * @return A list of ROD iterators for getting data from each ROD
     */
    protected void initializeRODs() {
        // set up reference ordered data
        rodIters = new ArrayList<Pair<String, RODIterator<? extends ReferenceOrderedDatum>>>();
        for (ReferenceOrderedData<? extends ReferenceOrderedDatum> data : rods) {
            rodIters.add(new Pair<String, RODIterator<? extends ReferenceOrderedDatum>>(data.getName(), data.iterator()));
        }
    }

    /**
     * An inappropriately placed testing of reading the reference
     */
    protected void testReference() {
        while (true) {
            ReferenceSequence ref = refFile.nextSequence();
            logger.debug(String.format("%s %d %d", ref.getName(), ref.length(), System.currentTimeMillis()));
            printProgress(true, "loci", new GenomeLoc("foo", 1));
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // shutdown
    //
    // --------------------------------------------------------------------------------------------------------------

    public boolean shutdown() {
        // todo: actually shutdown the resources

        return true;
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // dealing with reference ordered data
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Builds a list of the reference ordered datum at loc from each of the iterators.  This function
     * assumes you are accessing the data in order.  You can't use this function for random access.  Each
     * successive call moves you along the file, consuming all data before loc.
     *
     * @param loc      The location to get the rods at
     * @return A list of ReferenceOrderDatum at loc.  ROD without a datum at loc will be null in the list
     */
    protected RefMetaDataTracker getReferenceOrderedDataAtLocus(final GenomeLoc loc) {
        RefMetaDataTracker tracks = new RefMetaDataTracker();
        for (Pair<String, RODIterator<? extends ReferenceOrderedDatum>> pair : rodIters) {
            String name = pair.getFirst();
            tracks.bind(name, pair.getSecond().seekForward(loc));
        }

        return tracks;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // processing
    //
    // --------------------------------------------------------------------------------------------------------------

    public <M, T> T traverse(Walker<M, T> walker) {
        T sum = null;
        if ( hasLocations() && walker.isReduceByInterval() ) {
            logger.info("Doing reduce by interval processing");
            if ( ! hasLocations() )
                Utils.scareUser("ReduceByInterval requested by no interval was provided");
            List<Pair<GenomeLoc, T>> map = new ArrayList<Pair<GenomeLoc, T>>(locs.size());
            for ( GenomeLoc loc : locs ) {
                ArrayList<GenomeLoc> l = new ArrayList<GenomeLoc>();
                l.add(loc);
                T intervalSum = traverse(walker, l);
                sum = intervalSum;
                map.add(new Pair<GenomeLoc, T>(loc, intervalSum));
            }
            walker.onTraversalDone(map);
        } else {
            List<GenomeLoc> l = new ArrayList<GenomeLoc>();
            if ( hasLocations() )
                l = locs;
            sum = traverse(walker, l);
        }

        printOnTraversalDone("elements", sum);
        return sum;
    }

    public <M, T> T traverse(Walker<M, T> walker, List<GenomeLoc> locations) {
        return null;
    }

    public <M,T> T traverse( Walker<M,T> walker,
                             Shard shard,
                             ShardDataProvider dataProvider,
                             T sum ) {
        throw new UnsupportedOperationException( "This method is a required override for new traversal engines.  Please port your traversal engine to the new style." );
    }

    public void verifySortOrder(final boolean requiresSortedOrder) {
        if (beSafeP && !SORT_ON_FLY && samReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            final String msg = "SAM file is not sorted in coordinate order (according to header) Walker type with given arguments requires a sorted file for correct processing";
            if (requiresSortedOrder || strictness == SAMFileReader.ValidationStringency.STRICT)
                throw new RuntimeIOException(msg);
            else if (strictness == SAMFileReader.ValidationStringency.LENIENT)
                logger.warn(msg);
        }
    }
    @Deprecated
    public SAMFileReader getSamReader() { return this.samReader; }
}
