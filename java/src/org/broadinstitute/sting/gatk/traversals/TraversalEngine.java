package org.broadinstitute.sting.gatk.traversals;

import edu.mit.broad.picard.filter.FilteringIterator;
import edu.mit.broad.picard.filter.SamRecordFilter;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.sam.SamFileHeaderMerger;
import net.sf.functionalj.Function1;
import net.sf.functionalj.FunctionN;
import net.sf.functionalj.Functions;
import net.sf.functionalj.reflect.JdkStdReflect;
import net.sf.functionalj.reflect.StdReflect;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.RuntimeIOException;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.*;

import java.io.*;
import java.util.*;

public abstract class TraversalEngine {
    // list of reference ordered data objects
    protected List<ReferenceOrderedData> rods = null;

    // Iterator over rods
    List<ReferenceOrderedData.RODIterator> rodIters;

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
    protected File readsFile = null;                          // the name of the reads file
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

    // Number of records (loci, reads) we've processed
    protected long nRecords = 0;
    // How many reads have we processed, along with those skipped for various reasons
    protected int nReads = 0;
    protected int nSkippedReads = 0;
    protected int nUnmappedReads = 0;
    protected int nNotPrimary = 0;
    protected int nBadAlignments = 0;
    protected int nSkippedIndels = 0;

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

    /**
     * our log, which we want to capture anything from this class
     */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);


    // Locations we are going to process during the traversal
    private ArrayList<GenomeLoc> locs = null;

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
    public TraversalEngine(File reads, File ref, List<ReferenceOrderedData> rods) {
        readsFile = reads;
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

    public void setThreadedIO(final boolean threadedIO) {
        this.THREADED_IO = threadedIO;
    }

    public void setWalkOverAllSites(final boolean walkOverAllSites) {
        this.walkOverAllSites = walkOverAllSites;
    }

    public void setDebugging(final boolean d) {
        DEBUGGING = d;
    }

    public void setSafetyChecking(final boolean beSafeP) {
        if (!beSafeP)
            logger.warn("*** Turning off safety checking, I hope you know what you are doing.  Errors will result in debugging assert failures and other inscrutable messages...");
        this.beSafeP = beSafeP;
    }

    public void setFilterUnsortedReads(final boolean filterUnsorted) {
        if (!filterUnsorted)
            logger.warn("*** Turning on filtering of out of order reads, I *really* hope you know what you are doing, as you are removing data...");
        this.FILTER_UNSORTED_READS = filterUnsorted;
    }

    public void setSortOnFly(final int maxReadsToSort) {
        logger.info("Sorting read file on the fly: max reads allowed is " + maxReadsToSort);
        SORT_ON_FLY = true;
        maxOnFlySorts = maxReadsToSort;
    }

    public void setSortOnFly() { setSortOnFly(100000); }

    public void setDownsampleByFraction(final double fraction) {
        logger.info("Downsampling to approximately " + (fraction * 100.0) + "% of filtered reads");
        DOWNSAMPLE_BY_FRACTION = true;
        downsamplingFraction = fraction;
    }

    public void setDownsampleByCoverage(final int coverage) {
        logger.info("Downsampling to coverage " + coverage);
        DOWNSAMPLE_BY_COVERAGE = true;
        downsamplingCoverage = coverage;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing locations (areas of the genome we're traversing over)
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Parses the location string locStr and sets the traversal engine to only process
     * regions specified by the location string.  The string is of the form:
     * Of the form: loc1;loc2;...
     * Where each locN can be:
     * Ôchr2Õ, Ôchr2:1000000Õ or Ôchr2:1,000,000-2,000,000Õ
     *
     * @param locStr
     */
    public void setLocation(final String locStr) {
        this.locs = GenomeLoc.parseGenomeLocs(locStr);
    }

    /**
     * Read a file of genome locations to process.
     * regions specified by the location string.  The string is of the form:
     * Of the form: loc1;loc2;...
     * Where each locN can be:
     * Ôchr2Õ, Ôchr2:1000000Õ or Ôchr2:1,000,000-2,000,000Õ
     *
     * @param file_name
     */
    public void setLocationFromFile(final String file_name) {
        try {
            xReadLines reader = new xReadLines(new File(file_name));
            List<String> lines = reader.readLines();
            reader.close();
            String locStr = Utils.join(";", lines);
            logger.debug("locStr: " + locStr);
            setLocation(locStr);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }



    public boolean hasLocations() {
        return this.locs != null;
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
        final long nRecords = this.nRecords;
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
        logger.info(String.format("Traversal skipped %d reads out of %d total (%.2f%%)", nSkippedReads, nReads, (nSkippedReads * 100.0) / nReads));
        logger.info(String.format("  -> %d unmapped reads", nUnmappedReads));
        logger.info(String.format("  -> %d non-primary reads", nNotPrimary));
        logger.info(String.format("  -> %d reads with bad alignments", nBadAlignments));
        logger.info(String.format("  -> %d reads with indels", nSkippedIndels));
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

    protected Iterator<SAMRecord> initializeReads() {
        samReader = initializeSAMFile(readsFile);
        return WrapReadsIterator(getReadsIterator(samReader), true);
    }

    protected Iterator<SAMRecord> getReadsIterator(final SAMFileReader samReader) {
        // If the file has an index, querying functions are available.  Use them if possible...
        if ( samReader == null && readsFile.toString().endsWith(".list") ) {
            SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

            List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
            try {
                for ( String fileName : new xReadLines(readsFile) ) {
                    SAMFileReader reader = initializeSAMFile(new File(fileName));
                    readers.add(reader);
                }

                SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SORT_ORDER);
                return new MergingSamRecordIterator2(headerMerger);
            }
            catch ( FileNotFoundException e ) {
                logger.fatal("Couldn't open file in sam file list: " + readsFile);
            }
        }
        //if (samReader.hasIndex()) {
        //    return new SamQueryIterator(samReader, locs);
        //} else {
        return samReader.iterator();
        //}
    }

    protected Iterator<SAMRecord> WrapReadsIterator( final Iterator<SAMRecord> rawIterator, final boolean enableVerification ) {
        Iterator<SAMRecord> wrappedIterator = rawIterator;

        // NOTE: this (and other filtering) should be done before on-the-fly sorting
        //  as there is no reason to sort something that we will end of throwing away
        if (DOWNSAMPLE_BY_FRACTION)
            wrappedIterator = new DownsampleIterator(wrappedIterator, downsamplingFraction);

        if (SORT_ON_FLY)
            wrappedIterator = new SortSamIterator(wrappedIterator, maxOnFlySorts);

        if (beSafeP && enableVerification)
            wrappedIterator = new VerifyingSamIterator(wrappedIterator);

        if (THREADED_IO) {
            logger.info(String.format("Enabling threaded I/O with buffer of %d reads", THREADED_IO_BUFFER_SIZE));
            wrappedIterator = new ThreadedIterator<SAMRecord>(wrappedIterator, THREADED_IO_BUFFER_SIZE);
        }

        return wrappedIterator;
    }

    protected SAMFileReader initializeSAMFile(final File samFile) {
        // todo: fixme, this is a hack to try out dynamic merging
        if ( samFile.toString().endsWith(".list") ) {
            return null;
            // todo: omg, this is just scary, just it's just for testing purposes.  fix with the new DataSource system
        } else {
            SAMFileReader samReader = new SAMFileReader(samFile, true);
            samReader.setValidationStringency(strictness);

            final SAMFileHeader header = samReader.getFileHeader();
            logger.info(String.format("Sort order is: " + header.getSortOrder()));

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
//            if (!GenomeLoc.setupRefContigOrdering(this.refFile)) {
//                // We couldn't process the reference contig ordering, fail since we need it
//                logger.fatal(String.format("We couldn't load the contig dictionary associated with %s.  At the current time we require this dictionary file to efficiently access the FASTA file.  In the near future this program will automatically construct the dictionary for you and save it down.", refFileName));
//                throw new RuntimeException("We couldn't load the contig dictionary associated with " + refFileName + ".  At the current time we require this dictionary file to efficiently access the FASTA file.  In the near future this program will automatically construct the dictionary for you and save it down.");
//            }
        }
    }

    /**
     * Prepare the list of reference ordered data iterators for each of the rods
     *
     * @return A list of ROD iterators for getting data from each ROD
     */
    protected List<ReferenceOrderedData.RODIterator> initializeRODs() {
        // set up reference ordered data
        rodIters = new ArrayList<ReferenceOrderedData.RODIterator>();
        for (ReferenceOrderedData data : rods) {
            rodIters.add(data.iterator());
        }
        return rodIters;
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
     * @param rodIters Iterators to access the RODs
     * @param loc      The location to get the rods at
     * @return A list of ReferenceOrderDatum at loc.  ROD without a datum at loc will be null in the list
     */
    protected List<ReferenceOrderedDatum> getReferenceOrderedDataAtLocus(List<ReferenceOrderedData.RODIterator> rodIters,
                                                                         final GenomeLoc loc) {
        List<ReferenceOrderedDatum> data = new ArrayList<ReferenceOrderedDatum>();
        for (ReferenceOrderedData.RODIterator iter : rodIters) {
            data.add(iter.seekForward(loc));
        }
        return data;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // processing
    //
    // --------------------------------------------------------------------------------------------------------------

    public <M, T> T traverse(Walker<M, T> walker) {
        ArrayList<GenomeLoc> l = new ArrayList<GenomeLoc>();
        if ( hasLocations() )
            l = locs;

        return traverse(walker, l);
    }

    public <M, T> T traverse(Walker<M, T> walker, ArrayList<GenomeLoc> locations) {
        return null;
    }
    
    // --------------------------------------------------------------------------------------------------------------
    //
    // traversal by loci functions
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Class to filter out un-handle-able reads from the stream.  We currently are skipping
     * unmapped reads, non-primary reads, unaligned reads, and those with indels.  We should
     * really change this to handle indel containing reads.
     */
    class locusStreamFilterFunc implements SamRecordFilter {
        SAMRecord lastRead = null;
        public boolean filterOut(SAMRecord rec) {
            boolean result = false;
            String why = "";
            if (rec.getReadUnmappedFlag()) {
                nUnmappedReads++;
                result = true;
                why = "Unmapped";
            } else if (rec.getNotPrimaryAlignmentFlag()) {
                nNotPrimary++;
                result = true;
                why = "Not Primary";
            } else if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                nBadAlignments++;
                result = true;
                why = "No alignment start";
            } 
            else {
                result = false;
            }

            if (result) {
                nSkippedReads++;
                //System.out.printf("  [filter] %s => %b %s", rec.getReadName(), result, why);
            } else {
                nReads++;
            }
            return result;
        }
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
}
