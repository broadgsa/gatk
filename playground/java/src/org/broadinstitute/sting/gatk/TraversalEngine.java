package org.broadinstitute.sting.gatk;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.RuntimeIOException;
import edu.mit.broad.picard.filter.SamRecordFilter;
import edu.mit.broad.picard.filter.FilteringIterator;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequence;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;

import java.io.*;
import java.util.*;

import net.sf.functionalj.reflect.StdReflect;
import net.sf.functionalj.reflect.JdkStdReflect;
import net.sf.functionalj.FunctionN;
import net.sf.functionalj.Function1;
import net.sf.functionalj.Functions;
import net.sf.functionalj.util.Operators;

public class TraversalEngine {
    // list of reference ordered data objects
    private List<ReferenceOrderedData> rods = null;

    // Iterator over rods
    List<ReferenceOrderedData.RODIterator> rodIters;

    //private String regionStr = null;                        // String dec
    //private String traversalType = null;                    // String describing this traversal type

    // How strict should we be with SAM/BAM parsing?
    private ValidationStringency strictness = ValidationStringency.STRICT;

    // Time in milliseconds since we initialized this engine
    private long startTime = -1;
    private long lastProgressPrintTime = -1;                // When was the last time we printed our progress?

    // How long can we go without printing some progress info?
    private long MAX_PROGRESS_PRINT_TIME = 10 * 1000;        // 10 seconds in millisecs

    // Maximum number of reads to process before finishing
    private long maxReads = -1;

    // Name of the reads file, in BAM/SAM format
    private File readsFile = null;                          // the name of the reads file
    private SAMFileReader samReader = null;
    // iterator over the sam records in the readsFile
    private Iterator<SAMRecord> samReadIter = null;

    // The verifying iterator, it does checking
    VerifyingSamIterator verifyingSamReadIter = null;

    // The reference data -- filename, refSeqFile, and iterator
    private File refFileName = null;                        // the name of the reference file
    private ReferenceSequenceFile refFile = null;
    private ReferenceIterator refIter = null;

    // Number of records (loci, reads) we've processed
    private long nRecords = 0;
    // How many reads have we processed, along with those skipped for various reasons
    private int nReads = 0;
    private int nSkippedReads = 0;
    private int nUnmappedReads = 0;
    private int nNotPrimary = 0;
    private int nBadAlignments = 0;
    private int nSkippedIndels = 0;

    // Progress tracker for the sam file
    private FileProgressTracker samReadingTracker = null;

    public boolean DEBUGGING = false;
    public boolean beSafeP = true;
    public long N_RECORDS_TO_PRINT = 100000;
    public int THREADED_IO_BUFFER_SIZE = 10000;

    // Locations we are going to process during the traversal
    private GenomeLoc[] locs = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Setting up the engine
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Creates a new, uninitialized TraversalEngine
     *
     * @param reads SAM/BAM file of reads
     * @param ref Reference file in FASTA format, assumes a .dict file is also available
     * @param rods Array of reference ordered data sets
     */
    public TraversalEngine(File reads, File ref, ReferenceOrderedData[] rods ) {
        readsFile = reads;
        refFileName = ref;
        this.rods = Arrays.asList(rods);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Manipulating the underlying engine parameters
    //
    // --------------------------------------------------------------------------------------------------------------
    //public void setRegion(final String reg) { regionStr = regionStr; }
    //public void setTraversalType(final String type) { traversalType = type; }
    public void setStrictness( final ValidationStringency s ) { strictness = s; }
    public void setMaxReads( final int maxReads ) { this.maxReads = maxReads; }
    public void setDebugging( final boolean d ) { DEBUGGING = d; }
    public void setSafetyChecking( final boolean beSafeP ) {
        if ( ! beSafeP )
            System.out.printf("*** Turning off safety checking, I hope you know what you are doing.  Errors will result in debugging assert failures and other inscrutable messages...%n");
        this.beSafeP = beSafeP;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing locations (areas of the genome we're traversing over)
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Parses the location string locStr and sets the traversal engine to only process
     * regions specified by the location string.  The string is of the form:
     *   Of the form: loc1;loc2;...
     *   Where each locN can be:
     *     Ôchr2Õ, Ôchr2:1000000Õ or Ôchr2:1,000,000-2,000,000Õ
     *
     * @param locStr
     */
    public void setLocation( final String locStr ) {
        this.locs = parseGenomeLocs(locStr);
    }

    /**
     * Useful utility function that parses a location string into a coordinate-order sorted
     * array of GenomeLoc objects
     *
     * @param str
     * @return Array of GenomeLoc objects corresponding to the locations in the string, sorted by coordinate order
     */
    public static GenomeLoc[] parseGenomeLocs( final String str ) {
        // Of the form: loc1;loc2;...
        // Where each locN can be:
        // Ôchr2Õ, Ôchr2:1000000Õ or Ôchr2:1,000,000-2,000,000Õ
        StdReflect reflect = new JdkStdReflect();
        FunctionN<GenomeLoc> parseOne = reflect.staticFunction(GenomeLoc.class, "parseGenomeLoc", String.class);
        Function1<GenomeLoc, String> f1 = parseOne.f1();
        Collection<GenomeLoc> result = Functions.map(f1, Arrays.asList(str.split(";")));
        GenomeLoc[] locs = (GenomeLoc[])result.toArray(new GenomeLoc[0]);

        Arrays.sort(locs);
        for ( GenomeLoc l : locs )
            System.out.printf("  -> %s%n", l);

        System.out.printf("  Locations are: %s%n", Utils.join(" ", Functions.map( Operators.toString, Arrays.asList(locs) ) ) );

        return locs;
    }

    /**
     * A key function that returns true if the proposed GenomeLoc curr is within the list of
     * locations we are processing in this TraversalEngine
     *
     * @param curr
     * @return true if we should process GenomeLoc curr, otherwise false
     */
    public boolean inLocations( GenomeLoc curr ) {
        if ( this.locs == null )
            return true;
        else {
            for ( GenomeLoc loc : this.locs ) {
                //System.out.printf("  Overlap %s vs. %s => %b%n", loc, curr, loc.overlapsP(curr));
                if ( loc.overlapsP(curr) )
                    return true;
            }
            return false;
        }
    }

    /**
     * Returns true iff we have a specified series of locations to process AND we are past the last
     * location in the list.  It means that, in a serial processing of the genome, that we are done.
     *
     * @param curr Current genome Location
     * @return true if we are past the last location to process
     */
    private boolean pastFinalLocation( GenomeLoc curr ) {
        boolean r = locs != null && locs[locs.length-1].compareTo( curr ) == -1 && ! locs[locs.length-1].overlapsP(curr);
        //System.out.printf("  pastFinalLocation %s vs. %s => %d => %b%n", locs[locs.length-1], curr, locs[locs.length-1].compareTo( curr ), r);
        return r;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // printing
    //
    // --------------------------------------------------------------------------------------------------------------
    
    /**
     *
     * @param curTime (current runtime, in millisecs)
     * @return true if the maximum interval (in millisecs) has passed since the last printing
     */
    private boolean maxElapsedIntervalForPrinting(final long curTime) {
        return (curTime - this.lastProgressPrintTime) > MAX_PROGRESS_PRINT_TIME;
    }

    /**
     * Forward request to printProgress
     *
     * @param type
     * @param loc
     */
    public void printProgress(final String type, GenomeLoc loc) {
        printProgress( false, type, loc );
    }

    /**
     * Utility routine that prints out process information (including timing) every N records or
     * every M seconds, for N and M set in global variables.
     *
     * @param mustPrint If true, will print out info, regardless of nRecords or time interval
     * @param type String to print out describing our atomic traversal type ("read", "locus", etc)
     * @param loc Current location 
     */
    public void printProgress( boolean mustPrint, final String type, GenomeLoc loc ) {
        // If an index is enabled, file read progress is meaningless because a linear
        // traversal is not being performed.  For now, don't bother printing progress.
        // TODO: Create a sam indexed read tracker that tracks based on percentage through the query.
        if( samReadingTracker == null )
            return;

        final long nRecords = this.nRecords;
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;
        //System.out.printf("Cur = %d, last print = %d%n", curTime, lastProgressPrintTime);
                    
        if ( mustPrint || nRecords % N_RECORDS_TO_PRINT == 0 || maxElapsedIntervalForPrinting(curTime)) {
            this.lastProgressPrintTime = curTime;
            final double secsPer1MReads = (elapsed * 1000000.0) / nRecords;
            if ( loc != null )
                System.out.printf("[PROGRESS] Traversed to %s, processing %,d %s in %.2f secs (%.2f secs per 1M %s)%n", loc, nRecords, type, elapsed, secsPer1MReads, type);
            else
                System.out.printf("[PROGRESS] Traversed %,d %s in %.2f secs (%.2f secs per 1M %s)%n", nRecords, type, elapsed, secsPer1MReads, type);

            // Currently samReadingTracker will print misleading info if we're not processing the whole file
            if ( this.locs == null )
                System.out.printf("[PROGRESS]   -> %s%n", samReadingTracker.progressMeter());
        }
    }

    /**
     * Called after a traversal to print out information about the traversal process
     *
     * @param type String describing this type of traversal ("loci", "read")
     * @param sum The reduce result of the traversal
     * @param <T> ReduceType of the traversal
     */
    protected <T> void printOnTraversalDone( final String type, T sum ) {
        printProgress( true, type, null );
        System.out.println("Traversal reduce result is " + sum);
        System.out.printf("Traversal skipped %d reads out of %d total (%.2f%%)%n", nSkippedReads, nReads, (nSkippedReads * 100.0) / nReads);
        System.out.printf("  -> %d unmapped reads%n", nUnmappedReads );
        System.out.printf("  -> %d non-primary reads%n", nNotPrimary );
        System.out.printf("  -> %d reads with bad alignments%n", nBadAlignments );
        System.out.printf("  -> %d reads with indels%n", nSkippedIndels );
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
    public boolean initialize(final boolean THREADED_IO) {
        lastProgressPrintTime = startTime = System.currentTimeMillis();
        initializeReference();
        initializeReads(THREADED_IO);
        // Initial the reference ordered data iterators
        initializeRODs();

        //testReference();
        //loadReference();
        return true;
    }

    private void initializeReads(final boolean THREADED_IO) {

        Iterator<SAMRecord> samIterator;
        try {
            samReadIter = loadSAMFile( readsFile, THREADED_IO );
        }
        catch( IOException ex ) {
            // TODO: IOException should be a checked exception in this case.
            throw new RuntimeIOException(ex);
        }

        if ( beSafeP )
            samReadIter = new VerifyingSamIterator(samReadIter);

        if ( THREADED_IO ) {
            System.out.printf("Enabling threaded I/O with buffer of %d reads%n", THREADED_IO_BUFFER_SIZE);
            samReadIter = new ThreadedIterator<SAMRecord>(samReadIter, THREADED_IO_BUFFER_SIZE);
        }
    }

    protected Iterator<SAMRecord> loadSAMFile( final File samFile, final boolean threadedIO )
            throws IOException {
        Iterator<SAMRecord> iterator = null;

        samReader = new SAMFileReader(readsFile, true);
        samReader.setValidationStringency(strictness);

        final SAMFileHeader header = samReader.getFileHeader();
        System.err.println("Sort order is: " + header.getSortOrder());

        // If the file has an index, querying functions are available.  Use them if possible...
        if(samReader.hasIndex()) {
            iterator = new SamQueryIterator( samReader, locs );
        }
        else {
            // Ugh.  Close and reopen the file so that the file progress decorator can be assigned to the input stream.
            samReader.close();

            final FileInputStream samFileStream = new FileInputStream(readsFile);
            final InputStream bufferedStream= new BufferedInputStream(samFileStream);
            samReader = new SAMFileReader(readsFile, true);
            samReader.setValidationStringency(strictness);

            samReadingTracker = new FileProgressTracker<SAMRecord>( readsFile, samReader.iterator(), samFileStream.getChannel(), 1000 );
            iterator = samReadingTracker;
        }

        return iterator;
    }


    /**
     * Prepare the reference for stream processing
     *
     */
    protected void initializeReference() {
        if ( refFileName != null ) {
            this.refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFileName);
            this.refIter = new ReferenceIterator(this.refFile);
            if ( ! Utils.setupRefContigOrdering(this.refFile) ) {
                // We couldn't process the reference contig ordering, fail since we need it
                Utils.scareUser(String.format("We couldn't load the contig dictionary associated with %s.  At the current time we require this dictionary file to efficiently access the FASTA file.  In the near future this program will automatically construct the dictionary for you and save it down.", refFileName));
            }
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
        for ( ReferenceOrderedData data : rods ) {
            rodIters.add(data.iterator());
        }
        return rodIters;
    }

    protected void testReference() {
        while (true) {
            ReferenceSequence ref = refFile.nextSequence();
            System.out.printf("%s %d %d%n", ref.getName(), ref.length(), System.currentTimeMillis());
            printProgress(true, "loci", new GenomeLoc("foo", 1));
        }
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
     * @param loc The location to get the rods at
     * @return A list of ReferenceOrderDatum at loc.  ROD without a datum at loc will be null in the list
     */
    protected List<ReferenceOrderedDatum> getReferenceOrderedDataAtLocus(List<ReferenceOrderedData.RODIterator> rodIters,
                                                                        final GenomeLoc loc) {
        List<ReferenceOrderedDatum> data = new ArrayList<ReferenceOrderedDatum>();
        for ( ReferenceOrderedData.RODIterator iter : rodIters ) {
            data.add(iter.seekForward(loc));
        }
        return data;
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
     *
     */
    class locusStreamFilterFunc implements SamRecordFilter {
        public boolean filterOut(SAMRecord rec) {
            boolean result = false;
            String why = "";
            if ( rec.getReadUnmappedFlag() ) {
                nUnmappedReads++;
                result = true;
                why = "Unmapped";
            }
            else if ( rec.getNotPrimaryAlignmentFlag() ) {
                nNotPrimary++;
                result = true;
                why = "Not Primary";
            }
            else if ( rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ) {
                nBadAlignments++;
                result = true;
                why = "No alignment start";
            }
            else {
                result = false;
            }

            if ( result ) {
                nSkippedReads++;
                //System.out.printf("  [filter] %s => %b %s%n", rec.getReadName(), result, why);
            }
            else {
                nReads++;
            }
            return result;        
        }
    }

    public void verifySortOrder(final boolean requiresSortedOrder) {
        if ( beSafeP && samReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate ) {
            final String msg = "SAM file is not sorted in coordinate order (according to header) Walker type with given arguments requires a sorted file for correct processing";
            if ( requiresSortedOrder || strictness == SAMFileReader.ValidationStringency.STRICT )
                throw new RuntimeIOException(msg);
            else if ( strictness == SAMFileReader.ValidationStringency.LENIENT )
                Utils.warnUser(msg);
        }
    }

    /**
     * Traverse by loci -- the key driver of linearly ordered traversal of loci.  Provides reads, RODs, and
     * the reference base for each locus in the reference to the LocusWalker walker.  Supports all of the
     * interaction contract implied by the locus walker
     *
     * @param walker A locus walker object
     * @param <M> MapType -- the result of calling map() on walker
     * @param <T> ReduceType -- the result of calling reduce() on the walker
     * @return 0 on success
     */
    protected <M,T> int traverseByLoci(LocusWalker<M,T> walker) {
        verifySortOrder(true);

        // prepare the read filtering read iterator and provide it to a new locus iterator
        FilteringIterator filterIter = new FilteringIterator(samReadIter, new locusStreamFilterFunc());
        //LocusIterator iter = new SingleLocusIterator(filterIter);
        LocusIterator iter = new LocusIteratorByHanger(filterIter);

        // initialize the walker object
        walker.initialize();
        // Initialize the T sum using the walker
        T sum = walker.reduceInit();
        boolean done = false;
        GenomeLoc prevLoc = null;

        while ( iter.hasNext() && ! done ) {
            this.nRecords++;

            // actually get the read and hand it to the walker
            final LocusContext locus = iter.next();

            // Poor man's version of index LOL
            GenomeLoc curLoc = locus.getLocation();
            if ( inLocations(curLoc) ) {
                if ( prevLoc != null && curLoc.compareContigs(prevLoc) != 0 )
                    System.out.printf("Traversing to next chromosome...%n");

                // Jump forward in the reference to this locus location
                final ReferenceIterator refSite = refIter.seekForward(locus.getLocation());
                final char refBase = refSite.getBaseAsChar();
                locus.setReferenceContig(refSite.getCurrentContig());

                // Iterate forward to get all reference ordered data covering this locus
                final List<ReferenceOrderedDatum> rodData = getReferenceOrderedDataAtLocus(rodIters, locus.getLocation());

                if ( DEBUGGING )
                    System.out.printf("  Reference: %s:%d %c%n", refSite.getCurrentContig().getName(), refSite.getPosition(), refBase);

                //
                // Execute our contract with the walker.  Call filter, map, and reduce
                //
                final boolean keepMeP = walker.filter(rodData, refBase, locus);
                if ( keepMeP ) {
                    M x = walker.map(rodData, refBase, locus);
                    sum = walker.reduce(x, sum);
                }

                if ( this.maxReads > 0 && this.nRecords > this.maxReads ) {
                    System.out.println("Maximum number of reads encountered, terminating traversal " + this.nRecords);
                    done = true;
                }


            }

            printProgress("loci", locus.getLocation());
            if ( pastFinalLocation(locus.getLocation()) )
                done = true;
        }

        printOnTraversalDone("loci", sum);
        walker.onTraversalDone();
        return 0;
    }

    /**
     * Traverse by read -- the key driver of linearly ordered traversal of reads.  Provides a single read to
     * the walker object, in coordinate order.  Supports all of the
     * interaction contract implied by the read walker
     *                                  sor
     * @param walker A read walker object
     * @param <M> MapType -- the result of calling map() on walker
     * @param <T> ReduceType -- the result of calling reduce() on the walker
     * @return 0 on success
     */
    protected <M,R> int traverseByRead(ReadWalker<M,R> walker) {
        if ( refFileName == null && ! walker.requiresOrderedReads() && verifyingSamReadIter != null ) {
            System.out.println("STATUS: No reference file provided and unordered reads are tolerated, enabling out of order read processing.");
            if ( verifyingSamReadIter != null )
                verifyingSamReadIter.setCheckOrderP(false);
        }

        verifySortOrder( refFileName != null || walker.requiresOrderedReads());

        // Initialize the walker
        walker.initialize();

        // Initialize the sum
        R sum = walker.reduceInit();
        List<Integer> offsets = Arrays.asList(0);   // Offset of a single read is always 0

        boolean done = false;
        while ( samReadIter.hasNext() && ! done ) {
            this.nRecords++;

            // get the next read
            final SAMRecord read = samReadIter.next();
            final List<SAMRecord> reads = Arrays.asList(read);
            GenomeLoc loc = Utils.genomicLocationOf(read);

            // Jump forward in the reference to this locus location
            LocusContext locus = new LocusContext(loc, reads, offsets);
            if ( ! loc.isUnmapped() && refIter != null ) {
                final ReferenceIterator refSite = refIter.seekForward(loc);
                locus.setReferenceContig(refSite.getCurrentContig());
            }

            if ( inLocations(loc) ) {

                //
                // execute the walker contact
                //
                final boolean keepMeP = walker.filter(locus, read);
                if ( keepMeP ) {
                    M x = walker.map(locus, read);
                    sum = walker.reduce(x, sum);
                }

                if ( this.maxReads > 0 && this.nRecords > this.maxReads ) {
                    System.out.println("Maximum number of reads encountered, terminating traversal " + this.nRecords);
                    done = true;
                }
            }
            printProgress("reads", loc);

            if ( pastFinalLocation(loc) )
                done = true;
            //System.out.printf("Done? %b%n", done);
         }

        printOnTraversalDone("reads", sum);
        walker.onTraversalDone();
        return 0;
    }
}

