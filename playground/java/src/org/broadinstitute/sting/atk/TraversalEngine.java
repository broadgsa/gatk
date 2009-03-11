package org.broadinstitute.sting.atk;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.RuntimeIOException;
import edu.mit.broad.picard.filter.SamRecordFilter;
import edu.mit.broad.picard.filter.FilteringIterator;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.sting.utils.*;

import java.io.*;
import java.util.*;

import net.sf.functionalj.reflect.StdReflect;
import net.sf.functionalj.reflect.JdkStdReflect;
import net.sf.functionalj.FunctionN;
import net.sf.functionalj.Function1;
import net.sf.functionalj.Functions;
import net.sf.functionalj.util.Operators;

public class TraversalEngine {
    // Usage and parameters
    private File readsFile = null;
    private File refFileName = null;
    private List<ReferenceOrderedData> rods = null;

    private String regionStr = null;
    private String traversalType = null;
    private ValidationStringency strictness = ValidationStringency.STRICT;

    private long startTime = -1;
    private long lastProgressPrintTime = -1;
    private long MAX_PROGRESS_PRINT_TIME = 5 * 1000; // 10 seconds in millisecs
    private long maxReads = -1;
    private long nRecords = 0;
    private SAMFileReader samReader = null;
    private ReferenceSequenceFile refFile = null;
    private ReferenceIterator refIter = null;
    private SAMFileReader readStream;
    private Iterator<SAMRecord> samReadIter = null;

    private int nReads = 0;
    private int nSkippedReads = 0;
    private int nUnmappedReads = 0;
    private int nNotPrimary = 0;
    private int nBadAlignments = 0;
    private int nSkippedIndels = 0;
    private FileProgressTracker samReadingTracker = null;

    public boolean DEBUGGING = false;

    private GenomeLoc[] locs = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Setting up the engine
    //
    // --------------------------------------------------------------------------------------------------------------
    public TraversalEngine(File reads, File ref, ReferenceOrderedData[] rods ) {
        readsFile = reads;
        refFileName = ref;
        this.rods = Arrays.asList(rods);
    }

    protected int initialize() {
        lastProgressPrintTime = startTime = System.currentTimeMillis();
        loadReference();
        //testReference();
        //loadReference();
        try {
            final FileInputStream samFileStream = new FileInputStream(readsFile);
            final InputStream bufferedStream= new BufferedInputStream(samFileStream);
            //final InputStream bufferedStream= new BufferedInputStream(samInputStream, 10000000);
            final SAMFileReader samReader = new SAMFileReader(bufferedStream, true);
            samReader.setValidationStringency(strictness);

            final SAMFileHeader header = samReader.getFileHeader();
            System.err.println("Sort order is: " + header.getSortOrder());

            samReadingTracker = new FileProgressTracker<SAMRecord>( readsFile, samReader.iterator(), samFileStream.getChannel(), 1000 );
            samReadIter = samReadingTracker;
        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }

        return 0;
    }
    
    public void setRegion(final String reg) { regionStr = regionStr; }
    public void setTraversalType(final String type) { traversalType = type; }
    public void setStrictness( final ValidationStringency s ) { strictness = s; }
    public void setMaxReads( final int maxReads ) { this.maxReads = maxReads; }
    public void setDebugging( final boolean d ) { DEBUGGING = d; }


    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing locations (areas of the genome we're traversing over)
    //
    // --------------------------------------------------------------------------------------------------------------
    public void setLocation( final String locStr ) {
        this.locs = parseGenomeLocs(locStr);
    }

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

    public boolean pastFinalLocation( GenomeLoc curr ) {
        boolean r = locs != null && locs[locs.length-1].compareTo( curr ) == -1 && ! locs[locs.length-1].overlapsP(curr);
        //System.out.printf("  pastFinalLocation %s vs. %s => %d => %b%n", locs[locs.length-1], curr, locs[locs.length-1].compareTo( curr ), r);
        return r;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing with the reference sequence
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

    public void printProgress(final String type, GenomeLoc loc) { printProgress( false, type, loc ); }

    public void printProgress( boolean mustPrint, final String type, GenomeLoc loc ) {
        final long nRecords = this.nRecords;
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;
        //System.out.printf("Cur = %d, last print = %d%n", curTime, lastProgressPrintTime);
                    
        if ( mustPrint || nRecords % 100000 == 0 || maxElapsedIntervalForPrinting(curTime)) {
            this.lastProgressPrintTime = curTime;
            final double secsPer1MReads = (elapsed * 1000000.0) / nRecords;
            if ( loc != null )
                System.out.printf("[PROGRESS] Traversed to %s, processing %d %s %.2f secs (%.2f secs per 1M %s)%n", loc, nRecords, type, elapsed, secsPer1MReads, type);
            else
                System.out.printf("[PROGRESS] Traversed %d %s %.2f secs (%.2f secs per 1M %s)%n", nRecords, type, elapsed, secsPer1MReads, type);

            if ( this.locs == null )
                System.out.printf("[PROGRESS]   -> %s%n", samReadingTracker.progressMeter());
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing with the reference sequence
    //
    // --------------------------------------------------------------------------------------------------------------

    protected void loadReference() {
        if ( refFileName!= null ) {
            this.refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFileName);
            this.refIter = new ReferenceIterator(this.refFile);
            Utils.setupRefContigOrdering(this.refFile);
        }         
    }

    protected void testReference() {
        String line = "";
        refIter.seekForward("chr20", 79);
        for ( int i = 0; i < this.maxReads && refIter.hasNext(); i++ ) {
            final ReferenceIterator refSite = refIter.next();
            final char refBase = refSite.getBaseAsChar();
            line += refBase;
            if ( (i + 1) % 80 == 0 ) {
                System.out.println(line);
                line = "";
            }
            //System.out.printf("  Reference: %s:%d %c%n", refSite.getCurrentContig().getName(), refSite.getPosition(), refBase);
        }
        System.out.println(line);
        System.exit(1);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // dealing with reference ordered data
    //
    // --------------------------------------------------------------------------------------------------------------

    protected List<ReferenceOrderedData.RODIterator> initializeRODs() {
        // set up reference ordered data
        List<ReferenceOrderedData.RODIterator> rodIters = new ArrayList<ReferenceOrderedData.RODIterator>();
        for ( ReferenceOrderedData data : rods ) {
            rodIters.add(data.iterator());
        }
        return rodIters;
    }

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
            else if ( rec.getCigar().numCigarElements() > 1 ) {
                // FIXME -- deal with indels correctly!
                nSkippedIndels++;
                result = true;
                why = "Skipping indel: " + rec.getCigarString();
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

    protected <M,T> int traverseByLoci(LocusWalker<M,T> walker) {
        walker.initialize();
        FilteringIterator filterIter = new FilteringIterator(samReadIter, new locusStreamFilterFunc());
        CloseableIterator<LocusIterator> iter = new LocusIterator(filterIter);

        List<ReferenceOrderedData.RODIterator> rodIters = initializeRODs();

        T sum = walker.reduceInit();
        boolean done = false;
        while ( iter.hasNext() && ! done ) {
            this.nRecords++;

            // actually get the read and hand it to the walker
            final LocusIterator locus = iter.next();

            // Poor man's version of index LOL
            if ( inLocations(locus.getLocation()) ) {
                final ReferenceIterator refSite = refIter.seekForward(locus.getContig(), locus.getPosition());
                final char refBase = refSite.getBaseAsChar();
                final List<ReferenceOrderedDatum> rodData = getReferenceOrderedDataAtLocus(rodIters, locus.getLocation());

                if ( DEBUGGING )
                    System.out.printf("  Reference: %s:%d %c%n", refSite.getCurrentContig().getName(), refSite.getPosition(), refBase);

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

        printProgress( true, "loci", null );
        System.out.println("Traversal reduce result is " + sum);
        System.out.printf("Traversal skipped %d reads out of %d total (%.2f%%)%n", nSkippedReads, nReads, (nSkippedReads * 100.0) / nReads);
        System.out.printf("  -> %d unmapped reads%n", nUnmappedReads );
        System.out.printf("  -> %d non-primary reads%n", nNotPrimary );
        System.out.printf("  -> %d reads with bad alignments%n", nBadAlignments );
        System.out.printf("  -> %d reads with indels%n", nSkippedIndels );
        walker.onTraveralDone();
        return 0;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // traversal by read functions
    //
    // --------------------------------------------------------------------------------------------------------------
    protected <M,R> int traverseByRead(ReadWalker<M,R> walker) {
        walker.initialize();

        R sum = walker.reduceInit();
        boolean done = false;
        while ( samReadIter.hasNext() && ! done ) {
            this.nRecords++;

            // actually get the read and hand it to the walker
            final SAMRecord read = samReadIter.next();
            GenomeLoc loc = new GenomeLoc(read.getReferenceName(), read.getAlignmentStart());

            if ( inLocations(loc) ) {
                final boolean keepMeP = walker.filter(null, read);

                if ( keepMeP ) {
                    M x = walker.map(null, read);
                    sum = walker.reduce(x, sum);
                }

                if ( this.maxReads > 0 && this.nRecords > this.maxReads ) {
                    System.out.println("Maximum number of reads encountered, terminating traversal " + this.nRecords);
                    break;
                }

                printProgress("reads", loc);
            }
            
            if ( pastFinalLocation(loc) )
                done = true;
            //System.out.printf("Done? %b%n", done);
         }

        printProgress( true, "reads", null );
        System.out.println("Traversal reduce result is " + sum);
        walker.onTraveralDone();
        return 0;
    }
}
