package edu.mit.broad.sting.atk;

import edu.mit.broad.sam.*;
import edu.mit.broad.sam.SAMFileReader.ValidationStringency;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.util.RuntimeIOException;
import edu.mit.broad.picard.filter.SamRecordFilter;
import edu.mit.broad.picard.filter.FilteringIterator;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.sting.utils.ReferenceIterator;
import edu.mit.broad.sting.utils.ReferenceOrderedData;
import edu.mit.broad.sting.utils.ReferenceOrderedDatum;
import edu.mit.broad.sting.utils.Utils;

import java.io.*;
import java.util.List;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Arrays;

public class TraversalEngine {
    // Usage and parameters
    private File readsFile = null;
    private File refFileName = null;
    private List<ReferenceOrderedData> rods = null;

    private String regionStr = null;
    private String traversalType = null;
    private ValidationStringency strictness = ValidationStringency.STRICT;

    private long startTime = -1;
    private long maxReads = -1;
    private long nRecords = 0;
    private SAMFileReader samReader = null;
    private ReferenceSequenceFile refFile = null;
    private ReferenceIterator refIter = null;
    private SAMFileReader readStream;

    private int nReads = 0;
    private int nSkippedReads = 0;
    private int nUnmappedReads = 0;
    private int nNotPrimary = 0;
    private int nBadAlignments = 0;
    private int nSkippedIndels = 0;

    public boolean DEBUGGING = false;

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
    
    public void setRegion(final String reg) { regionStr = regionStr; }
    public void setTraversalType(final String type) { traversalType = type; }
    public void setStrictness( final ValidationStringency s ) { strictness = s; }
    public void setMaxReads( final int maxReads ) { this.maxReads = maxReads; }
    public void setDebugging( final boolean d ) { DEBUGGING = d; }

    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing with the reference sequence
    //
    // --------------------------------------------------------------------------------------------------------------
    public void printProgress(final String type) { printProgress( false, type ); }

    public void printProgress( boolean mustPrint, final String type ) {
        final long nRecords = this.nRecords;

        if ( mustPrint || nRecords % 100000 == 0 ) {
            final double elapsed = (System.currentTimeMillis() - startTime) / 1000.0;
            final double secsPer1MReads = (elapsed * 1000000.0) / nRecords;
            System.out.printf("Traversed %d %s %.2f secs (%.2f secs per 1M %s)%n", nRecords, type, elapsed, secsPer1MReads, type);
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
                                                                        final String contig, final int pos) {
        List<ReferenceOrderedDatum> data = new ArrayList<ReferenceOrderedDatum>();
        for ( ReferenceOrderedData.RODIterator iter : rodIters ) {
            data.add(iter.seekForward(contig, pos));
        }
        return data;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // traversal functions
    //
    // --------------------------------------------------------------------------------------------------------------
    protected int initialize() {
        startTime = System.currentTimeMillis();
        loadReference();
        //testReference();
        //loadReference();
        readStream = initializeReadStreams();
        return 0;
    }

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
        FilteringIterator filterIter = new FilteringIterator(readStream.iterator(), new locusStreamFilterFunc());
        CloseableIterator<LocusIterator> iter = new LocusIterator(filterIter);

        List<ReferenceOrderedData.RODIterator> rodIters = initializeRODs();

        T sum = walker.reduceInit();
        while ( iter.hasNext() ) {
            this.nRecords++;

            // actually get the read and hand it to the walker
            final LocusIterator locus = iter.next();
            final ReferenceIterator refSite = refIter.seekForward(locus.getContig(), locus.getPosition());
            final char refBase = refSite.getBaseAsChar();
            final List<ReferenceOrderedDatum> rodData = getReferenceOrderedDataAtLocus(rodIters, locus.getContig(), locus.getPosition());

            if ( DEBUGGING )
                System.out.printf("  Reference: %s:%d %c%n", refSite.getCurrentContig().getName(), refSite.getPosition(), refBase);

            final boolean keepMeP = walker.filter(rodData, refBase, locus);
            if ( keepMeP ) {
                M x = walker.map(rodData, refBase, locus);
                sum = walker.reduce(x, sum);
            }

            if ( this.maxReads > 0 && this.nRecords > this.maxReads ) {
                System.out.println("Maximum number of reads encountered, terminating traversal " + this.nRecords);
                break;
            }

            printProgress("loci");
        }

        printProgress( true, "loci" );
        System.out.println("Traversal reduce result is " + sum);
        System.out.printf("Traversal skipped %d reads out of %d total (%.2f%%)%n", nSkippedReads, nReads, (nSkippedReads * 100.0) / nReads);
        System.out.printf("  -> %d unmapped reads%n", nUnmappedReads );
        System.out.printf("  -> %d non-primary reads%n", nNotPrimary );
        System.out.printf("  -> %d reads with bad alignments%n", nBadAlignments );
        System.out.printf("  -> %d reads with indels%n", nSkippedIndels );
        walker.onTraveralDone();
        return 0;
    }

    protected <M,R> int traverseByRead(ReadWalker<M,R> walker) {
        walker.initialize();
        CloseableIterator<SAMRecord> iter = readStream.iterator();
        R sum = walker.reduceInit();
        while ( iter.hasNext() ) {
            this.nRecords++;

            // actually get the read and hand it to the walker
            final SAMRecord read = iter.next();
            final boolean keepMeP = walker.filter(null, read);
            if ( keepMeP ) {
                M x = walker.map(null, read);
                sum = walker.reduce(x, sum);
            }

            if ( this.maxReads > 0 && this.nRecords > this.maxReads ) {
                System.out.println("Maximum number of reads encountered, terminating traversal " + this.nRecords);
                break;
            }

            printProgress("reads");
        }

        printProgress( true, "reads" );
        System.out.println("Traversal reduce result is " + sum);
        walker.onTraveralDone();
        return 0;
    }

    //
    //
    // Prepare the input streams
    //
    //
    private SAMFileReader initializeReadStreams() {
        SAMFileReader reader = getSamReader(readsFile);
        return reader;
    }

    private SAMFileReader getSamReader(final File samFile) {
        try {
            final InputStream samInputStream = new FileInputStream(samFile);
            final InputStream bufferedStream= new BufferedInputStream(samInputStream);
            //final InputStream bufferedStream= new BufferedInputStream(samInputStream, 10000000);
            final SAMFileReader samReader = new SAMFileReader(bufferedStream, true);
            samReader.setValidationStringency(strictness);

            final SAMFileHeader header = samReader.getFileHeader();
            System.err.println("Sort order is: " + header.getSortOrder());

            return samReader;
        }
        catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }
}