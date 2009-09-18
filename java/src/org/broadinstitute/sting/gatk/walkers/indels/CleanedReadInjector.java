package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
/**
 * User: hanna                                                                                                                    
 * Date: Jun 10, 2009
 * Time: 2:40:19 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Copies reads from the input stream into the <code>outputBAM</code>, replacing those
 * reads which have been cleaned with their new clean copies.
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CleanedReadInjector extends ReadWalker<Integer,Integer> {

    /**
     * The source of all cleaned reads.
     */
    @Argument(fullName="cleaned_reads",shortName="cr",doc="BAM file of reads which have been cleaned.",required=true)
    SAMFileReader cleanedReadsSource = null;

    /**
     * Target file for BAM output.
     */
    @Argument(fullName="output_bam",shortName="ob",doc="Output BAM file",required=true)
    SAMFileWriter outputBAM = null;

    /**
     * The set of (sorted) cleaned reads
     */
    private Queue<SAMRecord> cleanedReads = new LinkedList<SAMRecord>();

    /**
     * The intervals specified by the user
     */
    private HashMap<String, ArrayList<GenomeLoc>> intervals = null;

    /**
     * A fast lookup table for uniquified read info
     */
    private HashSet<String> cleanedReadHash = new HashSet<String>();

    @Override
    public void initialize() {

        // For now, read the whole damn file into memory.  If this becomes a problem,
        // then we just need to read the hash into memory and the first read; we'd then
        // need to query the BAM file every time we needed to update the cleaned read iterator.
        CloseableIterator<SAMRecord> allReads = cleanedReadsSource.iterator();
        while ( allReads.hasNext() ) {
            SAMRecord read = allReads.next();
            cleanedReads.add(read);
            cleanedReadHash.add(getUniquifiedReadName(read));
        }
        allReads.close();

	// If there are intervals specified by the user,record them so we can make sure not
	// to emit reads outside the intervals.  For now, we'll group them by chromosome to
	// make lookup a bit faster.
        if ( this.getToolkit() != null &&
	     this.getToolkit().getArguments().intervals != null ) {
	    intervals = new HashMap<String, ArrayList<GenomeLoc>>();
	    List<GenomeLoc> locs = GenomeAnalysisEngine.parseIntervalRegion(this.getToolkit().getArguments().intervals);
	    Iterator<GenomeLoc> iter = GenomeLocSortedSet.createSetFromList(locs).iterator();
	    while ( iter.hasNext() ) {
		GenomeLoc loc = iter.next();
		if ( intervals.get(loc.getContig()) == null )
		    intervals.put(loc.getContig(), new ArrayList<GenomeLoc>());
		intervals.get(loc.getContig()).add(loc);
	    }
        }
    }

    /**
     * For every read, if it exists in the cleaned read set, ignore it; otherwise, emit it.
     * Also, if the head of the cleaned read list could be emitted here, do so.
     * @param ref Portion of the reference sequence overlapping this read.
     * @param read The read to inspect.
     * @return Number of reads cleaned by this iteration of the map function.
     */
    @Override
    public Integer map(char[] ref, SAMRecord read) {

        // first emit reads from the cleaned set if appropriate
        int cleanedReadCount = 0;
        SAMRecord firstCleanedRead = cleanedReads.peek();
        while ( firstCleanedRead != null &&
                firstCleanedRead.getReferenceIndex() <= read.getReferenceIndex() &&
                firstCleanedRead.getAlignmentStart() <= read.getAlignmentStart() ) {
	    if ( emit(firstCleanedRead) )
		cleanedReadCount++;
	    cleanedReads.remove();
            firstCleanedRead = cleanedReads.peek();
        }

        if ( !cleanedReadHash.contains(getUniquifiedReadName(read)) )
            outputBAM.addAlignment(read);
        return cleanedReadCount;
    }

    /**
     * Determine whether to emit the given read; if so, return true.
     */
    private boolean emit(SAMRecord read) {
	// if no intervals were specified, emit everything
	if ( intervals == null ) {
	    outputBAM.addAlignment(read);
	    return true;
	}

	ArrayList<GenomeLoc> intervalList = intervals.get(read.getReferenceName());
	if ( intervalList == null )
	    return false;

	GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);
	for ( GenomeLoc interval : intervalList ) {
	    // if it overlaps an interval, then we can emit it
	    if ( interval.overlapsP(readLoc) ) {
		outputBAM.addAlignment(read);
		return true;
	    }

	    // once we've passed any interval that could overlap it, just quit
	    if ( interval.isPast(readLoc) )
		return false;
	}

	// it didn't overlap an interval
	return false;
    }

    /**
     * Initialize traversal with number of reads which have been replaced with a clean version.
     * @return 0 to initialize the traversal.
     */
    @Override
    public Integer reduceInit() { return 0; }

    /**
     * Accumulates the number of reads which have been replaced by their clean equivalents.
     * @param value Number of reads cleaned by this iteration.
     * @param accum Prior number of reads cleaned.
     * @return Total number of reads cleaned, including this iteration.
     */
    @Override
    public Integer reduce( Integer value, Integer accum ) {
        return accum + value;
    }

    @Override
    public void onTraversalDone( Integer value ) {
        outputBAM.close();
    }

    /**
     * Gets the distinguishing elements of the read name from the given read.
     * @param read read to uniquify.
     * @return A (hopefully) completely unique name for the read.
     */
    private static String getUniquifiedReadName( SAMRecord read ) {
        return String.format("%s.%s.%s.%s",read.getAttribute("RG"),read.getReadName(),read.getFlags(),read.getReadString());            
    }
}
