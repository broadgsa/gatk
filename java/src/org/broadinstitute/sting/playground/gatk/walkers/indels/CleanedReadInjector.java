package org.broadinstitute.sting.playground.gatk.walkers.indels;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.Utils;
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
public class CleanedReadInjector extends ReadWalker<Integer,Integer> {

    /**
     * The source of all cleaned intervals.
     */
    @Argument(fullName="cleaned_intervals",shortName="ci",doc="Intervals which have been cleaned.",required=true)
    String intervalsSource = null;

    /**
     * The source of all cleaned reads.
     */
    @Argument(fullName="cleaned_reads",shortName="cr",doc="BAM file of reads which have been cleaned.",required=true)
    SAMFileReader cleanedReadsSource = null;

    /**
     * Target file for BAM output.
     */
    @Argument(fullName="output_bam",shortName="ob",doc="Output BAM file",required=true)
    String outputBAMFileName = null;

    /**
     * A stream of processed intervals.  The current head of the queue represents the next interval.
     */
    private Queue<GenomeLoc> intervals;

    /**
     * The interval currently in process.
     */
    private GenomeLoc interval = null;

    /**
     * A structure for fast lookup of all reads that have been cleaned in the current interval.
     */
    private Map<String,SAMRecord> cleanedReads = new HashMap<String,SAMRecord>();

    /**
     * The writer that handles writing of SAM files. 
     */
    SAMFileWriter outputBAM = null;

    @Override
    public void initialize() {
        intervals = parseIntervals( intervalsSource );
        interval = intervals.remove();
        loadCleanedReadsOverlappingInterval( interval );

        // HACK: The unit tests create their own output files.  Make sure this walker doesn't step
        //       on any toes.
        if( outputBAM == null ) {
            outputBAM = Utils.createSAMFileWriterWithCompression(getToolkit().getEngine().getSAMHeader(),
                                                                 false,
                                                                 outputBAMFileName,
                                                                 getToolkit().getBAMCompression());
        }
    }

    /**
     * For every read, inspect the set of cleaned reads that could possibly replace the current one.
     * If one exists, add the cleaned read to the output.  Otherwise, add the current read.
     * @param ref Portion of the reference sequence overlapping this read.
     * @param read The read to inspect.
     * @return Number of reads cleaned by this iteration of the map function.
     */
    @Override
    public Integer map(char[] ref, SAMRecord read) {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(read);

        while( loc.isPast(interval) && intervals.size() > 0 ) {
            interval = intervals.remove();
            loadCleanedReadsOverlappingInterval(interval);
        }

        String uniquifiedReadName = getUniquifiedReadName(read);
        int numCleaned = 0;

        if( loc.overlapsP(interval) && cleanedReads.containsKey(uniquifiedReadName) ) {
            SAMRecord cleanedRead = cleanedReads.get(uniquifiedReadName);
            cleanedRead.setReadName(read.getReadName());
            outputBAM.addAlignment(cleanedRead);
            numCleaned = 1;
        }
        else {
            outputBAM.addAlignment(read);
            numCleaned = 0;
        }

        return numCleaned;
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
     * Load the intervals directly from the command-line or from file, as appropriate.
     * Merge overlapping intervals.
     * @param intervalsSource Source of intervals.
     * @return a queue of sorted, merged intervals.
     */
    private Queue parseIntervals( String intervalsSource ) {
        List<GenomeLoc> parsedIntervals = GenomeAnalysisEngine.parseIntervalRegion(intervalsSource,false);
        GenomeLocSortedSet intervalSortedSet = new GenomeLocSortedSet();
        for( GenomeLoc parsedInterval: parsedIntervals )
            intervalSortedSet.addRegion(parsedInterval);

        return new LinkedList<GenomeLoc>( intervalSortedSet );        
    }

    /**
     * Load a list of all the reads overlapping the given interval into memory.
     * @param interval
     */
    private void loadCleanedReadsOverlappingInterval( GenomeLoc interval ) {
        CloseableIterator<SAMRecord> overlappingReads = cleanedReadsSource.queryOverlapping( interval.getContig(), (int)interval.getStart(), (int)interval.getStop() );
        // Load in all reads mapped to this region.  The cleaner will augment the read name in a way that uniquifies it.
        while( overlappingReads.hasNext() ) {
            SAMRecord read = overlappingReads.next();
            cleanedReads.put( getUniquifiedReadName(read), read );
        }
        overlappingReads.close();
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
