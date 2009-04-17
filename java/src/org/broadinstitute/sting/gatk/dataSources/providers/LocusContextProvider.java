package org.broadinstitute.sting.gatk.dataSources.providers;

import edu.mit.broad.picard.filter.FilteringIterator;
import edu.mit.broad.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.iterators.LocusIteratorByHanger;
import org.broadinstitute.sting.gatk.traversals.TraversalStatistics;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 8, 2009
 * Time: 3:00:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class LocusContextProvider {
    private Iterator<SAMRecord> reads;

    // What's the last locus accessed?  Used for sanity checking.
    private GenomeLoc lastLoc = null;
    private LocusIterator loci;
    private LocusContext locus;
    protected static Logger logger = Logger.getLogger(LocusContextProvider.class);

    public LocusContextProvider( Iterator<SAMRecord> reads ) {
        this.reads = new FilteringIterator(reads, new TraversalEngine.locusStreamFilterFunc());
        // prepare the iterator by loci from reads
        loci = new LocusIteratorByHanger(this.reads);
    }

    public LocusContext getLocusContext( GenomeLoc loc ) {
        // Precondition checks
        if( lastLoc != null && !loc.isPast( lastLoc ) )
            throw new RuntimeException( "Internal error: LocusContextProvider assumes that queries it receives are ordered." );

        if( (loc.getStop() - loc.getStart()) > 0 )
            throw new RuntimeException( "Internal error :LocusContextProviders currently require 1-base genomeLocs.");

        // jump to the first reference site
        LocusContext locusContext = advanceReadsToLoc( loci, loc );

        // if no locus context was found, create an empty locus
        if ( locusContext == null || locusContext.getLocation().compareTo( loc ) != 0 )
            locusContext = new LocusContext(loc, new ArrayList<SAMRecord>(), new ArrayList<Integer>());

        lastLoc = loc;

        return locusContext;
    }

    private LocusContext advanceReadsToLoc(LocusIterator locusIter, GenomeLoc target) {
        if ( ! locusIter.hasNext() )
            return null;

        if (locus == null) {
            locus = locusIter.next();
        }

        while (target.isPast(locus.getLocation()) && locusIter.hasNext() ) {
            logger.debug(String.format("  current locus is %s vs %s => %d", locus.getLocation(), target, locus.getLocation().compareTo(target)));
            locus = locusIter.next();
        }

        logger.debug(String.format("  returning %s", locus.getLocation()));
        return locus;
    }

}
