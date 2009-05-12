package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.LocusWindowWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;

import java.util.*;
import java.io.File;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Apr 23, 2009
 * Time: 10:26:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class TraverseByLocusWindows extends TraversalEngine {

    public TraverseByLocusWindows(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super(reads, ref, rods);
    }

    @Override
    public <M,T> T traverse(Walker<M,T> walker, List<GenomeLoc> locations) {
        if ( walker instanceof LocusWindowWalker ) {
            LocusWindowWalker<M, T> locusWindowWalker = (LocusWindowWalker<M, T>)walker;
            T sum = traverseByIntervals(locusWindowWalker, locations);
            return sum;
        } else {
            throw new IllegalArgumentException("Walker isn't an interval walker!");
        }
    }

    /**
     * Traverse by intervals -- the key driver of linearly ordered traversal of intervals.  Provides reads, RODs, and
     * the reference base for each interval in the reference to the intervalWalker walker.  Supports all of the
     * interaction contract implied by the interval walker
     *
     * @param walker An interval walker object
     * @param <M>    MapType -- the result of calling map() on walker
     * @param <T>    ReduceType -- the result of calling reduce() on the walker
     * @return 0 on success
     */
    protected <M, T> T traverseByIntervals(LocusWindowWalker<M, T> walker, List<GenomeLoc> locations) {
        logger.debug("Entering traverseByIntervals");

        if(readsFiles.size() > 1)
            throw new UnsupportedOperationException("Cannot do ByInterval traversal on file with multiple inputs.");        

        samReader = initializeSAMFile(readsFiles.get(0));

        verifySortOrder(true);

        walker.initialize();

        T sum = walker.reduceInit();

        if ( locations.isEmpty() ) {
            logger.debug("There are no intervals provided for the traversal");
        } else {
            if ( !samReader.hasIndex() )
                Utils.scareUser("Processing locations were requested, but no index was found for the input SAM/BAM file. This operation is potentially dangerously slow, aborting.");

            if ( walker.actOnNonIntervalReads() )
                sum = fullInputTraversal(walker, locations, sum);
            else
                sum = strictIntervalTraversal(walker, locations, sum);            
        }

        //printOnTraversalDone("intervals", sum);
        walker.onTraversalDone(sum);
        return sum;
    }

    protected <M, T> T strictIntervalTraversal(LocusWindowWalker<M, T> walker, List<GenomeLoc> locations, T sum) {
        LocusContext nextLocusToCarry = null;
        for ( GenomeLoc interval : locations ) {
            logger.debug(String.format("Processing interval %s", interval.toString()));

            CloseableIterator<SAMRecord> readIter = samReader.queryOverlapping( interval.getContig(),
                    (int)interval.getStart(),
                    (int)interval.getStop());
            Iterator<SAMRecord> wrappedIter = WrapReadsIterator(readIter, false);
            LocusContext locus = getLocusContext(wrappedIter, interval);
            readIter.close();

            if ( nextLocusToCarry == null ) {
                nextLocusToCarry = locus;
            } else if ( nextLocusToCarry.getLocation().overlapsP(locus.getLocation()) ) {
                nextLocusToCarry = merge(nextLocusToCarry, locus);
            } else {
                sum = carryWalkerOverInterval(walker, sum, nextLocusToCarry);
                nextLocusToCarry = locus;
            }
        }
        if ( nextLocusToCarry != null )
            sum = carryWalkerOverInterval(walker, sum, nextLocusToCarry);
        return sum;
    }

    protected <M, T> T fullInputTraversal(LocusWindowWalker<M, T> walker, List<GenomeLoc> locations, T sum) {
        ArrayList<LocusContext> nextLociToCarry = new ArrayList<LocusContext>();

        // set everything up
        GenomeLoc currentInterval = (locations.size() > 0 ? locations.get(0) : null);
        int locationsIndex = 0;
        ArrayList<SAMRecord> intervalReads = new ArrayList<SAMRecord>();
        Iterator<SAMRecord> readIter = samReader.iterator();

        while (readIter.hasNext()) {
            TraversalStatistics.nRecords++;
            SAMRecord read = readIter.next();

            // if there are no locations or we're past the last one or it's unmapped, then act on the read separately
            if ( currentInterval == null || read.getReadUnmappedFlag() ) {
                if ( nextLociToCarry.size() > 0 ) {
                    sum = carryWalkerOverInterval(walker, sum, nextLociToCarry.get(0));
                    for (int i=1; i < nextLociToCarry.size(); i++)
                        walker.nonIntervalReadAction(nextLociToCarry.get(i).getReads().get(0));
                    nextLociToCarry.clear();
                }
                walker.nonIntervalReadAction(read);
            }
            else {
                GenomeLoc loc = new GenomeLoc(read);
                // if we're in the current interval, add it to the list
                if ( currentInterval.overlapsP(loc) ) {
                    intervalReads.add(read);
                }
                // if we're not yet in the interval, act on the read separately
                else if ( currentInterval.isPast(loc) ) {
                    if ( nextLociToCarry.size() == 0 ) {
                        walker.nonIntervalReadAction(read);
                    } else {
                        ArrayList<SAMRecord> list = new ArrayList<SAMRecord>();
                        list.add(read);
                        nextLociToCarry.add(new LocusContext(loc, list, null));
                    }                        
                }
                // otherwise, we're past the interval so first deal with the collected reads and then this one
                else {
                    if ( intervalReads.size() > 0 ) {
                        Iterator<SAMRecord> wrappedIter = WrapReadsIterator(intervalReads.iterator(), false);
                        LocusContext locus = getLocusContext(wrappedIter, currentInterval);

                        if ( nextLociToCarry.size() == 0 ) {
                            nextLociToCarry.add(locus);
                        } else if ( nextLociToCarry.get(0).getLocation().overlapsP(locus.getLocation()) ) {
                            LocusContext newLocus = merge(nextLociToCarry.get(0), locus);
                            for (int i=1; i < nextLociToCarry.size(); i++)
                                newLocus = merge(newLocus, nextLociToCarry.get(i));    
                            nextLociToCarry.clear();
                            nextLociToCarry.add(newLocus);
                        } else {
                            sum = carryWalkerOverInterval(walker, sum, nextLociToCarry.get(0));
                            for (int i=1; i < nextLociToCarry.size(); i++)
                                walker.nonIntervalReadAction(nextLociToCarry.get(i).getReads().get(0));
                            nextLociToCarry.clear();
                            nextLociToCarry.add(locus);
                        }

                        // then prepare for the next interval
                        intervalReads.clear();
                    }
                    currentInterval = (++locationsIndex < locations.size() ? locations.get(locationsIndex) : null);

                    if ( nextLociToCarry.size() == 0 ) {
                        walker.nonIntervalReadAction(read);
                    } else {
                        ArrayList<SAMRecord> list = new ArrayList<SAMRecord>();
                        list.add(read);
                        nextLociToCarry.add(new LocusContext(loc, list, null));
                    }
                }
            }
        }
        // some cleanup
        if ( intervalReads.size() > 0 ) {
            Iterator<SAMRecord> wrappedIter = WrapReadsIterator(intervalReads.iterator(), false);
            LocusContext locus = getLocusContext(wrappedIter, currentInterval);
            if ( nextLociToCarry.size() == 0 ) {
                nextLociToCarry.add(locus);
            } else if ( nextLociToCarry.get(0).getLocation().overlapsP(locus.getLocation()) ) {
                LocusContext newLocus = merge(nextLociToCarry.get(0), locus);
                for (int i=1; i < nextLociToCarry.size(); i++)
                    newLocus = merge(newLocus, nextLociToCarry.get(i));
                nextLociToCarry.clear();
                nextLociToCarry.add(newLocus);
            } else {
                sum = carryWalkerOverInterval(walker, sum, nextLociToCarry.get(0));
                for (int i=1; i < nextLociToCarry.size(); i++)
                    walker.nonIntervalReadAction(nextLociToCarry.get(i).getReads().get(0));
                nextLociToCarry.clear();
                nextLociToCarry.add(locus);
            }
        }
        if ( nextLociToCarry.size() > 0 ) {
            sum = carryWalkerOverInterval(walker, sum, nextLociToCarry.get(0));
            for (int i=1; i < nextLociToCarry.size(); i++)
                walker.nonIntervalReadAction(nextLociToCarry.get(i).getReads().get(0));
        }

        return sum;
    }

    private LocusContext getLocusContext(Iterator<SAMRecord> readIter, GenomeLoc interval) {
        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        boolean done = false;
        long leftmostIndex = interval.getStart(),
                rightmostIndex = interval.getStop();
        while (readIter.hasNext() && !done) {
            TraversalStatistics.nRecords++;
            SAMRecord read = readIter.next();
            reads.add(read);
            if ( read.getAlignmentStart() < leftmostIndex )
                leftmostIndex = read.getAlignmentStart();
            if ( read.getAlignmentEnd() > rightmostIndex )
                rightmostIndex = read.getAlignmentEnd();
            if ( this.maxReads > 0 && TraversalStatistics.nRecords > this.maxReads ) {
                logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords));
                done = true;
            }
        }

        GenomeLoc window = new GenomeLoc(interval.getContig(), leftmostIndex, rightmostIndex);
        LocusContext locus = new LocusContext(window, reads, null);
        if ( DOWNSAMPLE_BY_COVERAGE )
            locus.downsampleToCoverage(downsamplingCoverage);

        return locus;
    }

    private LocusContext merge(LocusContext locus1, LocusContext locus2) {
        GenomeLoc loc = locus1.getLocation().merge(locus2.getLocation());
        TreeSet<SAMRecord> set = new TreeSet<SAMRecord>(new Comparator<SAMRecord>() {
            public int compare(SAMRecord obj1, SAMRecord obj2) {
                GenomeLoc myLoc = new GenomeLoc(obj1);
                GenomeLoc hisLoc = new GenomeLoc(obj2);
                int comparison = myLoc.compareTo(hisLoc);
                 // if the reads have the same start position, we must give a non-zero comparison
                 // (because java Sets often require "consistency with equals")
                if ( comparison == 0 )
                    comparison = obj1.getReadName().compareTo(obj2.getReadName());
                return comparison;
            }
        });

        set.addAll(locus1.getReads());
        set.addAll(locus2.getReads());
        return new LocusContext(loc, new ArrayList<SAMRecord>(set), null);
    }

    protected <M, T> T carryWalkerOverInterval(LocusWindowWalker<M, T> walker, T sum, LocusContext window) {
        ReferenceIterator refSite = refIter.seekForward(window.getLocation());
        StringBuffer refBases = new StringBuffer(refSite.getBaseAsString());
        int locusLength = (int)(window.getLocation().getStop() - window.getLocation().getStart());
        for ( int i = 0; i < locusLength; i++ ) {
            refSite = refSite.next();
            refBases.append(refSite.getBaseAsChar());
        }
        window.setReferenceContig(refSite.getCurrentContig());

        // Iterate forward to get all reference ordered data covering this interval
        final RefMetaDataTracker tracker = getReferenceOrderedDataAtLocus(window.getLocation());

        sum = walkAtinterval( walker, sum, window, refBases.toString(), tracker );

        printProgress("intervals", window.getLocation());

        return sum;
    }

    protected <M, T> T walkAtinterval( final LocusWindowWalker<M, T> walker,
                                       T sum,
                                       final LocusContext locus,
                                       final String refSeq,
                                       final RefMetaDataTracker tracker ) {

        //
        // Execute our contract with the walker.  Call filter, map, and reduce
        //
        final boolean keepMeP = walker.filter(tracker, refSeq, locus);
        if (keepMeP) {
            M x = walker.map(tracker, refSeq, locus);
            sum = walker.reduce(x, sum);
        }

        //printProgress("intervals", interval.getLocation());
        return sum;
    }
}