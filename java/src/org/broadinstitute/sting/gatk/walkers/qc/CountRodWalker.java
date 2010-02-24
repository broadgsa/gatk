package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import java.util.LinkedList;

/**
 * Prints out counts of the number of reference ordered data objects are
 * each locus for debugging RodWalkers.
 */
public class CountRodWalker extends RodWalker<CountRodWalker.Datum, Pair<ExpandingArrayList<Long>, Long>> {
    @Argument(fullName = "verbose", shortName = "v", doc="If true, Countrod will print out detailed information about the rods it finds and locations", required = false)
    public boolean verbose = false;

    @Argument(fullName = "showSkipped", shortName = "s", doc="If true, CountRod will print out the skippped locations", required = false)
    public boolean showSkipped = false;

    public class Datum {
        public long nRodsAtThisLocation = 0;
        public long nSkippedBases =0;
        public long nTotalBases = 0;

        public Datum( long nRodsAtThisLocation, long nSkippedBases, long nTotalBases ) {
            this.nRodsAtThisLocation = nRodsAtThisLocation;
            this.nSkippedBases = nSkippedBases;
            this.nTotalBases = nTotalBases;
        }

        public String toString() {
            return String.format("<%d %d %d>", nRodsAtThisLocation, nSkippedBases, nTotalBases);
        }
    }

    public Datum map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        GenomeLoc cur = context.getLocation();

        if ( verbose && showSkipped ) {
            for ( int i = 0; i < context.getSkippedBases(); i++ ) {
                out.printf("%s: skipped%n", GenomeLocParser.incPos(cur, i - context.getSkippedBases()));
            }
        }

        long nRodsHere = 0;
        long nTotalBases = 0;

        if ( ref == null ) {
            // we're getting the last skipped update
            if ( verbose )
                out.printf("Last position was %s: skipping %d bases%n",
                        context.getLocation(), context.getSkippedBases() );
            nRodsHere = -1; // don't update this
            nTotalBases = context.getSkippedBases();
        } else {
            Collection<RODRecordList> rods = new LinkedList<RODRecordList>();
            for ( RODRecordList rod : tracker.getBoundRodTracks() ) {
                //System.out.printf("Considering rod %s%n", rod);
                if ( rod.getLocation().getStart() == context.getLocation().getStart() && ! rod.getName().equals("interval") ) {
                    // only consider the first element
                    //System.out.printf("adding it%n");
                    rods.add(rod);
                }
            }

            nRodsHere = rods.size();

            if ( nRodsHere > 0 ) {
                if ( verbose ) {
                    List<String> names = new ArrayList<String>();
                    for ( RODRecordList rod : rods ) {
                        names.add(rod.getName());
                    }

                    //System.out.printf("context is %s", context.getSkippedBases());
                    out.printf("At %s: found %d rod(s) [%s] after skipping %d bases%n",
                            context.getLocation(), nRodsHere, Utils.join(",", names), context.getSkippedBases() );
                }
            }

            nTotalBases = context.getSkippedBases() + 1;
        }

        return new Datum(nRodsHere, context.getSkippedBases(), nTotalBases);
    }

    public Pair<ExpandingArrayList<Long>, Long> reduceInit() {
        return new Pair<ExpandingArrayList<Long>, Long>(new ExpandingArrayList<Long>(), 0l);
    }

    private void updateCounts(ExpandingArrayList<Long> counts, long nRods, long nObs) {
        if ( nRods >= 0 ) {
            long prev = counts.get((int)nRods) == null ? 0l : counts.get((int)nRods);
            counts.set((int)nRods, nObs + prev);
        }
    }

    public Pair<ExpandingArrayList<Long>, Long> reduce(Datum point, Pair<ExpandingArrayList<Long>, Long> sum) {
        ExpandingArrayList<Long> counts = sum.getFirst();
        updateCounts(counts, point.nRodsAtThisLocation, 1);
        updateCounts(counts, 0, point.nSkippedBases);

        Pair<ExpandingArrayList<Long>, Long> r = new Pair<ExpandingArrayList<Long>, Long>(counts, point.nTotalBases + sum.getSecond());

        //System.out.printf("Reduce: %s %s => %s%n", point, sum, r);
        return r;
    }
}