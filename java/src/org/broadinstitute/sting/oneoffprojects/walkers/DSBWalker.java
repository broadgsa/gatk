package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Dec 3, 2009
 * Time: 11:54:35 AM
 * To change this template use File | Settings | File Templates.
 */
@Requires({DataSource.READS, DataSource.REFERENCE})

public class DSBWalker extends LocusWalker<Integer,Integer> {
    @Argument(fullName="coverage",shortName="C",doc="Regions with coverage above specified threshold will be reported",required=true)
    int COV_CUTOFF = 0;
    @Argument(fullName="minLength",shortName="ml",doc="Only regions longer than the specified value will be reported",required=false)
    int MINLENGTH_CUTOFF = 0;

    private int MERGE_DIST = 300; // merge intervals that are closer than this distance from one another

    private long maxcov = 0;
    private long maxz = 0;
    private long mergedmaxcov = 0;
    private long mergedmaxz = 0;
    GenomeLoc mergedInterval = null;
    GenomeLoc currentInterval = null;

    private long nIntervals = 0;

    private void emit(GenomeLoc l) {
        if ( mergedInterval == null ) {
            mergedInterval = l.clone();
            mergedmaxcov = maxcov;
            mergedmaxz = maxz;
            return;
        }

        if ( mergedInterval.getContigIndex() != l.getContigIndex() ) {
            long length = mergedInterval.getStop()-mergedInterval.getStart()+1;
            if ( length >= MINLENGTH_CUTOFF ) {
                out.println(mergedInterval+"\t"+length+"\t"+mergedmaxcov+"\t"+mergedmaxz); // eject old interval
                nIntervals++;
            }
            mergedInterval = l.clone();
            mergedmaxcov = maxcov;
            mergedmaxz = maxz;
            return;
        }

        // merged interval exists and new interval is on the same contig. Check if the new interval
        // is close enough so we got to merge and keep waiting:

        if ( l.getStart() - mergedInterval.getStop() < MERGE_DIST ) {
            mergedInterval = GenomeLocParser.setStop(mergedInterval,l.getStop());
            if ( maxcov > mergedmaxcov) mergedmaxcov = maxcov;
            if ( maxz > mergedmaxz ) mergedmaxz = maxz;
            return;
        }

        // nope, new interval is far enough. Print old one and keep current one.

        long length = mergedInterval.getStop()-mergedInterval.getStart()+1;
        if ( length >= MINLENGTH_CUTOFF ) {
            out.println(mergedInterval+"\t"+length+"\t"+mergedmaxcov+"\t"+mergedmaxz); // eject old interval
            nIntervals++;
        }
        mergedInterval = l.clone();
        mergedmaxcov = maxcov;
        mergedmaxz = maxz;

    }

    public void onTraversalDone() {
        if ( mergedInterval != null ) {
            long length = mergedInterval.getStop()-mergedInterval.getStart()+1;
            if ( length >= MINLENGTH_CUTOFF ) {
                out.println(mergedInterval+"\t"+length+"\t"+mergedmaxcov+"\t"+mergedmaxz); // eject old interval
                nIntervals++;
            }
        }
        System.out.println(nIntervals+" intervals detected.");
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        ReadBackedPileup pileup = context.getPileup();
        List<SAMRecord> reads = pileup.getReads();

        int nZero = pileup.getNumberOfMappingQualityZeroReads();

        int nonZCoverage = reads.size() - nZero;

        if ( nonZCoverage >= COV_CUTOFF ) {

            // if we were not inside an interval, start one:
            if ( currentInterval == null ) {
                maxcov = nonZCoverage;
                maxz = nZero;
                currentInterval = context.getLocation().clone();
//                System.out.println("Setting current to "+currentInterval);
                return 0;
            }

            // if we were inside an interval and we just jumped onto a new contig, get rid of the old interval
            if ( currentInterval.compareContigs(context.getLocation()) != 0 ) {
                // we just moved to a new contig
                System.out.println("On contig "+context.getLocation().getContig());
                emit(currentInterval);
                maxcov = nonZCoverage;
                maxz = nZero;
                currentInterval = context.getLocation().clone();
                return 0;
            }

            // we are on the same contig, we are within the interval, so we need to extend the current interval:
            currentInterval = GenomeLocParser.setStop(currentInterval,context.getLocation().getStop()); // still within the interval, adjust stop
            //System.out.println("Extending current to "+currentInterval +" ("+context.getLocation()+", "+context.getLocation().getStop()+")");
            if ( nonZCoverage > maxcov ) maxcov = nonZCoverage; // adjust maxcov
            if ( nZero > maxz ) maxz = nZero; // adjust maxz
        } else {
            // low coverage, if we were inside an interval, it stops now:
            if ( currentInterval != null ) {
 //              System.out.println("Emitting current as "+currentInterval);
               emit(currentInterval);
               currentInterval = null;
               maxcov = 0;
               maxz = 0;
            }
        }

        return 0;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum+value;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
