package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.AlignmentBlock;

import java.util.List;

// A walker to determine intervals within which reads should be cleaned and realigned
// because they contain one or more indels (which possibly caused them not to align perfectly).
// Note that the reduce step assumes that reductions occur in order (i.e. no hierarchical reductions),
// although this can easily be changed if necessary.

@WalkerName("IndelIntervals")
public class IndelIntervalWalker extends ReadWalker<IndelIntervalWalker.Interval, IndelIntervalWalker.Interval> {
    @Argument(fullName="maxReadLength", shortName="maxRead", required=false, defaultValue="-1")
    public int maxReadLength;
    @Argument(fullName="minIndelsPerInterval", shortName="minIndels", required=false, defaultValue="1")
    public int minIntervalIndelCount;

    public void initialize() {}

    public boolean filter(LocusContext context, SAMRecord read) {
        return ( !read.getReadUnmappedFlag() &&            // mapped
                 read.getReadLength() <= maxReadLength &&  // not too big
                 read.getAlignmentBlocks().size() > 1);    // indel
    }

    public Interval map(LocusContext context, SAMRecord read) {
        List<AlignmentBlock> blocks = read.getAlignmentBlocks();
        long indelLeftEdge = read.getAlignmentStart() + blocks.get(0).getLength() - 1;
        long indelRightEdge = read.getAlignmentEnd() - blocks.get(blocks.size()-1).getLength() + 1;
        GenomeLoc indelLoc = new GenomeLoc(context.getLocation().getContigIndex(), indelLeftEdge, indelRightEdge);
        return new Interval(context.getLocation(), indelLoc);
    }

    public Interval reduceInit() {
        return null;
    }

    public Interval reduce(Interval value, Interval sum) {
        // if there is no interval to the left, then this is the first one
        if ( sum == null )
            return value;

        //System.out.println(sum + " ==> " + value);

        // if the intervals don't overlap, print out the leftmost one and start a new one
        if ( !sum.readLoc.overlapsP(value.readLoc) ) {
            if ( sum.indelCount >= minIntervalIndelCount )
                out.println(sum);
            return value;
        }
        // otherwise, merge them
        return sum.merge(value);
    }

    public void onTraversalDone(Interval result) {
        if ( result != null && result.indelCount >= minIntervalIndelCount )
            out.println(result);        
    }

    public class Interval {
        public GenomeLoc readLoc = null;
        public GenomeLoc indelLoc = null;
        public int indelCount;

        Interval(GenomeLoc readLoc, GenomeLoc indelLoc) {
            this.indelLoc = indelLoc;
            this.readLoc = readLoc;
            indelCount = 1;
        }

        public Interval merge(Interval i) {
            long indelLeftEdge = Math.min(this.indelLoc.getStart(), i.indelLoc.getStart());
            long indelRightEdge = Math.max(this.indelLoc.getStop(), i.indelLoc.getStop());
            GenomeLoc mergedIndelLoc = new GenomeLoc(this.indelLoc.getContigIndex(), indelLeftEdge, indelRightEdge);
            Interval mergedInterval = new Interval(this.readLoc.merge(i.readLoc), mergedIndelLoc);
            mergedInterval.indelCount = this.indelCount + i.indelCount;
            return mergedInterval;
        }

        public String toString() {
            return indelLoc.toString();
        }
    }
}