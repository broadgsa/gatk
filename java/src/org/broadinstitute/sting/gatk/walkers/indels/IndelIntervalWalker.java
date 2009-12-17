package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;

import java.util.List;

/**
 * Emits intervals consisting of indels from the aligned reads.
 */
@WalkerName("IndelIntervals")
@Requires({DataSource.READS, DataSource.REFERENCE})
@ReadFilters({Platform454Filter.class, ZeroMappingQualityReadFilter.class})
public class IndelIntervalWalker extends ReadWalker<IndelIntervalWalker.Interval, IndelIntervalWalker.Interval> {
    @Argument(fullName="allow454Reads", shortName="454", doc="process 454 reads", required=false)
    boolean allow454 = false;
    @Argument(fullName="minIndelsPerInterval", shortName="minIndels", doc="min indels per interval", required=false)
    int minIntervalIndelCount = 1;

    public void initialize() {}

    public boolean filter(char[] ref, SAMRecord read) {
        return ( !read.getReadUnmappedFlag() &&            // mapped
                 read.getMappingQuality() != 0 &&          // positive mapping quality
                 read.getCigarLength() > 1 &&   // indel
                 (allow454 || !Utils.is454Read(read)) );
    }

    public Interval map(char[] ref, SAMRecord read) {
//        List<AlignmentBlock> blocks = read.getAlignmentBlocks();
//        long indelLeftEdge = read.getAlignmentStart() + blocks.get(0).getLength() - 1;
//        long indelRightEdge = read.getAlignmentEnd() - blocks.get(blocks.size()-1).getLength() + 1;

        long indelLeftEdge = -1;
        long indelRightEdge = -1;
        int lengthOnRef = 0; // length of the event on the reference (this preset value is right for insertion)
        int pos = read.getAlignmentStart();
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            switch( ce.getOperator() ) {
            case S :
            case H : break; // ignore clipping completely: alignment start is set to the first NON-clipped base
            case P : // skip padding and matching bases, advance position on the ref:
            case M : pos += ce.getLength(); break;
            case D :
                  lengthOnRef = ce.getLength(); // deletion occupies non-zero number of bases on the ref
            case I :
                    // for insertion we do not have to set lengthOnRef as it was already initialized to 0
                if ( indelLeftEdge == -1 ) { // it's the first indel in the current read
                    indelLeftEdge = pos - 1; // set to the last (matching) base before the indel
                }
                pos += lengthOnRef;
                indelRightEdge = pos; // whether it is the first indel in the read or not, we set right edge immediately after it
                break;
            default :
                throw new StingException("Unrecognized cigar element '"+ce.getOperator()+"' in read "+read.getReadName());
            }
        }

        // at this point indelLeftEdge == -1 or [indelLeftEdge,indelRightEdge] is the interval encompassing ALL indels in the current read:
        if ( indelLeftEdge == -1 ) return null;

        // we've seen some bam files coming from other centers that have reads (and indels!!) hanging over the contig ends;
        // we do not want to deal with this stuff (what is it, anyway??)
        if ( indelRightEdge > GenomeLocParser.getContigInfo(read.getReferenceName()).getSequenceLength() ) {
            System.out.println("WARNING: read " + read.getReadName()+" contains indel(s) spanning past the contig end (read ignored).");
            return null;
        }

//        if ( indelLeftEdge == 49563377 ) System.out.println("read: " +read.getReadName());
        GenomeLoc indelLoc = GenomeLocParser.createGenomeLoc(read.getReferenceIndex(), indelLeftEdge, indelRightEdge);
        GenomeLoc refLoc = GenomeLocParser.createGenomeLoc(read);

 //       if ( indelLeftEdge == 10313124 || indelLeftEdge == 10313170 ) System.out.println("read: " +read.getReadName() + " ; " + refLoc + " ; " + indelLoc);

        return new Interval(refLoc, indelLoc);
    }

    public Interval reduceInit() {
        return null;
    }

    // Note that the reduce step assumes that reductions occur in order (i.e. no hierarchical reductions),
    // although this can easily be changed if necessary.
    public Interval reduce(Interval value, Interval sum) {
        // if there is no interval to the left, then this is the first one
        if ( sum == null )
            return value;

        if ( value == null ) return sum; // nothing to do, wait for the next
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
            GenomeLoc mergedIndelLoc = GenomeLocParser.createGenomeLoc(this.indelLoc.getContigIndex(), indelLeftEdge, indelRightEdge);
            Interval mergedInterval = new Interval(this.readLoc.merge(i.readLoc), mergedIndelLoc);
            mergedInterval.indelCount = this.indelCount + i.indelCount;
            return mergedInterval;
        }

        public String toString() {
            return indelLoc.toString();
        }
    }
}
