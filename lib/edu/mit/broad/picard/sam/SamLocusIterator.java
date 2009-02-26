package edu.mit.broad.picard.sam;

import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.NotPrimarySkippingIterator;
import edu.mit.broad.picard.directed.GenomeMask;

import java.util.*;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public class SamLocusIterator implements Iterable<SamLocusIterator.LocusInfo>, CloseableIterator<SamLocusIterator.LocusInfo> {

    /**
     * The unit of iteration.  Holds the locus, plus the base, quality and strand for each read at that locus.
     */
    public static class LocusInfo {
        protected final int sequenceIndex;
        protected final int position;
        protected final List<Byte> bases = new ArrayList<Byte>(100);
        protected final List<Byte> qualities = new ArrayList<Byte>(100);
        protected final List<Boolean> negativeStrandFlags = new ArrayList<Boolean>(100);

        LocusInfo(final int sequenceIndex, final int position) {
            this.sequenceIndex = sequenceIndex;
            this.position = position;
        }

        /**
         * Accumulate info for one read at the locus.
         */
        public void add(final Byte readBase, final Byte baseQuality, final boolean strand) {
            bases.add(readBase);
            qualities.add(baseQuality);
            negativeStrandFlags.add(strand);
        }

        public int getSequenceIndex() { return sequenceIndex; }
        public int getPosition() { return position; }
        public List<Byte> getBases() { return bases; }
        public List<Byte> getQualities() { return qualities; }
        public List<Boolean> getNegativeStrandFlags() { return negativeStrandFlags; }

        public String getBasesAsString() { return bytesToString(bases); }

        private static String bytesToString(final List<Byte> data) {
            if (data == null || data.size() == 0) {
                return "";
            }

            final char[] chars = new char[data.size()];
            for (int i = 0; i < data.size(); i++) {
                chars[i] = (char) (data.get(i) & 0xFF);
            }
            return new String(chars);
        }
    }




    private final CloseableIterator<SAMRecord> underlyingIterator;
    private final NotPrimarySkippingIterator it;
    private final LinkedList<LocusInfo> complete = new LinkedList<LocusInfo>();
    private final LinkedList<LocusInfo> accumulator = new LinkedList<LocusInfo>();

    private boolean includeNonPfReads = false;
    private boolean includeDuplicates = false;
    private int qualityScoreCutoff = -Integer.MAX_VALUE;
    
    private GenomeMask mask;
    private int lastContig = 0;
    private int lastPosition = 0;

    private boolean finishedAlignedReads = false;


    // this should probably take a SAM
    public SamLocusIterator(final CloseableIterator<SAMRecord> samIterator) {
        this.underlyingIterator = samIterator;
        this.it = new NotPrimarySkippingIterator(samIterator);
    }

    public Iterator<LocusInfo> iterator() {
        return this;
    }

    public void close() {
        this.underlyingIterator.close();
    }

    private boolean samHasMore() {
        return !finishedAlignedReads && it.hasCurrent();
    }
    public boolean hasNext() {
        return ((complete.size() > 0) || (accumulator.size() > 0) || (samHasMore()) || hasRemainingMaskBases());
    }

    private boolean hasRemainingMaskBases() {
        if (mask == null) return false;

        // if there are more contigs in the mask, by definition some of them must have
        // marked bases otherwise if we're in the last contig, but we're not at the last marked position,
        // there is also more in the mask
        return (lastContig <= mask.getMaxContig() ||
               (lastContig == mask.getMaxContig() && lastPosition <= mask.get(lastContig).nextSetBit(lastPosition+1)));
    }

    public LocusInfo next() {

        // if we don't have any completed entries to return, try and make some!
        while(complete.size() == 0 && samHasMore()) {
            final SAMRecord rec = it.getCurrent();
            final String cigar = rec.getCigarString();
            
            // as soon as we hit our first non-aligned read, we can stop!
            if (cigar.equals("*")) {
                this.finishedAlignedReads = true;
                continue;
            }

            // skip dupe reads, if so requested
            if (!isIncludeDuplicates() && rec.getDuplicateReadFlag()) { it.advance(); continue; }

            // skip non-PF reads, if so requested
            if (!isIncludeNonPfReads() && rec.getReadFailsVendorQualityCheckFlag()) { it.advance(); continue; }
            
            // when we switch contigs, emit everything in the accumulator
            if (accumulator.size() > 0 && accumulator.getFirst().sequenceIndex != rec.getReferenceIndex()) {
                while (accumulator.size() > 0) {
                    popLocus();
                }
            }

            // pop off things we're not going to accumulate more coverage at the locus in question
            while(accumulator.size() > 0 && accumulator.getFirst().position < rec.getAlignmentStart()) {
                popLocus();
            }

            // check that it's a non-gapped alignment for now!
            // TODO: handle gapped and clipped alignments
            if (!cigar.matches("[0-9]+M")) {
                System.out.println("Cannot deal with clipped or gapped alignments. CIGAR="+cigar);
                System.exit(1);
            }

            // at this point, either the list is empty or the head should
            // be the same position as the first base of the read

            // interpret the CIGAR string and add the base info
            for(int j=0; j < rec.getReadBases().length; j++) {
                // if the position is empty, initialize it
                if (j > accumulator.size() - 1) {
                    accumulator.add(new LocusInfo(rec.getReferenceIndex(), rec.getAlignmentStart() + j));
                }

                // if the quality score cutoff is met, accumulate the base info
                if (rec.getBaseQualities()[j] >= getQualityScoreCutoff()) {
                    accumulator.get(j).add(rec.getReadBases()[j], rec.getBaseQualities()[j], rec.getReadNegativeStrandFlag());
                }
            }


            it.advance();
        }

        // if we have nothing to return to the user, and we're at the end of the SAM iterator,
        // push everything into the complete queue
        if (complete.size() == 0 && !samHasMore()) {
            while(accumulator.size() > 0) {
                popLocus();
            }
        }

        // if there are completed entries, return those
        if (complete.size() > 0) {
            return complete.removeFirst();
        } else {

            // In this case... we're past the last read from SAM so see if we can
            // fill out any more (zero coverage) entries from the mask
            LocusInfo zeroResult = null;
            while (zeroResult == null && lastContig <= mask.getMaxContig()) {
                final int nextbit = mask.get(lastContig).nextSetBit(lastPosition+1);

                // try the next contig
                if (nextbit == -1) {
                    lastContig++;
                    lastPosition = 0;
                } else {
                    lastPosition = nextbit;
                    zeroResult = new LocusInfo(lastContig, lastPosition);
                }
            }

            return zeroResult;
        }
    }

    /**
     * Pop the first entry from the LocusInfo accumulator into the complete queue.  In addition,
     * check the GenomeMask and if there are intervening mask positions between the last popped base and the one
     * about to be popped, put those on the complete queue as well.
     */
    private void popLocus() {
        final LocusInfo li = accumulator.removeFirst();

        // fill in any gaps based on our genome mask
        final int liContig = li.getSequenceIndex();

        // if we're not on the same contig, fill in the rest of the bits for the previous contig first...
        if (lastContig < liContig) {
            while (lastContig < liContig) {
                int nextbit = 0;

                if (mask != null && mask.get(lastContig) != null) {
                    while (nextbit != -1) {
                        nextbit = mask.get(lastContig).nextSetBit(lastPosition + 1);
                        if (nextbit > -1) {
                            complete.addLast(new LocusInfo(lastContig, nextbit));
                            lastPosition = nextbit;
                        }
                    }
                }
                lastPosition=0;
                lastContig++;
            }
        }

        // now that we're on the same contig, fill in any unfilled positions
        // if we have some bits in the mask to fill in...
        if (mask != null && mask.get(lastContig) != null && lastPosition + 1 < li.getPosition()) {
            while (lastPosition + 1 < li.getPosition()) {

                final int nextbit = mask.get(lastContig).nextSetBit(lastPosition + 1);

                // if there are no more mask bits, or the next mask bit is
                // at or after the current data, just continue on
                if (nextbit == -1 || nextbit >= li.getPosition()) { break; }

                // otherwise, pop on the desired empty locus info
                complete.addLast(new LocusInfo(lastContig, nextbit));
                lastPosition = nextbit;
            }
        }

        // only add to the complete queue if it's in the mask (or we have no mask!)
        if (mask == null || mask.get(li.getSequenceIndex(), li.getPosition())) {
            complete.addLast(li);
        }

        lastContig = liContig;
        lastPosition = li.getPosition();


    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    // --------------------------------------------------------------------------------------------
    // Helper methods below this point...
    // --------------------------------------------------------------------------------------------

    public void setGenomeMask(final GenomeMask mask) { this.mask = mask; }
    public GenomeMask getGenomeMask() { return this.mask; }

    public boolean isIncludeNonPfReads() { return includeNonPfReads; }
    public void setIncludeNonPfReads(final boolean includeNonPfReads) { this.includeNonPfReads = includeNonPfReads; }

    public boolean isIncludeDuplicates() { return includeDuplicates; }
    public void setIncludeDuplicates(final boolean includeDuplicates) { this.includeDuplicates = includeDuplicates; }

    public int getQualityScoreCutoff() { return qualityScoreCutoff; }
    public void setQualityScoreCutoff(final int qualityScoreCutoff) { this.qualityScoreCutoff = qualityScoreCutoff; }


}
