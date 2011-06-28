package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.Iterator;

/**
 * Verifies that the incoming stream of reads is correctly sorted
 */
public class VerifyingSamIterator implements StingSAMIterator {
    private GenomeLocParser genomeLocParser;
    StingSAMIterator it;
    SAMRecord last = null;
    boolean checkOrderP = true;

    public VerifyingSamIterator(GenomeLocParser genomeLocParser,StingSAMIterator it) {
        this.genomeLocParser = genomeLocParser;
        this.it = it;
    }

    public boolean hasNext() { return this.it.hasNext(); }
    public SAMRecord next() {

        SAMRecord cur = it.next();
        if ( last != null )
            verifyRecord(last, cur);
        if ( ! cur.getReadUnmappedFlag() )
            last = cur;
        return cur;
    }

    private void verifyRecord( final SAMRecord last, final SAMRecord cur ) {
        if ( checkOrderP && isOutOfOrder(last, cur) ) {
            this.last = null;
            throw new UserException.MissortedBAM(String.format("reads are out of order:%nlast:%n%s%ncurrent:%n%s%n", last.format(), cur.format()) );
        }
    }

    private boolean isOutOfOrder( final SAMRecord last, final SAMRecord cur ) {
        if ( last == null || cur.getReadUnmappedFlag() )
            return false;
        else {
            if(last.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX || last.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START)
                throw new UserException.MalformedBAM(last,String.format("read %s has inconsistent mapping information.",last.format()));
            if(cur.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX || cur.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START)
                throw new UserException.MalformedBAM(last,String.format("read %s has inconsistent mapping information.",cur.format()));

            GenomeLoc lastLoc = genomeLocParser.createGenomeLoc( last );
            GenomeLoc curLoc = genomeLocParser.createGenomeLoc( cur );
            return curLoc.compareTo(lastLoc) == -1;
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    public void close() {
        it.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
