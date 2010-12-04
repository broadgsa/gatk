package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.RuntimeIOException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

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
            throw new RuntimeIOException(String.format("Reads are out of order:%nlast:%n%s%ncurrent:%n%s%n", last.format(), cur.format()) );
        }
    }

    private boolean isOutOfOrder( final SAMRecord last, final SAMRecord cur ) {
        if ( last == null || cur.getReadUnmappedFlag() )
            return false;
        else {
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
