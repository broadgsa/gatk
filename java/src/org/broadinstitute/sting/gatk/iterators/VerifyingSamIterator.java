package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.RuntimeIOException;

import java.util.Iterator;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 15, 2009
 * Time: 6:02:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class VerifyingSamIterator implements Iterator<SAMRecord> {
    Iterator<SAMRecord> it;
    SAMRecord last = null;
    boolean checkOrderP = true;
    long nOutOfOrderReads = 0; 

    public VerifyingSamIterator(Iterator<SAMRecord> it) {
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

    /**
     * If true, enables ordered checking of the reads in the file.  By default this is enabled.
     * @param checkP If true, sam records will be checked to insure they come in order
     */
    public void setCheckOrderP( boolean checkP ) {
        checkOrderP = checkP;
    }

    public void verifyRecord( final SAMRecord last, final SAMRecord cur ) {
        if ( checkOrderP && isOutOfOrder(last, cur) ) {
            this.last = null;
            throw new RuntimeIOException(String.format("Reads are out of order:%nlast:%n%s%ncurrent:%n%s%n", last.format(), cur.format()) );
        }
    }

    public static boolean isOutOfOrder( final SAMRecord last, final SAMRecord cur ) {
        if ( last == null || cur.getReadUnmappedFlag() )
            return false;
        else {
            GenomeLoc lastLoc = Utils.genomicLocationOf( last );
            GenomeLoc curLoc = Utils.genomicLocationOf( cur );
            return curLoc.compareTo(lastLoc) == -1;
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}
