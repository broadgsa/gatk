package org.broadinstitute.sting.utils.recalibration;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 2/13/12
 */

public class BQSRSamIterator implements StingSAMIterator {
    private final StingSAMIterator it;
    private final BaseRecalibration bqsr;

    /**
     * Creates a new BQSRSamIterator and applies BQSR on the fly to incoming reads.
     *
     * @param it    The incoming SamIterator to wrap
     * @param bqsr  The object which holds the BQSR table information and knows how to apply it
     */
    @Requires({
            "it != null",
            "bqsr != null"})
    public BQSRSamIterator(StingSAMIterator it, BaseRecalibration bqsr) {
        if ( bqsr == null ) throw new ReviewedStingException("BUG: shouldn't create BQSRSamIterator with null recalibration object");

        this.it = it;
        this.bqsr = bqsr;
    }

    @Requires("hasNext()")
    @Ensures("result != null")
    public SAMRecord next()     {
        SAMRecord read = it.next();
        bqsr.recalibrateRead((GATKSAMRecord) read);
        return read;
    }

    public boolean hasNext()    { return this.it.hasNext(); }
    public void remove()        { throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!"); }
    public void close()         { it.close(); }
    public Iterator<SAMRecord> iterator() { return this; }
}
