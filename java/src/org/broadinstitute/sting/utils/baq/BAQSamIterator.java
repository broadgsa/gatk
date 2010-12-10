package org.broadinstitute.sting.utils.baq;

import net.sf.samtools.SAMRecord;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;

import java.util.Iterator;

/**
 * Simple iterator that applies Heng's BAQ calculation to a stream of incoming reads
 */
public class BAQSamIterator implements StingSAMIterator {
    StingSAMIterator it;
    BAQ baqHMM = new BAQ();         // creates a BAQ creator with default parameters
    IndexedFastaSequenceFile refReader = null;
    BAQ.CalculationMode cmode;
    BAQ.QualityMode qmode;

    /**
     * Creates a new BAMSamIterator using the reference getter refReader and applies the BAM to the reads coming
     * in from it.  See BAQ docs for baqType information.
     *
     * @param refReader
     * @param it
     * @param cmode
     * @param qmode
     */
    public BAQSamIterator(IndexedFastaSequenceFile refReader, StingSAMIterator it, BAQ.CalculationMode cmode, BAQ.QualityMode qmode) {
        if ( cmode == BAQ.CalculationMode.OFF) throw new ReviewedStingException("BUG: shouldn't create BAQSamIterator with calculation mode OFF");
        if ( qmode == BAQ.QualityMode.DONT_MODIFY ) throw new ReviewedStingException("BUG: shouldn't create BAQSamIterator with quailty mode DONT_MODIFY");

        this.refReader = refReader;
        this.it = it;
        this.cmode = cmode;
        this.qmode = qmode;
    }

    public SAMRecord next()     {
        //System.out.printf("BAQing during input%n");
        SAMRecord read = it.next();
        baqHMM.baqRead(read, refReader, cmode, qmode);
        return read;
    }

    public boolean hasNext()    { return this.it.hasNext(); }
    public void remove()        { throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!"); }
    public void close()         { it.close(); }
    public Iterator<SAMRecord> iterator() { return this; }
}
