package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;

import java.util.Iterator;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.Reads;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Nov 4, 2009
 * Time: 10:41:02 PM
 */
public class PlusOneFixIterator implements StingSAMIterator {

    /**
     * this holds the value that we use to correct our query overlapping calls with,
     * but more importantly it ties the ReadStreamPointer code to this so we don't loose
     * track of this adjustment across the code.
     */
    public static final Integer PLUS_ONE_FIX_CONSTANT = 1;

    /**
     * our interval region
     */
    private final GenomeLoc mInterval;

    /**
     * our iterator
     */
    private final StingSAMIterator mIterator;

    /**
     * our next record
     */
    private SAMRecord mNextRecord;

    /**
     * This iterator filters all reads that have their start < PLUS_ONE_FIX_CONSTANT from the end of the
     * interval:
     * <p/>
     * <p/>
     * (this won't make sense in a non-fixed width font)
     * --------------
     * interval     |
     * ------------------------
     *            |   good read|
     * --------------------------------------------
     *             | overlap < PLUS_ONE_FIX_CONSTANT, we drop|
     *             -------------------------------------------
     * <p/>
     * The problem is that we had to put a +1 on our queryOverlapping() stop position value to SAMTools
     * due to an indexing problem.  Other iterators don't know about this adjustment though, so calculations
     * like overlap fail.  This iterator eliminates them before the read is seen by other iterators.
     * <p/>
     * create a plus one fix iterator, given:
     *
     * @param inteveral the interval we're checking
     * @param iterator  the iterator to draw reads from
     */
    public PlusOneFixIterator(GenomeLoc inteveral, StingSAMIterator iterator) {
        if (inteveral == null || iterator == null)
            throw new StingException("Both parameters to PlusOneFixIterator have to != null");

        // save off the iterator and the interval
        mInterval = inteveral;
        mIterator = iterator;

        // pop off the first read
        if (mIterator.hasNext()) {
            next();
        }
    }


    @Override
    public boolean hasNext() {
        return (mNextRecord != null);
    }

    @Override
    public SAMRecord next() {
        SAMRecord ret = mNextRecord;
        while (mIterator.hasNext()) {
            mNextRecord = mIterator.next();
            if (!(mNextRecord.getAlignmentStart() > (mInterval.getStop() - PLUS_ONE_FIX_CONSTANT))) return ret;
        }
        mNextRecord = null;
        return ret;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("REMOVE on PlusOneFixIterator is not supported");   
    }

    @Override
    public Reads getSourceInfo() {
        return this.mIterator.getSourceInfo();
    }

    @Override
    public void close() {
        this.mIterator.close();
    }

    @Override
    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
