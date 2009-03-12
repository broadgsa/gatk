package org.broadinstitute.sting.atk;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.PushbackIterator;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Predicate;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public abstract class LocusIterator implements Iterable<LocusContext>, CloseableIterator<LocusContext> {
    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public Iterator<LocusContext> iterator() {
        return this;
    }

    public void close() {
        //this.it.close();
    }

    public abstract boolean hasNext();
    public abstract LocusContext next();

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }
}
