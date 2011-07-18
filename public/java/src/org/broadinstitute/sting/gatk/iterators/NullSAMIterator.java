package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;

import java.util.Iterator;
import java.util.NoSuchElementException;
/**
 * User: hanna
 * Date: May 19, 2009
 * Time: 6:47:16 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A placeholder for an iterator with no data.
 */
public class NullSAMIterator implements StingSAMIterator {
    public NullSAMIterator() {}

    public Iterator<SAMRecord> iterator() { return this; }
    public void close() { /* NO-OP */ }

    public boolean hasNext() { return false; }
    public SAMRecord next() { throw new NoSuchElementException("No next element is available."); }
    public void remove() { throw new UnsupportedOperationException("Cannot remove from a StingSAMIterator"); }
}
