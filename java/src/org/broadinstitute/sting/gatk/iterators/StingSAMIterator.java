package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
/**
 *
 * User: aaron
 * Date: May 6, 2009
 * Time: 5:30:41 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */

/**
 * @author aaron
 * @version 1.0
 * @date May 6, 2009
 * <p/>
 * Interface ClosableGetHeaderIterator
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */
public interface StingSAMIterator extends CloseableIterator<SAMRecord>, Iterable<SAMRecord> {

    /**
     * gets the header from the iterator
     * @return the samfileheader for the iterator, null if one is not available
     */
    public SAMFileHeader getHeader();
}
