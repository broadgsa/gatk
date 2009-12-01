package org.broadinstitute.sting.playground.indels;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 19, 2009
 * Time: 7:28:40 PM
 * To change this template use File | Settings | File Templates.
 */

/** This interface abstracts processing of SAM records. Its only receive() method should be called to send a record
 * to the implementation.
 */
public interface RecordReceiver {
    public void receive(SAMRecord r);
}
