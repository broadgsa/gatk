package org.broadinstitute.sting.playground.indels;

import net.sf.samtools.SAMRecord;

import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 19, 2009
 * Time: 7:51:47 PM
 * To change this template use File | Settings | File Templates.
 */

/** This interface abstracts processing of piles (collections) of SAM records.
 * Its only receive() method should be called to send a collection of records
 * to the implementation.
 */

public interface RecordPileReceiver {
    public void receive(Collection<SAMRecord> c) ;
}
